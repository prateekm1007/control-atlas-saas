"""
Diagnostic calibration runner with evidence-bound acquisition and output.

The runner measures IN-PROCESS through the same engine functions used by the
API path. It performs no authenticated API calls, reads no secrets, and never
enables open mode.

Boundaries enforced here:

* every coordinate file is hashed, and a supplied pin binds the exact bytes;
* local files must match a supplied prior digest when one is provided;
* AlphaFold acquisition is metadata-first and entity-exact; authoritative
  inputs require an exact entity/version/URL/format/size/digest pin and are
  never promoted to a newer runtime version; unpinned inputs are marked
  non-authoritative diagnostics;
* every redirect hop is validated before it is requested;
* governance registry and law-semantic locks are verified before output;
* existing run directories are never overwritten;
* a dirty engine tree aborts the run before a directory is created;
* incomplete acquisition or measurement produces a non-zero exit code.

This module does not choose thresholds, tune cohorts, or evaluate a held-out
split. It writes diagnostic evidence only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin, urlparse

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import tos  # noqa: E402
from tos.engine.adjudicator import AdjudicationInput, adjudicate_laws  # noqa: E402
from tos.engine.tier1_measurements import Tier1Measurements  # noqa: E402
from tos.ingestion.processor import IngestionError, IngestionProcessor  # noqa: E402
from tos.governance.station_sop import (  # noqa: E402
    DETERMINISTIC_COUNT,
    LAW_CANON,
    LAW_CANON_HASH,
    LAW_METHOD_CLASSIFICATIONS,
    STATION_METADATA,
)
from tos.governance.modality_matrix import compute_matrix_hash  # noqa: E402

RUNNER_VERSION = "v25-runner/1.2.0"
RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
HTTP_TIMEOUT_S = 60
MAX_REDIRECTS = 5
MAX_STRUCTURE_BYTES = 50 * 1024 * 1024
MIN_STRUCTURE_BYTES = 200
USER_AGENT = "toscanini-calibration-runner/1.2 (diagnostic evidence run)"
ALLOWED_HOSTS = {"files.rcsb.org", "alphafold.ebi.ac.uk"}
ALLOWED_PORTS = {None, 443}
REDIRECT_STATUSES = {301, 302, 303, 307, 308}
HEX64 = re.compile(r"^[0-9a-f]{64}$")
PDB_ID = re.compile(r"^[0-9A-Za-z]{4}$")


class CalibrationRunError(RuntimeError):
    """A calibration run precondition or evidence boundary failed."""


class AcquisitionError(RuntimeError):
    """A coordinate acquisition failed with a recordable reason."""


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256(Path(path).read_bytes())


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _git_output(*args: str) -> str | None:
    try:
        engine_repo = Path(tos.__file__).resolve().parents[2]
        return subprocess.run(
            ["git", "-C", str(engine_repo), *args],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except Exception:
        return None


def _engine_provenance() -> dict:
    code_sha = _git_output("rev-parse", "HEAD")
    if code_sha is None:
        raise CalibrationRunError(
            "calibration runs require a git checkout of the engine "
            "(could not resolve the imported tos package repository)"
        )
    status = _git_output("status", "--porcelain") or ""
    return {
        "code_sha": code_sha,
        "code_sha_source": "git rev-parse HEAD of the imported tos package repository",
        "working_tree_clean": not status,
        "tree_status": status.splitlines()[:50],
    }


def _law_semantic_hash() -> str:
    payload = {
        "laws": {
            law_id: {
                "title": data["title"],
                "threshold": data["threshold"],
                "unit": data["unit"],
                "operator": data["operator"],
                "scope": data["scope"],
                "type": data["type"],
            }
            for law_id, data in sorted(LAW_CANON.items())
        },
        "method_classifications": dict(sorted(LAW_METHOD_CLASSIFICATIONS.items())),
    }
    return _sha256(json.dumps(payload, sort_keys=True).encode("utf-8"))


def _verify_policy_locks(registry_path: Path) -> dict:
    governance = registry_path.parent
    registry_lock_path = governance / "law_registry.sha256"
    semantic_lock_path = governance / "law_semantic.sha256"
    _require_file(registry_path, "law registry")
    _require_file(registry_lock_path, "law registry lock")
    _require_file(semantic_lock_path, "law semantic lock")
    registry_lock = registry_lock_path.read_text(encoding="utf-8").strip().lower()
    semantic_lock = semantic_lock_path.read_text(encoding="utf-8").strip().lower()
    if not HEX64.fullmatch(registry_lock):
        raise CalibrationRunError("law_registry.sha256 is not a 64-hex digest")
    if not HEX64.fullmatch(semantic_lock):
        raise CalibrationRunError("law_semantic.sha256 is not a 64-hex digest")
    actual_registry = _sha256_file(registry_path)
    actual_semantic = _law_semantic_hash()
    if actual_registry != registry_lock:
        raise CalibrationRunError(
            "law registry digest does not match governance/law_registry.sha256")
    if actual_semantic != semantic_lock:
        raise CalibrationRunError(
            "law semantic hash does not match governance/law_semantic.sha256")
    return {
        "law_registry_sha256": actual_registry,
        "law_semantic_sha256": actual_semantic,
        "registry_lock_verified": True,
        "semantic_lock_verified": True,
    }


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise CalibrationRunError(f"{label} is missing: {path}")


def _assert_https(url: str) -> None:
    parsed = urlparse(url)
    if parsed.scheme != "https":
        raise AcquisitionError(f"non-HTTPS URL refused: {url}")
    if parsed.username is not None or parsed.password is not None:
        raise AcquisitionError(f"credentials in URL refused: {url}")
    try:
        port = parsed.port
    except ValueError as exc:
        raise AcquisitionError(f"invalid port in URL: {url}") from exc
    if port not in ALLOWED_PORTS:
        raise AcquisitionError(f"port not allowed: {port}")
    if parsed.hostname not in ALLOWED_HOSTS:
        raise AcquisitionError(f"host not allowed: {parsed.hostname}")


def _http_get(url: str, timeout: int = HTTP_TIMEOUT_S):
    return requests.get(
        url, timeout=timeout, headers={"User-Agent": USER_AGENT},
        allow_redirects=False,
    )


def _redirect_target(response, current_url: str) -> str | None:
    status = getattr(response, "status_code", None)
    if status not in REDIRECT_STATUSES:
        return None
    headers = getattr(response, "headers", None) or {}
    location = headers.get("Location") or headers.get("location")
    if not location:
        raise AcquisitionError(f"HTTP {status} redirect without Location")
    return urljoin(current_url, location)


def _download(url: str, fetcher=None, min_bytes: int = MIN_STRUCTURE_BYTES) -> bytes:
    get = fetcher or _http_get
    current_url = url
    for hop in range(MAX_REDIRECTS + 1):
        _assert_https(current_url)
        try:
            response = get(current_url, timeout=HTTP_TIMEOUT_S)
        except Exception as exc:
            raise AcquisitionError(
                f"request failed: {type(exc).__name__}: {exc}") from exc
        redirect = _redirect_target(response, current_url)
        if redirect is not None:
            if hop == MAX_REDIRECTS:
                raise AcquisitionError(
                    f"too many redirects (>{MAX_REDIRECTS}) for {url}")
            _assert_https(redirect)
            current_url = redirect
            continue

        status = getattr(response, "status_code", None)
        if status != 200:
            raise AcquisitionError(f"HTTP {status} from {current_url}")
        final_url = getattr(response, "url", current_url) or current_url
        _assert_https(final_url)
        if urlparse(final_url).scheme != "https":
            raise AcquisitionError(f"non-HTTPS final URL refused: {final_url}")
        content = getattr(response, "content", b"")
        if not isinstance(content, (bytes, bytearray)):
            raise AcquisitionError("response body is not bytes")
        content = bytes(content)
        if len(content) > MAX_STRUCTURE_BYTES:
            raise AcquisitionError(f"structure exceeds {MAX_STRUCTURE_BYTES} bytes")
        if len(content) < min_bytes:
            raise AcquisitionError(f"payload is implausibly small ({len(content)} bytes)")
        return content
    raise AcquisitionError(f"too many redirects (>{MAX_REDIRECTS}) for {url}")


def detect_structure_format(content: bytes) -> str:
    """Detect PDB vs mmCIF from content. JSON is never a structure."""
    if not content or not content.strip():
        raise ValueError("empty content: not a structure")
    head = content.lstrip()[:512]
    if head.startswith(b"{") or head.startswith(b"["):
        raise ValueError("JSON payload is not a structure")
    if head.startswith(b"data_") or head.startswith(b"#"):
        return "cif"
    if b"_atom_site." in content[:16384] or b"loop_" in head:
        return "cif"
    first = content.lstrip().splitlines()[0][:6]
    if first in (b"ATOM  ", b"HETATM", b"HEADER", b"TITLE ", b"REMARK",
                 b"COMPND", b"CRYST1", b"SEQRES", b"MODEL ", b"END   "):
        return "pdb"
    raise ValueError("unrecognized structure format")


def normalize_accession(candidate_id: str) -> str:
    clean = candidate_id.strip().upper()
    if clean.startswith("AF-"):
        clean = clean[3:]
    if "-F" in clean:
        head, _, tail = clean.rpartition("-F")
        if tail.isdigit():
            clean = head
    return clean


def _afdb_entity_id(candidate_id: str) -> str:
    clean = candidate_id.strip().upper()
    if clean.startswith("AF-"):
        match = re.fullmatch(r"AF-([A-Z0-9]+)-F(\d+)", clean)
        if match:
            return f"AF-{match.group(1)}-F{int(match.group(2))}"
        return f"AF-{normalize_accession(clean)}-F1"
    return f"AF-{normalize_accession(clean)}-F1"


def _afdb_record(pdb_id: str) -> dict:
    accession = normalize_accession(pdb_id)
    entity = _afdb_entity_id(pdb_id)
    return {
        "pdb_id": pdb_id,
        "acquisition_status": "failed",
        "authoritative": False,
        "pinned": False,
        "source": None,
        "sha256": None,
        "bytes": None,
        "fetched_at_utc": None,
        "detail": None,
        "api_url": AFDB_API.format(accession=accession),
        "requested_entity_id": entity,
    }


def _validate_afdb_expectation(record: dict, expectation: dict | None) -> bool:
    if expectation is None:
        record["diagnostic_only"] = True
        return True
    required = ("entity_id", "version", "source_url", "format", "sha256", "bytes")
    missing = [field for field in required if field not in expectation]
    if missing:
        record["detail"] = (
            f"incomplete AlphaFold pin; missing {missing}"
        )
        return False
    if expectation.get("entity_id") != record["requested_entity_id"]:
        record["detail"] = (
            f"entity pin {expectation.get('entity_id')!r} does not match "
            f"requested entity {record['requested_entity_id']!r}"
        )
        return False
    if not isinstance(expectation.get("version"), int) or isinstance(
            expectation.get("version"), bool) or expectation["version"] < 1:
        record["detail"] = "AlphaFold version pin must be a positive integer"
        return False
    if expectation.get("format") != "pdb":
        record["detail"] = "AlphaFold pin must bind the PDB coordinate format"
        return False
    if not isinstance(expectation.get("bytes"), int) or expectation["bytes"] < 0:
        record["detail"] = "AlphaFold byte pin must be a non-negative integer"
        return False
    if not isinstance(expectation.get("source_url"), str):
        record["detail"] = "AlphaFold source_url pin must be a string"
        return False
    pinned_digest = expectation.get("sha256")
    if not isinstance(pinned_digest, str) or not HEX64.fullmatch(pinned_digest.lower()):
        record["detail"] = "AlphaFold digest pin must be 64 hex characters"
        return False
    try:
        _assert_https(expectation["source_url"])
    except AcquisitionError as exc:
        record["detail"] = f"invalid AlphaFold source_url pin: {exc}"
        return False
    record["pinned"] = True
    record["pinned_version"] = expectation["version"]
    record["pinned_source_url"] = expectation["source_url"]
    return True


def _acquire_afdb(pdb_id: str, expectation: dict | None = None,
                  fetcher=None) -> tuple[dict, bytes | None]:
    record = _afdb_record(pdb_id)
    if not _validate_afdb_expectation(record, expectation):
        return record, None
    accession = normalize_accession(pdb_id)
    entity = record["requested_entity_id"]
    try:
        raw = _download(record["api_url"], fetcher=fetcher, min_bytes=2)
        record["metadata_sha256"] = _sha256(raw)
        entries = json.loads(raw.decode("utf-8"))
    except (AcquisitionError, ValueError, UnicodeDecodeError) as exc:
        record["detail"] = f"metadata acquisition failed: {exc}"[:300]
        if "HTTP 404" in str(exc):
            record["detail"] = "no AlphaFold prediction returned (HTTP 404)"
        return record, None

    if not isinstance(entries, list) or not entries:
        record["detail"] = "no predictions returned"
        return record, None
    entries = [e for e in entries if isinstance(e, dict)]
    matches = [e for e in entries if e.get("modelEntityId") == entity]
    if not matches:
        served = [e.get("modelEntityId") for e in entries]
        record["served_entity_ids"] = served
        record["migrated_identifier"] = served[0] if served else None
        record["detail"] = (
            f"requested entity {entity} not served; API returned {served}"
        )[:300]
        return record, None

    try:
        versions = {
            entry["modelEntityId"]: int(entry.get("latestVersion", 0))
            for entry in matches
        }
    except (TypeError, ValueError):
        record["detail"] = "latestVersion is not an integer"
        return record, None

    if record["pinned"]:
        pinned_version = record["pinned_version"]
        version_matches = [
            entry for entry in matches
            if int(entry.get("latestVersion", 0)) == pinned_version
        ]
        if not version_matches:
            record["served_versions"] = versions
            record["detail"] = (
                f"pinned version {pinned_version} for {entity} is not served; "
                "no runtime version promotion performed"
            )
            return record, None
        if len(version_matches) != 1:
            record["detail"] = "multiple API entries claim the pinned entity/version"
            return record, None
        best = version_matches[0]
        latest_version = pinned_version
    else:
        try:
            best = max(matches, key=lambda e: int(e.get("latestVersion", 0)))
            latest_version = int(best.get("latestVersion", 0))
        except (TypeError, ValueError):
            record["detail"] = "latestVersion is not an integer"
            return record, None

    fields = tuple(sorted(best))
    record.update({
        "selected_entity_id": best.get("modelEntityId"),
        "selected_version": latest_version,
        "uniprot_accession": best.get("uniprotAccession", accession),
        "global_metric_value": best.get("globalMetricValue"),
        "fields_present": list(fields),
    })
    for required in ("modelEntityId", "latestVersion", "uniprotAccession"):
        if required not in best:
            record["detail"] = f"API entry lacks required field {required!r}"
            return record, None
    if str(best.get("uniprotAccession", "")).upper() != accession.upper():
        record["detail"] = (
            f"served accession {best.get('uniprotAccession')!r} does not match "
            f"requested accession {accession!r}"
        )
        return record, None
    if record["pinned"]:
        if best.get("pdbUrl") != record["pinned_source_url"]:
            record["detail"] = (
                "pinned source_url is not the served PDB URL for the exact "
                "entity/version"
            )
            return record, None
        source = record["pinned_source_url"]
    else:
        source = best.get("pdbUrl") or best.get("cifUrl")
    if not source:
        record["detail"] = "API entry has neither pdbUrl nor cifUrl"
        return record, None
    try:
        content = _download(source, fetcher=fetcher)
        fmt = detect_structure_format(content)
    except (AcquisitionError, ValueError) as exc:
        record["source"] = source
        record["detail"] = f"structure acquisition failed: {exc}"[:300]
        return record, None
    if fmt != "pdb":
        record["source"] = source
        record["detail"] = (
            "mmCIF structure acquired but this candidate's ingestion path "
            "supports PDB only; no fabricated conversion performed"
        )
        return record, None
    if record["pinned"]:
        if len(content) != expectation["bytes"]:
            record["source"] = source
            record["detail"] = (
                f"pinned byte length {expectation['bytes']} does not match "
                f"served content {len(content)}"
            )
            return record, None
        digest = _sha256(content)
        if digest != expectation["sha256"].lower():
            record["source"] = source
            record["detail"] = (
                f"pinned coordinate digest {expectation['sha256'].lower()} "
                f"does not match served content {digest}"
            )
            return record, None

    record.update({
        "acquisition_status": "acquired",
        "authoritative": bool(record["pinned"]),
        "source": source,
        "sha256": _sha256(content),
        "bytes": len(content),
        "fetched_at_utc": _now(),
        "format": fmt,
    })
    return record, content


def _validate_expectation(record: dict, expectation: dict | None) -> bool:
    if expectation is None:
        record["diagnostic_only"] = True
        return True
    if not isinstance(expectation, dict):
        record["detail"] = "invalid prior digest entry (not an object)"
        return False
    if record.get("requested_entity_id"):
        return _validate_afdb_expectation(record, expectation)
    expected = expectation.get("sha256")
    if not isinstance(expected, str) or not HEX64.fullmatch(expected.lower()):
        record["detail"] = "invalid prior digest (sha256 must be 64 hex chars)"
        return False
    record["pinned"] = True
    return True


def acquire_coordinates(pdb_id: str, source_dir: Path | None,
                        digest_expectation: dict | None = None,
                        fetcher=None) -> tuple[dict, bytes | None]:
    """Acquire one structure and pin its digest. Never fabricates content."""
    record = {
        "pdb_id": pdb_id,
        "acquisition_status": "failed",
        "authoritative": False,
        "pinned": False,
        "source": None,
        "sha256": None,
        "bytes": None,
        "fetched_at_utc": None,
        "detail": None,
        "requested_entity_id": _afdb_entity_id(pdb_id) if pdb_id.startswith("AF-") else None,
    }
    if not _validate_expectation(record, digest_expectation):
        return record, None

    content: bytes | None = None
    if source_dir is not None:
        source_dir = Path(source_dir)
        if not source_dir.is_dir():
            record["detail"] = f"source directory does not exist: {source_dir}"
            return record, None
        if pdb_id.startswith("AF-"):
            af_record, content = _acquire_afdb(
                pdb_id, digest_expectation, fetcher=fetcher)
            if content is None:
                return af_record, None
            if digest_expectation is not None:
                expected = digest_expectation["sha256"].lower()
                actual = hashlib.sha256(content).hexdigest()
                if actual != expected:
                    af_record.update({
                        "acquisition_status": "failed",
                        "detail": (
                            f"coordinate digest mismatch: expected {expected}, "
                            f"got {actual}"
                        ),
                    })
                    return af_record, None
                expected_bytes = digest_expectation.get("bytes")
                if expected_bytes is not None and expected_bytes != len(content):
                    af_record.update({
                        "acquisition_status": "failed",
                        "detail": (
                            f"coordinate byte-length mismatch: expected "
                            f"{expected_bytes}, got {len(content)}"
                        ),
                    })
                    return af_record, None
            return af_record, content
        for suffix in (".pdb", ".cif"):
            candidate = source_dir / f"{pdb_id}{suffix}"
            if candidate.is_file():
                try:
                    content = candidate.read_bytes()
                except OSError as exc:
                    record["detail"] = f"cannot read {candidate}: {exc}"[:300]
                    return record, None
                record["source"] = f"local:{candidate.name}"
                record["fetched_at_utc"] = _now()
                break
        if content is None:
            record["detail"] = f"no local file {pdb_id}.pdb/.cif in {source_dir}"
    elif pdb_id.startswith("AF-"):
        af_record, content = _acquire_afdb(
            pdb_id, digest_expectation, fetcher=fetcher)
        if content is None:
            return af_record, None
        if digest_expectation is not None:
            expected = digest_expectation["sha256"].lower()
            actual = hashlib.sha256(content).hexdigest()
            if actual != expected:
                af_record.update({
                    "acquisition_status": "failed",
                    "detail": f"coordinate digest mismatch: expected {expected}, got {actual}",
                })
                return af_record, None
            expected_bytes = digest_expectation.get("bytes")
            if expected_bytes is not None and expected_bytes != len(content):
                af_record.update({
                    "acquisition_status": "failed",
                    "detail": (
                        f"coordinate byte-length mismatch: expected {expected_bytes}, "
                        f"got {len(content)}"
                    ),
                })
                return af_record, None
        return af_record, content
    else:
        if not PDB_ID.fullmatch(pdb_id):
            record["detail"] = f"invalid PDB identifier: {pdb_id!r}"
            return record, None
        url = RCSB_URL.format(pdb_id=pdb_id)
        record["source"] = url
        try:
            content = _download(url, fetcher=fetcher)
            record["fetched_at_utc"] = _now()
        except AcquisitionError as exc:
            record["detail"] = str(exc)[:300]

    if content is None:
        return record, None
    if len(content) < MIN_STRUCTURE_BYTES:
        record["detail"] = f"structure is implausibly small ({len(content)} bytes)"
        return record, None
    try:
        fmt = detect_structure_format(content)
    except ValueError as exc:
        record["detail"] = str(exc)[:300]
        return record, None
    if fmt != "pdb":
        record["detail"] = "mmCIF is not supported by this ingestion path"
        return record, None

    digest = _sha256(content)
    if digest_expectation is not None:
        expected = digest_expectation["sha256"].lower()
        if digest != expected:
            record["detail"] = f"coordinate digest mismatch: expected {expected}, got {digest}"
            return record, None
        expected_bytes = digest_expectation.get("bytes")
        if expected_bytes is not None and expected_bytes != len(content):
            record["detail"] = (
                f"coordinate byte-length mismatch: expected {expected_bytes}, "
                f"got {len(content)}"
            )
            return record, None

    record.update({
        "acquisition_status": "acquired",
        "authoritative": bool(record["pinned"]),
        "sha256": digest,
        "bytes": len(content),
        "format": fmt,
    })
    return record, content


def measure_target(content: bytes, ingest_mode: str) -> dict:
    fmt = detect_structure_format(content)
    if fmt != "pdb":
        raise ValueError(f"unsupported structure format for measurement: {fmt}")
    try:
        structure = IngestionProcessor.run(content, f"structure.{fmt}", "Audit", ingest_mode)
    except IngestionError as exc:
        raise ValueError(f"no atoms parsed from structure: {exc}") from exc
    if not structure.atoms:
        raise ValueError("no atoms parsed from structure")

    if ingest_mode == "predicted":
        structure.confidence.is_experimental = False
        structure.confidence.method = "predicted"
        structure.confidence.resolution = None

    is_experimental = bool(getattr(structure.confidence, "is_experimental", False))
    method_type = getattr(structure.confidence, "method", "predicted") \
        if is_experimental else "predicted"
    raw, coverage, _ = Tier1Measurements.run_full_audit(structure)
    result = adjudicate_laws(AdjudicationInput(
        raw_measurements=raw,
        coverage=coverage,
        is_experimental=is_experimental,
        method_type=method_type,
        architecture="NONE",
    ))

    law_rows: dict[str, dict] = {}
    for row in result.law_rows:
        law_rows[row["law_id"]] = {
            "observed": row["observed"],
            "threshold": row["threshold"],
            "operator": row["operator"],
            "units": row["units"],
            "status": row["status"],
            "raw_status": row.get("raw_status"),
            "method": row["method"],
            "enforced": row.get("enforced", True),
            "sample_size": row["sample_size"],
            "reason": row.get("reason"),
            "granularity": row.get("granularity"),
            "plddt_advisory": row.get("plddt_advisory"),
        }

    indeterminate_causes = [
        lid for lid, row in law_rows.items()
        if row["method"] == "deterministic"
        and row["status"] in ("NOT_MEASURED", "MEASUREMENT_FAILED")
    ]
    blocking_laws = [
        lid for lid, row in law_rows.items()
        if row["method"] == "deterministic"
        and row["status"] in ("FAIL", "VETO")
    ]

    return {
        "verdict": result.verdict,
        "verdict_reason": result.verdict_reason,
        "coverage": round(coverage, 2),
        "coverage_basis": "residues",
        "is_experimental": is_experimental,
        "method_type": method_type,
        "resolution": getattr(structure.confidence, "resolution", None),
        "atoms": len(structure.atoms),
        "det_passed": result.det_passed,
        "det_total": DETERMINISTIC_COUNT,
        "blocking_laws": blocking_laws,
        "indeterminate_causes": indeterminate_causes,
        "laws": law_rows,
    }


def _dataset_ids(dataset: dict) -> list[str]:
    if not isinstance(dataset, dict) or not dataset:
        raise CalibrationRunError("dataset must be a non-empty object")
    ids: list[str] = []
    for category, entries in dataset.items():
        if not isinstance(category, str) or not isinstance(entries, list):
            raise CalibrationRunError("dataset categories must map to lists")
        for entry in entries:
            if not isinstance(entry, dict) or not isinstance(entry.get("pdb_id"), str):
                raise CalibrationRunError("each target requires a string pdb_id")
            ids.append(entry["pdb_id"])
    if len(ids) != len(set(ids)):
        raise CalibrationRunError("duplicate pdb_id in dataset")
    return ids


def run_calibration(dataset: dict, out_dir: Path,
                    source_dir: Path | None = None,
                    prior_digests: dict | None = None,
                    fetcher=None) -> tuple[dict, int]:
    """Execute a diagnostic run. Returns (run_document, exit_code)."""
    ids = _dataset_ids(dataset)
    prior_digests = prior_digests or {}
    out_dir = Path(out_dir)
    if out_dir.exists():
        raise CalibrationRunError(
            f"refusing to overwrite existing run directory: {out_dir}"
        )

    provenance = _engine_provenance()
    if not provenance["working_tree_clean"]:
        shown = "\n  ".join(provenance["tree_status"][:20])
        raise CalibrationRunError(
            f"refusing to run on a dirty engine tree at {provenance['code_sha'][:8]}:\n"
            f"  {shown}"
        )

    registry_path = Path(__file__).resolve().parents[2] / "governance" / "law_registry.yaml"
    policy_locks = _verify_policy_locks(registry_path)
    try:
        out_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError as exc:
        raise CalibrationRunError(
            f"refusing to overwrite existing run directory: {out_dir}") from exc
    run_manifest = {
        "runner_version": RUNNER_VERSION,
        "run_type": "diagnostic",
        "purpose": (
            "Document corrected-metric distributions. This is not a "
            "threshold validation or held-out evaluation."
        ),
        "run_timestamp_utc": _now(),
        **provenance,
        "backend_version": STATION_METADATA["version"],
        "policy_binding": {
            "canon_hash": LAW_CANON_HASH,
            "law_semantic_hash": policy_locks["law_semantic_sha256"],
            "modality_matrix_sha256": compute_matrix_hash(),
            "law_registry_sha256": policy_locks["law_registry_sha256"],
        },
        "policy_locks": policy_locks,
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "dataset_definition_sha256": _sha256(
            json.dumps(dataset, sort_keys=True).encode("utf-8")),
        "config": {
            "ingest_modes": {"rcsb": "experimental", "alphafold": "predicted"},
            "network_sources": {
                "rcsb": RCSB_URL,
                "alphafold": "metadata-first entity selection; PDB payload required",
            },
            "http_timeout_s": HTTP_TIMEOUT_S,
            "max_structure_bytes": MAX_STRUCTURE_BYTES,
        },
        "expected_target_count": len(ids),
        "digest_verification": {
            "prior_manifest_supplied": bool(prior_digests),
            "pinned_targets": len(prior_digests),
        },
    }
    (out_dir / "run_manifest.json").write_text(
        json.dumps(run_manifest, indent=2, sort_keys=True), encoding="utf-8")

    coord_manifest: dict[str, dict] = {}
    results: dict[str, dict] = {}
    failures: list[str] = []

    for category, entries in dataset.items():
        for entry in entries:
            pdb_id = entry["pdb_id"]
            expectation = prior_digests.get(pdb_id, {}) if prior_digests else None
            acq, content = acquire_coordinates(
                pdb_id, source_dir, expectation, fetcher=fetcher)
            coord_manifest[pdb_id] = acq
            target = {
                "category": category,
                "resolution": entry.get("resolution"),
                "desc": entry.get("desc"),
                "acquisition": acq,
            }
            if acq["acquisition_status"] != "acquired" or content is None:
                target.update({
                    "measurement_status": "NOT_MEASURED",
                    "verdict": None,
                    "laws": None,
                    "detail": acq.get("detail"),
                })
                failures.append(f"{pdb_id}: acquisition failed ({acq.get('detail')})")
            else:
                mode = "predicted" if pdb_id.startswith("AF-") else "experimental"
                try:
                    target.update({
                        "measurement_status": "measured",
                        **measure_target(content, mode),
                    })
                except Exception as exc:
                    target.update({
                        "measurement_status": "MEASUREMENT_FAILED",
                        "verdict": None,
                        "laws": None,
                        "detail": f"{type(exc).__name__}: {exc}"[:300],
                    })
                    failures.append(f"{pdb_id}: measurement failed ({exc})")
            results[pdb_id] = target

    (out_dir / "coordinate_manifest.json").write_text(
        json.dumps(coord_manifest, indent=2, sort_keys=True), encoding="utf-8")

    acquired = sum(
        1 for record in results.values()
        if record["acquisition"]["acquisition_status"] == "acquired"
    )
    measured = sum(
        1 for record in results.values()
        if record["measurement_status"] == "measured"
    )
    not_measured = sum(
        1 for record in results.values()
        if record["measurement_status"] == "NOT_MEASURED"
    )
    measurement_failed = sum(
        1 for record in results.values()
        if record["measurement_status"] == "MEASUREMENT_FAILED"
    )
    unpinned = sum(
        1 for record in results.values()
        if record["acquisition"].get("diagnostic_only") is True
    )
    complete = measured == len(ids) and not failures
    run_document = {
        "meta": {
            **run_manifest,
            "completed_at_utc": _now(),
        },
        "summary": {
            "expected_targets": len(ids),
            "acquired": acquired,
            "measured": measured,
            "not_measured": not_measured,
            "measurement_failed": measurement_failed,
            "unpinned_diagnostic_targets": unpinned,
            "failure_count": len(failures),
            "acquisition_failures": sum(
                1 for record in results.values()
                if record["measurement_status"] == "NOT_MEASURED"
            ),
            "measurement_failures": measurement_failed,
            "run_status": "complete" if complete else "incomplete",
            "verdicts": _verdict_histogram(results),
            "failures": failures,
        },
        "targets": results,
    }
    (out_dir / "calibration_results.json").write_text(
        json.dumps(run_document, indent=2, sort_keys=True), encoding="utf-8")
    # Generate run_integrity_manifest.json for evaluator compatibility
    integrity_manifest = {
        "run_status": "complete" if complete else "incomplete",
        "run_type": "diagnostic",
        "expected_target_count": len(ids),
        "acquired_count": acquired,
        "measured_count": measured,
        "policy_binding": run_manifest["policy_binding"],
        "code_sha": provenance["code_sha"],
        "artifact_digests": {
            "run_manifest.json": _sha256_file(out_dir / "run_manifest.json"),
            "coordinate_manifest.json": _sha256_file(out_dir / "coordinate_manifest.json"),
            "calibration_results.json": _sha256_file(out_dir / "calibration_results.json"),
        },
    }
    (out_dir / "run_integrity_manifest.json").write_text(
        json.dumps(integrity_manifest, indent=2, sort_keys=True), encoding="utf-8")

    write_manifest_sha256(out_dir)
    return run_document, 0 if complete else 1


def _verdict_histogram(results: dict) -> dict:
    histogram: dict[str, int] = {}
    for record in results.values():
        verdict = record.get("verdict") or f"not:{record['measurement_status']}"
        histogram[verdict] = histogram.get(verdict, 0) + 1
    return dict(sorted(histogram.items()))


def write_manifest_sha256(run_dir: Path) -> None:
    run_dir = Path(run_dir)
    lines = []
    for path in sorted(run_dir.rglob("*")):
        if not path.is_file() or path.name == "MANIFEST.sha256":
            continue
        lines.append(
            f"{_sha256(path.read_bytes())}  {path.relative_to(run_dir).as_posix()}"
        )
    (run_dir / "MANIFEST.sha256").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dataset", default="v23", choices=["v23"])
    parser.add_argument("--out", default=None,
                        help="new output run directory (never overwritten)")
    parser.add_argument("--source-dir", default=None,
                        help="local directory of pre-fetched <pdb_id>.pdb files")
    parser.add_argument("--verify-digests-from", default=None,
                        help="prior coordinate_manifest.json to verify digests against")
    args = parser.parse_args(argv)

    if args.dataset == "v23":
        from calibration.dataset_v23 import CALIBRATION_DATASET
        dataset = CALIBRATION_DATASET
    else:
        raise SystemExit(2)

    prior = None
    if args.verify_digests_from:
        prior_path = Path(args.verify_digests_from)
        try:
            prior = json.loads(prior_path.read_text(encoding="utf-8"))
        except Exception as exc:
            print(f"cannot read prior digest manifest: {exc}", file=sys.stderr)
            return 2
        if not isinstance(prior, dict):
            print("prior digest manifest must be an object", file=sys.stderr)
            return 2

    out = Path(args.out) if args.out else (
        Path(__file__).resolve().parent / "runs" / "v25")
    try:
        _, code = run_calibration(
            dataset, out,
            source_dir=Path(args.source_dir) if args.source_dir else None,
            prior_digests=prior,
        )
    except CalibrationRunError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    return code


if __name__ == "__main__":
    raise SystemExit(main())
