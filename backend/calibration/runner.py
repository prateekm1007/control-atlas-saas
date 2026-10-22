"""
Diagnostic calibration runner with evidence-bound acquisition and output.

The runner measures IN-PROCESS through the same engine functions used by the
API path. It performs no authenticated API calls, reads no secrets, and never
enables open mode.

Boundaries enforced here:

* every coordinate file is hashed and pinned before measurement;
* local files must match a supplied prior digest when one is provided;
* AlphaFold acquisition is metadata-first and entity-exact; a moved or
  ambiguous entity fails closed instead of being silently substituted;
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
from urllib.parse import urlparse

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import tos  # noqa: E402
from tos.engine.adjudicator import AdjudicationInput, adjudicate_laws  # noqa: E402
from tos.engine.tier1_measurements import Tier1Measurements  # noqa: E402
from tos.ingestion.processor import IngestionProcessor  # noqa: E402
from tos.governance.station_sop import (  # noqa: E402
    DETERMINISTIC_COUNT,
    LAW_CANON,
    LAW_CANON_HASH,
    LAW_METHOD_CLASSIFICATIONS,
    STATION_METADATA,
)
from tos.governance.modality_matrix import compute_matrix_hash  # noqa: E402

RUNNER_VERSION = "v25-runner/1.1.0"
RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
HTTP_TIMEOUT_S = 60
MAX_STRUCTURE_BYTES = 50 * 1024 * 1024
MIN_STRUCTURE_BYTES = 200
USER_AGENT = "toscanini-calibration-runner/1.1 (diagnostic evidence run)"
ALLOWED_HOSTS = {"files.rcsb.org", "alphafold.ebi.ac.uk"}
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


def _assert_https(url: str) -> None:
    parsed = urlparse(url)
    if parsed.scheme != "https":
        raise AcquisitionError(f"non-HTTPS URL refused: {url}")
    if parsed.hostname not in ALLOWED_HOSTS:
        raise AcquisitionError(f"host not allowed: {parsed.hostname}")


def _http_get(url: str, timeout: int = HTTP_TIMEOUT_S):
    return requests.get(
        url, timeout=timeout, headers={"User-Agent": USER_AGENT},
        allow_redirects=True,
    )


def _download(url: str, fetcher=None, min_bytes: int = MIN_STRUCTURE_BYTES) -> bytes:
    _assert_https(url)
    get = fetcher or _http_get
    try:
        response = get(url, timeout=HTTP_TIMEOUT_S)
    except Exception as exc:
        raise AcquisitionError(f"request failed: {type(exc).__name__}: {exc}") from exc
    status = getattr(response, "status_code", None)
    if status != 200:
        raise AcquisitionError(f"HTTP {status} from {url}")
    final_url = getattr(response, "url", url) or url
    _assert_https(final_url)
    content = getattr(response, "content", b"")
    if not isinstance(content, (bytes, bytearray)):
        raise AcquisitionError("response body is not bytes")
    content = bytes(content)
    if len(content) > MAX_STRUCTURE_BYTES:
        raise AcquisitionError(f"structure exceeds {MAX_STRUCTURE_BYTES} bytes")
    if len(content) < min_bytes:
        raise AcquisitionError(f"payload is implausibly small ({len(content)} bytes)")
    return content


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


def _afdb_record(pdb_id: str) -> dict:
    accession = normalize_accession(pdb_id)
    entity = f"AF-{accession}-F1"
    return {
        "pdb_id": pdb_id,
        "acquisition_status": "failed",
        "source": None,
        "sha256": None,
        "bytes": None,
        "fetched_at_utc": None,
        "detail": None,
        "api_url": AFDB_API.format(accession=accession),
        "requested_entity_id": entity,
    }


def _acquire_afdb(pdb_id: str, fetcher=None) -> tuple[dict, bytes | None]:
    record = _afdb_record(pdb_id)
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
        best = max(matches, key=lambda e: int(e.get("latestVersion", 0)))
        latest_version = int(best.get("latestVersion", 0))
    except (TypeError, ValueError):
        record["detail"] = "latestVersion is not an integer"
        return record, None

    fields = tuple(sorted(best))
    record.update({
        "selected_entity_id": best.get("modelEntityId"),
        "latest_version": latest_version,
        "uniprot_accession": best.get("uniprotAccession", accession),
        "global_metric_value": best.get("globalMetricValue"),
        "fields_present": list(fields),
    })
    for required in ("modelEntityId", "latestVersion", "uniprotAccession"):
        if required not in best:
            record["detail"] = f"API entry lacks required field {required!r}"
            return record, None

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

    record.update({
        "acquisition_status": "acquired",
        "source": source,
        "sha256": _sha256(content),
        "bytes": len(content),
        "fetched_at_utc": _now(),
        "format": fmt,
    })
    return record, content


def _validate_expectation(record: dict, expectation: dict | None) -> bool:
    if expectation is None:
        return True
    if not isinstance(expectation, dict):
        record["detail"] = "invalid prior digest entry (not an object)"
        return False
    expected = expectation.get("sha256")
    if not isinstance(expected, str) or not HEX64.fullmatch(expected.lower()):
        record["detail"] = "invalid prior digest (sha256 must be 64 hex chars)"
        return False
    return True


def acquire_coordinates(pdb_id: str, source_dir: Path | None,
                        digest_expectation: dict | None = None,
                        fetcher=None) -> tuple[dict, bytes | None]:
    """Acquire one structure and pin its digest. Never fabricates content."""
    record = {
        "pdb_id": pdb_id,
        "acquisition_status": "failed",
        "source": None,
        "sha256": None,
        "bytes": None,
        "fetched_at_utc": None,
        "detail": None,
    }
    if not _validate_expectation(record, digest_expectation):
        return record, None

    content: bytes | None = None
    if source_dir is not None:
        source_dir = Path(source_dir)
        if not source_dir.is_dir():
            record["detail"] = f"source directory does not exist: {source_dir}"
            return record, None
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
        af_record, content = _acquire_afdb(pdb_id, fetcher=fetcher)
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
        "sha256": digest,
        "bytes": len(content),
        "format": fmt,
    })
    return record, content


def measure_target(content: bytes, ingest_mode: str) -> dict:
    fmt = detect_structure_format(content)
    if fmt != "pdb":
        raise ValueError(f"unsupported structure format for measurement: {fmt}")
    structure = IngestionProcessor.run(content, f"structure.{fmt}", "Audit", ingest_mode)
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

    out_dir.mkdir(parents=True)
    registry_path = Path(__file__).resolve().parents[2] / "governance" / "law_registry.yaml"
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
            "law_semantic_hash": _law_semantic_hash(),
            "modality_matrix_sha256": compute_matrix_hash(),
            "law_registry_sha256": _sha256_file(registry_path) if registry_path.is_file() else None,
        },
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
            "acquisition_failures": len(failures),
            "run_status": "complete" if complete else "incomplete",
            "verdicts": _verdict_histogram(results),
            "failures": failures,
        },
        "targets": results,
    }
    (out_dir / "calibration_results.json").write_text(
        json.dumps(run_document, indent=2, sort_keys=True), encoding="utf-8")
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
