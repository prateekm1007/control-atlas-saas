"""
WP-SCI-1 v24 Diagnostic Calibration Runner (directive §4)
==========================================================

Executes a versioned, evidence-bound calibration run of the corrected
law engine and writes a canonical evidence layout under
`backend/calibration/runs/<version>/`.

Design constraints (from the 2026-10-19 scientific-validity directive):

- Secret handling: the runner measures IN-PROCESS through the same
  functions the production API path calls
  (IngestionProcessor.run -> Tier1Measurements.run_full_audit ->
  adjudicate_laws). It performs NO authenticated network API calls and
  therefore needs no API key; it never reads, prints, or logs secret
  material and never enables open mode.
- Coordinate pinning: every acquired coordinate file is bound by its
  SHA-256 digest in `coordinate_manifest.json`. Re-runs verify the
  digest of locally supplied files; remote sources are recorded with the
  exact URL and fetch timestamp.
- Historical isolation: the v23 artifacts (backend/calibration/
  calibration_results.json, calibration_report_v23.md) and the legacy
  root-level calibration_results.json are NEVER overwritten.
- Fail-closed completeness: the run exits non-zero unless every expected
  target was acquired AND measured. Failed acquisition is recorded
  explicitly and never counted as a successful calibration case.
- No threshold decisions are made here. This run is DIAGNOSTIC: its
  distributions document the corrected metrics; thresholds may only be
  re-validated through the independent evaluation protocol under
  backend/calibration/protocol/.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import tos  # noqa: E402  (provenance: the ENGINE actually imported)
from tos.engine.adjudicator import adjudicate_laws, AdjudicationInput
from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import IngestionProcessor
from tos.governance.station_sop import (
    LAW_CANON,
    LAW_METHOD_CLASSIFICATIONS,
    STATION_METADATA,
    LAW_CANON_HASH,
    compute_law_semantic_hash,
)
from tos.governance.modality_matrix import compute_matrix_hash
from tos.generation.afdb import (  # noqa: E402
    AFDBUnavailable,
    detect_structure_format,
    fetch_prediction,
    uniprot_entry_type,
)

RUNNER_VERSION = "v25-runner/1.0.0"
RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"
HTTP_TIMEOUT_S = 60
USER_AGENT = ("toscanini-calibration-runner/1.0 (scientific evidence run; "
              "contact: repository owner)")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _engine_repo_root() -> Path | None:
    """Repo root of the git checkout that OWNS the imported tos package.

    Provenance v2 (Round-12 direction 4, package-install audit): a
    repository that merely CONTAINS the installed package (for example a
    virtualenv or site-packages directory nested inside an unrelated
    checkout) is NOT the source of the imported code and must never be
    reported as such. Ownership requires BOTH:
      1. the package file lives inside the candidate repo's work tree, and
      2. git tracks the package file in that repo (ls-files --error-unmatch).
    Returns None for installed packages (wheel/site-packages), packages
    with no Git metadata, and untracked copies — the fail-closed answer.
    """
    try:
        pkg_file = Path(tos.__file__).resolve()
    except Exception:
        return None
    for parent in pkg_file.parents:
        if not (parent / ".git").exists():
            continue
        try:
            root = subprocess.run(
                ["git", "-C", str(parent), "rev-parse", "--show-toplevel"],
                capture_output=True, text=True, check=True,
            ).stdout.strip()
        except Exception:
            return None
        root_path = Path(root).resolve()
        try:
            rel = pkg_file.relative_to(root_path)
        except ValueError:
            return None  # package file is outside this repo's work tree
        tracked = subprocess.run(
            ["git", "-C", str(root_path), "ls-files", "--error-unmatch",
             str(rel)],
            capture_output=True, text=True,
        )
        if tracked.returncode != 0:
            return None  # repo contains the dir but does not own the package
        return root_path
    return None


def _installed_distribution() -> dict | None:
    """Best-effort installed-distribution record for the tos package
    (Provenance v2: known installed-artifact route). Returns None when
    the imported package is not pip-installed (normal source checkout)."""
    try:
        from importlib import metadata
    except ImportError:  # pragma: no cover
        return None
    try:
        dist_names = metadata.packages_distributions().get("tos", [])
    except Exception:
        return None
    for name in dist_names:
        try:
            dist = metadata.distribution(name)
        except Exception:
            continue
        record = {"name": dist.metadata["Name"], "version": dist.version}
        try:
            direct = dist.read_text("direct_url.json")
            if direct:
                parsed = json.loads(direct)
                # never carry raw URLs with credentials into artifacts
                if "url" in parsed and isinstance(parsed["url"], str):
                    parsed["url"] = _sanitize_git_url(parsed["url"])
                record["direct_url"] = parsed
        except Exception:
            pass
        return record
    return None


def _git_output(*args: str) -> str | None:
    """Best-effort git metadata of the git checkout that OWNS the
    imported tos package; returns None outside such a checkout
    (installed packages, no Git metadata, or untracked copies)."""
    root = _engine_repo_root()
    if root is None:
        return None
    try:
        return subprocess.run(
            ["git", "-C", str(root), *args],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except Exception:
        return None


def _engine_provenance(strict: bool = True) -> dict:
    """P0 item 15: record the commit of the code that will ACTUALLY run.

    The v24 run recorded the runner checkout's HEAD while the working
    tree contained uncommitted engine changes, so its report numbers
    were untethered from any reproducible commit. Provenance is now
    taken from the git checkout that OWNS the imported tos package
    (verified by work-tree containment AND git tracking — a checkout
    that merely contains an installed copy is rejected), and the tree
    state is RECORDED (working_tree_clean + tree_status).

    Provenance v2 routes:
    - git-checkout: code_sha + tree state (above).
    - installed-package: the distribution metadata is recorded and
      code_sha is None — an installed artifact is never attributed a
      commit inferred from a containing checkout.
    - unavailable: neither route resolved.

    strict=True (authoritative measurement): anything but the
    git-checkout route FAILS CLOSED (SystemExit) — no evidence run on
    unverifiable source provenance.
    strict=False (records): the unavailable/installed state is recorded
    explicitly and never guessed.
    """
    code_sha = _git_output("rev-parse", "HEAD")
    if code_sha is not None:
        status = _git_output("status", "--porcelain") or ""
        return {
            "code_sha": code_sha,
            "code_sha_source": (
                "git rev-parse HEAD of the repository that OWNS the imported "
                "tos engine package (work-tree containment + git tracking "
                "verified)"
            ),
            "source_provenance": "git-checkout",
            "working_tree_clean": not status,
            "tree_status": status.splitlines()[:50],
        }
    # No git checkout owns the imported package: try the installed-
    # artifact route, then record unavailability explicitly. Never
    # infer a commit from a checkout that does not own the package.
    dist = _installed_distribution()
    if dist is not None:
        record = {
            "code_sha": None,
            "code_sha_source": (
                "unavailable — the imported tos package is pip-installed; "
                "no git checkout owns it. Distribution metadata recorded; "
                "a commit is never inferred from a containing checkout."
            ),
            "source_provenance": "installed-package",
            "installed_distribution": dist,
            "working_tree_clean": None,
            "tree_status": None,
        }
        if strict:
            raise SystemExit(
                "authoritative calibration runs require a git checkout "
                "that OWNS the imported tos engine package; an installed "
                "package has no verifiable source commit")
        return record
    record = {
        "code_sha": None,
        "code_sha_source": (
            "unavailable — the imported tos package is neither owned by a "
            "git checkout (work-tree containment + git tracking) nor "
            "installed with distribution metadata"
        ),
        "source_provenance": "unavailable",
        "working_tree_clean": None,
        "tree_status": None,
    }
    if strict:
        raise SystemExit(
            "calibration runs require a git checkout of the engine "
            "(could not resolve a repository that OWNS the imported tos "
            "package)")
    return record


def _sanitize_git_url(url: str) -> str:
    """Strip credentials from a git remote URL (EC §3: no secrets in
    reports, artifacts, or logs). Handles scheme URLs with userinfo
    (https://user:pass@host/..., https://token@host/...) and scp-style
    URLs (git@host:path)."""
    import re
    from urllib.parse import urlsplit, urlunsplit
    # scp-style first (no scheme): user@host:path
    url = re.sub(r"^[^/@\s]+@(?=[^/]+:)", "", url)
    try:
        parts = urlsplit(url)
    except ValueError:
        return "<unparseable remote url>"
    if parts.username is None and parts.password is None:
        return url
    host = parts.hostname or ""
    if parts.port:
        host = f"{host}:{parts.port}"
    return urlunsplit((parts.scheme, host, parts.path,
                       parts.query, parts.fragment))


def _git_context() -> dict:
    """Sanitized git context of the engine repository. Remote URLs are
    credential-stripped before they enter any artifact."""
    branch = _git_output("rev-parse", "--abbrev-ref", "HEAD")
    describe = _git_output("describe", "--tags", "--always", "--dirty")
    remotes_raw = _git_output("remote", "-v") or ""
    remotes: dict[str, str] = {}
    for line in remotes_raw.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] not in remotes:
            remotes[parts[0]] = _sanitize_git_url(parts[1])
    return {
        "available": _engine_repo_root() is not None,
        "branch": branch,
        "describe": describe,
        "remotes": remotes,
        "note": (
            "git context of the repository that OWNS the imported tos "
            "package; null when the package is installed or has no Git "
            "metadata — never inferred from a containing checkout"
        ),
    }


def _dependency_versions() -> dict:
    """Best-effort versions of dependencies material to measurement
    (EC §6: dependency/runtime versions where material). A dependency
    that cannot be imported is recorded explicitly as unavailable —
    never silently omitted."""
    out: dict[str, str] = {}
    for name in ("numpy", "gemmi"):
        try:
            module = __import__(name)
            out[name] = str(getattr(module, "__version__", "unknown"))
        except Exception as exc:  # pragma: no cover - environment-specific
            out[name] = f"unavailable ({type(exc).__name__})"
    return out


def build_provenance_manifest(calibration_executed: bool = False) -> dict:
    """Provenance v2: complete standalone provenance RECORD.

    Composes the engine provenance (git checkout / installed package /
    explicitly-unavailable) with the sanitized git context, the runtime,
    the material dependency versions, and the policy-binding hashes.
    Safe to persist anywhere: it contains no credentials (remotes are
    sanitized; no environment variables or tokens are recorded).

    calibration_executed is recorded explicitly (default False): a
    provenance record must never be mistakable for a calibration run —
    the demonstration manifest records False even though it identifies
    the runner version.
    """
    return {
        "schema": "toscanini-runner-provenance/2.0.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "runner_version": RUNNER_VERSION,
        "calibration_executed": bool(calibration_executed),
        "record_type": (
            "calibration-run" if calibration_executed
            else "provenance-record-only"),
        "engine": _engine_provenance(strict=False),
        "git_context": _git_context(),
        "runtime": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "platform": platform.platform(),
        },
        "dependencies": _dependency_versions(),
        "policy_binding": {
            "canon_hash": LAW_CANON_HASH,
            "law_semantic_hash": compute_law_semantic_hash(),
            "modality_matrix_sha256": compute_matrix_hash(),
        },
    }


def write_provenance_manifest(out_path: Path) -> dict:
    """Write the standalone provenance manifest as JSON; returns it."""
    manifest = build_provenance_manifest()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def target_source_url(pdb_id: str) -> str:
    if pdb_id.startswith("AF-"):
        # P0 item 11: metadata-first acquisition; the file URL is only
        # known after selecting the entry by modelEntityId/latestVersion.
        accession = pdb_id.replace("AF-", "").split("-F")[0]
        return AFDB_API.format(uniprot=accession)
    return RCSB_URL.format(pdb_id=pdb_id)


def _acquire_afdb_remote(pdb_id: str, record: dict) -> bytes | None:
    """Metadata-first AF acquisition; updates the record in place."""
    try:
        pred = fetch_prediction(pdb_id)
    except AFDBUnavailable as exc:
        record["detail"] = f"{type(exc).__name__}: {exc}"[:300]
        record["afdb_record"] = exc.record
        if exc.record.get("migrated_identifier"):
            record["migrated_identifier"] = exc.record["migrated_identifier"]
        # Item 14: record WHY (informational UniProt status, e.g. the
        # P62988 merged-entry case) so cohort decisions are auditable.
        accession = pdb_id.replace("AF-", "").split("-F")[0]
        entry_type = uniprot_entry_type(accession)
        if entry_type is not None:
            record["uniprot_entry_type"] = entry_type
        return None
    record.update({
        "acquisition_status": "acquired",
        "source": pred.source_url,
        "sha256": pred.sha256,
        "bytes": len(pred.content),
        "afdb": {
            "model_entity_id": pred.model_entity_id,
            "latest_version": pred.latest_version,
            "requested_entity_id": pred.requested_entity_id,
            "uniprot_accession": pred.uniprot_accession,
            "global_metric_value": pred.global_metric_value,
            "format": pred.fmt,
            "fields_present": list(pred.fields_present),
        },
    })
    return pred.content


def acquire_coordinates(pdb_id: str, source_dir: Path | None,
                        digest_expectation: dict | None) -> tuple[dict, bytes | None]:
    """Fetch (or read locally) one coordinate file and pin its digest.

    Returns (acquisition_record, content_or_None); never raises for
    missing/failed acquisition — the caller records the failure
    explicitly.
    """
    record = {"pdb_id": pdb_id, "acquisition_status": "failed",
              "source": None, "sha256": None, "bytes": None,
              "fetched_at_utc": None, "detail": None}
    data: bytes | None = None

    if source_dir is not None:
        for candidate in (f"{pdb_id}.pdb", f"{pdb_id}.cif"):
            p = source_dir / candidate
            if p.exists():
                data = p.read_bytes()
                record["source"] = f"local:{p.name}"
                record["fetched_at_utc"] = datetime.now(
                    timezone.utc).isoformat()
                break
        if data is None:
            record["detail"] = f"no local file {pdb_id}.pdb/.cif in {source_dir}"
            return record, None
    elif pdb_id.startswith("AF-"):
        # P0 item 11: metadata-first AFDB acquisition (no hardcoded
        # model_v4 URLs). The record carries URL, entity, version,
        # SHA-256, and — on failure — migrated identifiers and the
        # UniProt status (item 14).
        data = _acquire_afdb_remote(pdb_id, record)
        if data is None:
            return record, None
        record["fetched_at_utc"] = datetime.now(timezone.utc).isoformat()
    else:
        url = target_source_url(pdb_id)
        record["source"] = url
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=HTTP_TIMEOUT_S) as resp:
                data = resp.read()
            record["fetched_at_utc"] = datetime.now(timezone.utc).isoformat()
        except Exception as exc:  # network failure, HTTP error, timeout
            record["detail"] = f"{type(exc).__name__}: {exc}"[:300]
            return record, None

    digest = _sha256(data)
    if digest_expectation and digest_expectation.get("sha256") not in (None, digest):
        record["detail"] = (
            f"coordinate digest mismatch: expected "
            f"{digest_expectation.get('sha256')}, got {digest}"
        )
        return record, None

    record.update({"acquisition_status": "acquired",
                   "sha256": digest, "bytes": len(data)})
    return record, data


def measure_target(content: bytes, ingest_mode: str) -> dict:
    """Run the production measurement path for one structure.

    Mirrors backend/main.py::_run_physics_sync (stages 1-4) without any
    HTTP, persistence, PDF, or narrative side effects.

    P0 item 15: the input filename is derived from the CONTENT-detected
    format (the hardcoded "origin.pdb" mislabeled mmCIF input).
    """
    fmt = detect_structure_format(content)
    structure = IngestionProcessor.run(
        content, f"structure.{fmt}", "Audit", ingest_mode)
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

    law_rows = {}
    for row in result.law_rows:
        law_rows[row["law_id"]] = {
            "observed": row["observed"],
            "threshold": row["threshold"],
            "operator": row["operator"],
            "units": row["units"],
            "status": row["status"],
            "method": row["method"],
            "sample_size": row["sample_size"],
            **({"contained": True, "containment_reason": row["containment_reason"]}
               if row.get("contained") else {}),
            **({"reason": row["reason"]} if row.get("reason") else {}),
        }

    return {
        "verdict": result.verdict,
        "coverage": round(coverage, 2),
        "is_experimental": is_experimental,
        "method_type": method_type,
        "resolution": getattr(structure.confidence, "resolution", None),
        "atoms": len(structure.atoms),
        "det_passed": result.det_passed,
        "det_total": len(LAW_CANON),
        "contained_law_failures": result.contained_law_failures,
        "indeterminate_causes": result.indeterminate_causes,
        "laws": law_rows,
    }


def run_calibration(dataset: dict, out_dir: Path,
                    source_dir: Path | None = None,
                    prior_digests: dict | None = None) -> tuple[dict, int]:
    """Execute the full calibration. Returns (run_document, exit_code)."""
    # P0 item 15: never overwrite an existing run — evidence runs are
    # append-only. A re-run must choose a new directory.
    if out_dir.exists():
        raise SystemExit(
            f"refusing to overwrite existing run directory: {out_dir}. "
            "Calibration runs are evidence artifacts; use a new directory.")
    out_dir.mkdir(parents=True)
    coord_manifest_path = out_dir / "coordinate_manifest.json"
    prior_digests = prior_digests or {}

    provenance = _engine_provenance()
    run_manifest = {
        "runner_version": RUNNER_VERSION,
        "run_type": "diagnostic",
        "calibration_executed": True,
        "purpose": (
            "Document corrected-metric distributions. NOT a "
            "threshold validation: thresholds may only be re-validated via "
            "backend/calibration/protocol/evaluation_protocol_v1.md."
        ),
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        **provenance,
        "backend_version": STATION_METADATA["version"],
        "policy_binding": {
            "canon_hash": LAW_CANON_HASH,
            "law_semantic_hash": compute_law_semantic_hash(),
            "modality_matrix_sha256": compute_matrix_hash(),
            "law_registry_sha256": (_sha256(
                (Path(__file__).resolve().parents[2] / "governance" /
                 "law_registry.yaml").read_bytes()
            ) if (_git_output("rev-parse", "HEAD") is not None) else None),
        },
        "environment": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "platform": platform.platform(),
        },
        "git": _git_context(),
        "dependencies": _dependency_versions(),
        "dataset_definition_sha256": _sha256(
            json.dumps(dataset, sort_keys=True).encode("utf-8")),
        "config": {
            "ingest_modes": {"rcsb": "experimental", "alphafold": "predicted"},
            "network_sources": {
                "rcsb": RCSB_URL,
                "alphafold": (
                    "metadata-first: " + AFDB_API +
                    " -> select by modelEntityId + latestVersion -> cifUrl"
                ),
            },
            "http_timeout_s": HTTP_TIMEOUT_S,
        },
        "expected_target_count": sum(len(v) for v in dataset.values()),
    }
    (out_dir / "run_manifest.json").write_text(
        json.dumps(run_manifest, indent=2), encoding="utf-8")

    coord_manifest: dict[str, dict] = {}
    results: dict[str, dict] = {}
    failures: list[str] = []

    for category, entries in dataset.items():
        for entry in entries:
            pdb_id = entry["pdb_id"]
            acq, content = acquire_coordinates(
                pdb_id, source_dir, prior_digests.get(pdb_id))
            coord_manifest[pdb_id] = acq

            target_record = {
                "category": category,
                "resolution": entry.get("resolution"),
                "desc": entry.get("desc"),
                "acquisition": acq,
            }
            if acq["acquisition_status"] != "acquired":
                target_record.update({
                    "measurement_status": "NOT_MEASURED",
                    "verdict": None,
                    "laws": None,
                    "detail": acq.get("detail"),
                })
                failures.append(f"{pdb_id}: acquisition failed ({acq.get('detail')})")
            else:
                ingest_mode = "predicted" if pdb_id.startswith("AF-") else "experimental"
                try:
                    measurement = measure_target(content, ingest_mode)
                    target_record.update({
                        "measurement_status": "measured",
                        **measurement,
                    })
                except Exception as exc:
                    target_record.update({
                        "measurement_status": "MEASUREMENT_FAILED",
                        "verdict": None,
                        "laws": None,
                        "detail": f"{type(exc).__name__}: {exc}"[:300],
                    })
                    failures.append(f"{pdb_id}: measurement failed ({exc})")
            results[pdb_id] = target_record

    coord_manifest_path.write_text(
        json.dumps(coord_manifest, indent=2, sort_keys=True), encoding="utf-8")

    acquired = sum(1 for r in results.values()
                   if r["acquisition"]["acquisition_status"] == "acquired")
    measured = sum(1 for r in results.values() if r["measurement_status"] == "measured")
    # P0 item 15: complete ONLY at the full expected target count with
    # zero failures — e.g. 50 of 50 for the v25 cohort.
    all_measured = (measured == run_manifest["expected_target_count"])
    run_document = {
        "meta": {
            **run_manifest,
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        },
        "summary": {
            "expected_targets": run_manifest["expected_target_count"],
            "acquired": acquired,
            "measured": measured,
            "acquisition_failures": len(failures),
            "run_status": "complete" if (all_measured and not failures) else "incomplete",
            "verdicts": _verdict_histogram(results),
            "failures": failures,
        },
        "targets": results,
    }
    (out_dir / "calibration_results.json").write_text(
        json.dumps(run_document, indent=2), encoding="utf-8")

    write_manifest_sha256(out_dir)

    exit_code = 0 if not failures else 1
    return run_document, exit_code


def _verdict_histogram(results: dict) -> dict:
    hist: dict[str, int] = {}
    for r in results.values():
        v = r.get("verdict") or f"not:{r['measurement_status']}"
        hist[v] = hist.get(v, 0) + 1
    return dict(sorted(hist.items()))


def write_manifest_sha256(run_dir: Path) -> None:
    """SHA-256 manifest over every file in the run dir except itself,
    sorted by relative path (deterministic)."""
    lines = []
    for p in sorted(run_dir.rglob("*")):
        if not p.is_file():
            continue
        if p.name == "MANIFEST.sha256":
            continue
        digest = _sha256(p.read_bytes())
        lines.append(f"{digest}  {p.relative_to(run_dir).as_posix()}")
    (run_dir / "MANIFEST.sha256").write_text("\n".join(lines) + "\n",
                                             encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dataset", default="v23",
                        choices=["v23"], help="dataset definition to run")
    parser.add_argument("--out", default=str(here / "runs" / "v25"),
                        help="output run directory (must not already exist)")
    parser.add_argument("--source-dir", default=None,
                        help="local directory of pre-fetched coordinate files "
                             "(<pdb_id>.pdb); skips network acquisition")
    parser.add_argument("--verify-digests-from", default=None,
                        help="prior coordinate_manifest.json to verify digests against")
    parser.add_argument("--emit-provenance", default=None, metavar="PATH",
                        help="write a standalone provenance manifest to PATH "
                             "and exit 0 without running a calibration "
                             "(records the tree state faithfully instead "
                             "of refusing; remotes are credential-stripped)")
    args = parser.parse_args(argv)

    # U9 (Round-11): provenance emission is a RECORD, not an evidence run —
    # it never refuses a dirty tree, it records it. Useful for CI stamps
    # and for demonstrating provenance capture without a full run.
    if args.emit_provenance:
        manifest = write_provenance_manifest(Path(args.emit_provenance))
        print(json.dumps({
            "written": args.emit_provenance,
            "schema": manifest["schema"],
            "code_sha": manifest["engine"]["code_sha"],
            "working_tree_clean": manifest["engine"]["working_tree_clean"],
        }))
        return 0

    # P0 item 15: evidence runs refuse a dirty engine tree — the recorded
    # code_sha must describe the code that actually measured. (This is
    # exactly the v24 provenance hole: post-fix numbers under a pre-fix SHA.)
    provenance = _engine_provenance()
    if not provenance["working_tree_clean"]:
        shown = "\n  ".join(provenance["tree_status"][:20])
        print(
            "refusing to run a calibration on a dirty working tree:\n"
            f"  engine repo @ {provenance['code_sha'][:8]}:\n  {shown}\n"
            "  commit or stash, then re-run.", file=sys.stderr)
        return 2

    if args.dataset == "v23":
        from calibration.dataset_v23 import CALIBRATION_DATASET
        dataset = CALIBRATION_DATASET

    prior = None
    if args.verify_digests_from:
        prior = json.loads(Path(args.verify_digests_from).read_text())
        prior = {k: v for k, v in prior.items() if isinstance(v, dict)}

    _, exit_code = run_calibration(
        dataset, Path(args.out),
        source_dir=Path(args.source_dir) if args.source_dir else None,
        prior_digests=prior,
    )
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
