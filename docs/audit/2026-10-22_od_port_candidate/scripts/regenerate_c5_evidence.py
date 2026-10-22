#!/usr/bin/env python3
"""
C5 evidence regeneration (round-26 direction 5, candidate review branch).

Deterministically re-derives the v24 verdict-change evidence with PORTABLE
provenance. No machine-specific paths are written into any output artifact;
only revision hashes, pinned URL sources, byte counts and SHA-256 digests.

What it does
------------

1. Materializes two engine checkouts byte-exactly from git objects
   (``git archive``): the OLD base revision (``c78ef15``) and the NEW
   candidate revision (default: current HEAD).
2. Reads the v24 manifests from the immutable delivery ref (``97d62a8``)
   and hard-verifies their SHA-256 pins before any use:
   coordinate_manifest  9722c3e45d316ff58ed326df074eec138b81ceb520191a492088cc0ddef63762
   run_manifest         f602351f4844531d409cfc256939508103e3c2d4f57fe21d7f77c34494aa4036
3. Re-fetches the 40 v24-acquired coordinates from their pinned source
   URLs and verifies each against the recorded byte count and SHA-256.
   The 10 AlphaFold ``model_v4`` targets are probed once and recorded as
   unavailable (expected HTTP 404).
4. Fetches the supplementary predicted cohort via the AlphaFold DB
   metadata-first API (current model version; expected v6). mmCIF bytes
   are converted to PDB records by a minimal deterministic converter
   inside this script: BOTH engines take byte-identical converted input,
   so verdict deltas can only come from engine differences.
5. Evaluates every structure through an identical ingestion/adjudication
   path mirroring ``backend/main.py``'s canonical call (_run_physics_sync):
   ``IngestionProcessor.run(content, "origin.pdb", "Audit", "Audit")`` ->
   ``Tier1Measurements.run_full_audit(structure, user_intent="NONE")`` ->
   ``adjudicate_laws(AdjudicationInput(raw_measurements, coverage,
   is_experimental, method_type, architecture="NONE"))``.
   Each evaluation runs in an isolated subprocess with exactly one engine
   checkout on the import path.
6. Writes the regenerated evidence set (sorted-key JSON, no machine paths):
   ``c5/v24-pinned-old-c78.json``, ``c5/v24-pinned-new-<rev12>.json``,
   ``c5/af-v6-supp-old-c78.json``, ``c5/af-v6-supp-new-<rev12>.json``,
   ``c5/v24-coordinate-verification.json``,
   ``c5/AF-<target>.json`` (acquisition records), and
   ``C5-verdict-change-<rev12>.md`` (regenerated narrative).

Status label (unchanged from the frozen delivery record): diagnostic
evidence for owner review; nothing here validates thresholds, cohorts or
the historical verdict-change table (OD-14 regime; no data/threshold edits).

Usage
-----
Run from inside the candidate work tree (network required for the fetch
phase; re-runs are idempotent thanks to the sha-pinned download cache)::

    python docs/audit/2026-10-22_od_port_candidate/scripts/regenerate_c5_evidence.py \
        --work-root <abs scratch dir>

Add ``--probe`` to materialize the two engines and evaluate one small
experimental target plus one supplementary target for pipeline validation
before the full run.
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import urllib.error
import urllib.request

OLD_REV = "c78ef15f696a2b1e2658e7ac7b54461f70c6e040"
# The candidate branch's own base. The round-25 delta is exactly
# base..new, so base-vs-new is the attributed comparison. OLD_REV stays
# as a cross-lineage documentary check (the od-scope lineage already
# carries the P0/WP-SCI-1 items, which are NOT part of the ported scope).
BASE_REV = "63da0e347002ed9b77071ff1ed0cefda7c0ce88a"
DEFAULT_V24_REF = "97d62a89f0efb2917e80a82808a930d34af7f492"
V24_COORD_MANIFEST_PIN = (
    "9722c3e45d316ff58ed326df074eec138b81ceb520191a492088cc0ddef63762")
V24_RUN_MANIFEST_PIN = (
    "f602351f4844531d409cfc256939508103e3c2d4f57fe21d7f77c34494aa4036")
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
HTTP_UA = ("toscanini-c5-regen/1.0 "
           "(diagnostic evidence regeneration; sha256-verified pinned cohort)")

# ---------------------------------------------------------------------------
# Evaluator: one isolated python process per (engine checkout, structure).
# Identical ingestion/adjudication path for both engines; every attribute is
# getattr-guarded because the two revisions ship different result dataclasses
# (candidate: indeterminate_causes, no containment; old rev: containment).
# ---------------------------------------------------------------------------
_EVALUATOR = r'''
import json, sys, re

content_path, sha = sys.argv[1], sys.argv[2]

from tos.ingestion.processor import IngestionProcessor
from tos.engine.tier1_measurements import Tier1Measurements
from tos.engine.adjudicator import AdjudicationInput, adjudicate_laws

content = open(content_path, "rb").read()
structure = IngestionProcessor.run(content, "origin.pdb", "Audit", "Audit")
results, coverage, third = Tier1Measurements.run_full_audit(
    structure, user_intent="NONE")

_is_exp = bool(getattr(structure.confidence, "is_experimental", False))
_method_type = getattr(structure.confidence, "method", "predicted")
adj = adjudicate_laws(AdjudicationInput(
    raw_measurements=results,
    coverage=coverage,
    is_experimental=_is_exp,
    method_type=_method_type,
    architecture="NONE",
))

resolution = None
for _line in content.decode(errors="ignore").splitlines():
    if _line.startswith("REMARK   2 RESOLUTION."):
        _m = re.search(r"RESOLUTION\.\s+(\d+\.\d+)\s+ANGSTROM", _line)
        if _m:
            resolution = float(_m.group(1))
        break

out = {}
if isinstance(third, (int, float)) and not isinstance(third, bool):
    out["atom_coverage_return"] = third
laws = {}
if getattr(adj, "law_rows", None):
    for entry in adj.law_rows:
        entry = entry if isinstance(entry, dict) else {}
        law_id = entry.get("law_id") or entry.get("id")
        if law_id is None:
            continue
        laws[law_id] = {
            "status": entry.get("status"),
            "observed": entry.get("observed"),
            "sample_size": entry.get("sample_size",
                                     entry.get("sample")),
            "threshold": entry.get("threshold"),
            "method": entry.get("method"),
            "reason": entry.get("reason"),
        }
elif results:
    for law_id, row in results.items():
        laws[law_id] = {
            "status": row.get("status"),
            "observed": row.get("observed"),
            "sample_size": row.get("sample_size", row.get("sample")),
            "threshold": row.get("threshold"),
            "method": row.get("method"),
            "reason": row.get("reason"),
        }
det_total = len([1 for row in laws.values()
                 if row.get("method") == "deterministic"])
out.update({
    "atoms": len(structure.atoms),
    "sha256": sha,
    "bytes": len(content),
    "coverage": coverage,
    "det_passed": getattr(adj, "det_passed", None),
    "det_total": det_total,
    "heur_passed": getattr(adj, "heur_passed", None),
    "deterministic_score": getattr(adj, "deterministic_score", None),
    "advisory_score": getattr(adj, "advisory_score", None),
    "contained_law_failures": list(
        getattr(adj, "contained_law_failures", [])),
    "indeterminate_causes": list(getattr(adj, "indeterminate_causes", [])),
    "failing_deterministic": list(getattr(adj, "failing_deterministic", [])),
    "is_experimental": _is_exp,
    "method_type": _method_type,
    "resolution": resolution,
    "verdict": getattr(adj, "verdict", None),
    "laws": laws,
})
print(json.dumps(out, sort_keys=True))
'''


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git(root: str, *args: str) -> bytes:
    return subprocess.check_output(["git", *args], cwd=root)


def materialize(git_root: str, rev: str, dst: str) -> None:
    os.makedirs(dst, exist_ok=True)
    proc = subprocess.Popen(["git", "archive", "--format=tar", rev],
                            cwd=git_root, stdout=subprocess.PIPE)
    subprocess.run(["tar", "-xC", dst], stdin=proc.stdout, check=True)
    proc.stdout.close()
    if proc.wait() != 0:
        raise RuntimeError(f"git archive failed for {rev}")


def http_get(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": HTTP_UA})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def fetch_cached(url: str, cache_dir: str, expected_sha: str = None,
                 expected_bytes: int = None, write_payload: bool = True):
    """Idempotent sha-pinned fetch. Returns a record with the cache path."""
    out_path = os.path.join(
        cache_dir, hashlib.sha256(url.encode()).hexdigest() + ".download")
    if os.path.exists(out_path):
        digest = sha256_file(out_path)
        if expected_sha is None or digest == expected_sha:
            size = os.path.getsize(out_path)
            return {"sha256": digest, "bytes": size, "path": out_path,
                    "bytes_match": (expected_bytes is None
                                    or size == expected_bytes),
                    "sha_match": (expected_sha is None
                                  or digest == expected_sha)}
        os.remove(out_path)
    try:
        data = http_get(url)
    except urllib.error.HTTPError as exc:
        return {"sha256": None, "bytes": None, "path": None,
                "bytes_match": False, "sha_match": False,
                "detail": f"HTTPError: HTTP Error {exc.code}"}
    except Exception as exc:
        return {"sha256": None, "bytes": None, "path": None,
                "bytes_match": False, "sha_match": False,
                "detail": str(exc)[:160]}
    if write_payload:
        with open(out_path, "wb") as fh:
            fh.write(data)
    return {"sha256": sha256_bytes(data), "bytes": len(data),
            "path": out_path if write_payload else None,
            "bytes_match": (expected_bytes is None
                            or len(data) == expected_bytes),
            "sha_match": (expected_sha is None
                          or sha256_bytes(data) == expected_sha)}


def probe_url(url: str):
    try:
        http_get(url, 60)
        return 200
    except urllib.error.HTTPError as exc:
        return exc.code
    except Exception as exc:
        return str(exc)[:120]


# ---------------------------------------------------------------------------
# mmCIF -> PDB conversion (deterministic; applied identically to both engines)
# ---------------------------------------------------------------------------
def cif_to_pdb_bytes(cif_bytes: bytes) -> bytes:
    """Minimal deterministic mmCIF to PDB conversion for AlphaFold
    output files: parses the ``_atom_site`` loop only. Fixed-width PDB
    columns, occupancy 1.00, B-factor = mmCIF ``B_iso_or_equiv`` (the
    AlphaFoldDB pLDDT), chain 'A'."""
    text = cif_bytes.decode("utf-8", errors="replace")
    columns = []
    data_rows = []
    state = "other"
    for raw in text.splitlines():
        stripped = raw.strip()
        if stripped.startswith("loop_"):
            state = "scan"
            continue
        if stripped.startswith("_atom_site."):
            columns.append(stripped.split(".", 1)[1])
            state = "header"
            continue
        if state == "scan":
            continue
        if state == "header":
            if stripped.startswith("_") or not stripped:
                state = "other"
                continue
            if stripped.startswith("#"):
                state = "other"
                continue
            data_rows.append(stripped.split())
            state = "rows"
            continue
        if state == "rows":
            if stripped.startswith("#") or stripped.startswith("_"):
                break
            if stripped:
                data_rows.append(stripped.split())
    if not columns:
        raise RuntimeError("mmCIF _atom_site header not found")

    idx = {name: i for i, name in enumerate(columns)}
    need = ("label_atom_id", "label_comp_id", "Cartn_x", "Cartn_y",
            "Cartn_z", "B_iso_or_equiv", "type_symbol")
    for field in need:
        if field not in idx:
            raise RuntimeError(f"mmCIF column {field} missing")
    atom_i = idx["label_atom_id"]
    comp_i = idx["label_comp_id"]
    seq_i = idx.get("label_seq_id", idx.get("auth_seq_id"))
    if seq_i is None:
        raise RuntimeError("mmCIF sequence id column missing")
    x_i, y_i, z_i = idx["Cartn_x"], idx["Cartn_y"], idx["Cartn_z"]
    b_i, el_i = idx["B_iso_or_equiv"], idx["type_symbol"]

    pdb_lines = ["HEADER    CONVERTED FROM ALPHAFOLDB mmCIF"]
    serial = 1
    for row in data_rows:
        atom = row[atom_i].strip('"')
        comp = row[comp_i]
        seq = int(row[seq_i])
        x, y, z = float(row[x_i]), float(row[y_i]), float(row[z_i])
        b_iso = float(row[b_i])
        element = row[el_i].strip().upper()
        pdb_lines.append(
            f"ATOM  {serial:5d} {atom:>4s} {comp:3s} A{seq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
            f"{1.00:6.2f}{b_iso:6.2f}          {element:>2s}  ")
        serial += 1
    return ("\n".join(pdb_lines) + "\n").encode("ascii")


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------
def evaluate(engine_backend_dir: str, pdb_path: str, sha: str) -> dict:
    comp = subprocess.run(
        [sys.executable, "-c", _EVALUATOR, pdb_path, sha],
        cwd=engine_backend_dir, capture_output=True, text=True,
        timeout=1800)
    if comp.returncode != 0:
        raise RuntimeError(
            f"evaluator failed ({engine_backend_dir}):\n"
            f"{comp.stderr.strip()[-1200:]}")
    return json.loads(comp.stdout)


def engine_code_fingerprint(engine_backend_dir: str) -> str:
    """Cache-key fingerprint of the engine sources the evaluator depends on;
    any engine edit invalidates cached evaluations."""
    h = hashlib.sha256()
    for rel in ("tos/engine/tier1_measurements.py",
                "tos/engine/adjudicator.py",
                "tos/ingestion/processor.py"):
        path = os.path.join(engine_backend_dir, rel)
        h.update(rel.encode())
        h.update(sha256_file(path).encode())
    return h.hexdigest()[:12]


def evaluate_cached(engine_key: str, engine_backend_dir: str, pdb_path: str,
                    sha: str, cache_root: str) -> dict:
    """Idempotent per-(engine, target) evaluation cache keyed by target name
    and guarded by the input file's SHA-256 plus the engine code fingerprint;
    lets interrupted runs resume and invalidates stale engine results."""
    fingerprint = engine_code_fingerprint(engine_backend_dir)
    out_dir = os.path.join(cache_root, "evals", fingerprint)
    os.makedirs(out_dir, exist_ok=True)
    target = os.path.basename(pdb_path)
    if target.endswith(".converted.pdb"):
        target = target.split(".converted.pdb")[0]
    else:
        target = target + ".v24"
    cache_path = os.path.join(out_dir, target + ".json")
    if os.path.exists(cache_path):
        with open(cache_path, "r", encoding="utf-8") as fh:
            cached = json.load(fh)
        if cached.get("input_sha256") == sha:
            return cached["result"]
    result = evaluate(engine_backend_dir, pdb_path, sha)
    with open(cache_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump({"input_sha256": sha, "result": result}, fh,
                  indent=0, sort_keys=True)
        fh.write("\n")
    return result


def law_status_deltas(laws_old: dict, laws_new: dict):
    deltas = []
    for law_id in sorted(set(laws_old) | set(laws_new)):
        o = laws_old.get(law_id, {})
        n = laws_new.get(law_id, {})
        if (o.get("status") or "ABSENT") != (n.get("status") or "ABSENT"):
            deltas.append({
                "law": law_id,
                "old_status": o.get("status") or "ABSENT",
                "new_status": n.get("status") or "ABSENT",
                "observed_old": o.get("observed"),
                "observed_new": n.get("observed"),
                "sample_old": o.get("sample_size"),
                "sample_new": n.get("sample_size")})
    return deltas


def verdict_bucket(verdicts):
    order = ("PASS", "VETO", "INDETERMINATE")
    return " / ".join(f"{sum(1 for v in verdicts if v == k)} {k}"
                      for k in order)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--old-rev", default=OLD_REV)
    ap.add_argument("--base-rev", default=BASE_REV,
                    help="candidate branch base; primary comparison "
                         "base-vs-new isolates the reviewed round-25 delta")
    ap.add_argument("--new-rev", default=None, help="default: current HEAD")
    ap.add_argument("--v24-ref", default=DEFAULT_V24_REF)
    ap.add_argument("--work-root", required=True,
                    help="scratch dir (download cache + engine checkouts)")
    ap.add_argument("--out-root", default=None,
                    help="default: evidence/engine in this tree")
    ap.add_argument("--skip-fetch", action="store_true",
                    help="use the cached bytes only (no network)")
    ap.add_argument("--probe", action="store_true",
                    help="materialize + evaluate one small target per "
                         "engine, then exit (pipeline validation)")
    ap.add_argument("--phase", choices=["pinned", "support", "all", "emit"],
                    default="all",
                    help="work chunking: pinned (40 v24 evaluations), "
                         "support (AFDB v6 evaluations), all (default), "
                         "emit (assemble files + narrative only)")
    args = ap.parse_args()

    git_root = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"], check=True,
        capture_output=True, text=True).stdout.strip()
    old_rev = args.old_rev
    base_rev = args.base_rev
    new_rev = args.new_rev or git(git_root, "rev-parse", "HEAD").decode().strip()
    new_rev12 = new_rev[:12]

    work_root = os.path.abspath(args.work_root)
    out_root = os.path.abspath(args.out_root or os.path.join(
        git_root, "evidence", "engine"))
    c5_out = os.path.join(out_root, "c5")
    os.makedirs(c5_out, exist_ok=True)
    dl_dir = os.path.join(work_root, "downloads")
    os.makedirs(dl_dir, exist_ok=True)
    old_backend = os.path.join(work_root, "-old", "backend")
    base_backend = os.path.join(work_root, "-base", "backend")
    new_backend = os.path.join(work_root, "-new", "backend")
    materialize(git_root, old_rev, os.path.join(work_root, "-old"))
    materialize(git_root, base_rev, os.path.join(work_root, "-base"))
    materialize(git_root, new_rev, os.path.join(work_root, "-new"))

    # ---- v24 manifests (byte-pins verified against the immutable ref) -----
    coord_raw = git(git_root, "show", args.v24_ref + ":backend/calibration/runs/v24/coordinate_manifest.json")
    run_raw = git(git_root, "show", args.v24_ref + ":backend/calibration/runs/v24/run_manifest.json")
    coord_sha, run_sha = sha256_bytes(coord_raw), sha256_bytes(run_raw)
    if coord_sha != V24_COORD_MANIFEST_PIN or run_sha != V24_RUN_MANIFEST_PIN:
        raise SystemExit("v24 manifest pins do not match the immutable "
                         f"ref ({coord_sha} vs {V24_COORD_MANIFEST_PIN}; "
                         f"{run_sha} vs {V24_RUN_MANIFEST_PIN})")
    coord_manifest = json.loads(coord_raw)

    acquisition = {"v24_manifest_total": len(coord_manifest),
                   "v24_acquired_verified": 0, "v24_unavailable": 0,
                   "afdb_v6_attempted": 0, "afdb_v6_acquired": 0,
                   "afdb_v6_missing": []}
    verification = {
        "coordinate_manifest": {
            "path_in_v24_ref": "backend/calibration/runs/v24/coordinate_manifest.json",
            "sha256": coord_sha},
        "run_manifest": {
            "path_in_v24_ref": "backend/calibration/runs/v24/run_manifest.json",
            "sha256": run_sha},
        "targets": {}}

    # ---- v24 cohort: fetch + verify + evaluate ----------------------------
    if args.probe:
        targets_probe = sorted(k for k, m in coord_manifest.items()
                               if m["acquisition_status"] == "acquired")[0:1]
        meta = coord_manifest[targets_probe[0]]
        rec = fetch_cached(meta["source"], dl_dir,
                           expected_sha=meta["sha256"],
                           expected_bytes=meta["bytes"])
        old_out = evaluate(old_backend, rec["path"], rec["sha256"])
        base_out = evaluate(base_backend, rec["path"], rec["sha256"])
        new_out = evaluate(new_backend, rec["path"], rec["sha256"])
        print("PROBE old:", json.dumps(
            {k: old_out[k] for k in ("atoms", "coverage", "verdict",
                                     "det_passed", "det_total",
                                     "contained_law_failures",
                                     "indeterminate_causes")}, sort_keys=True))
        print("PROBE base:", json.dumps(
            {k: base_out[k] for k in ("atoms", "coverage", "verdict",
                                      "det_passed", "det_total",
                                      "contained_law_failures",
                                      "indeterminate_causes")}, sort_keys=True))
        print("PROBE new:", json.dumps(
            {k: new_out[k] for k in ("atoms", "coverage", "verdict",
                                     "det_passed", "det_total",
                                     "contained_law_failures",
                                     "indeterminate_causes")}, sort_keys=True))
        print("PROBE base LAw rows vs new:", json.dumps(
            {k: [base_out["laws"].get(k, {}).get("status"),
                 new_out["laws"].get(k, {}).get("status")]
             for k in sorted(set(base_out["laws"]) | set(new_out["laws"]))
             if base_out["laws"].get(k, {}).get("status")
             != new_out["laws"].get(k, {}).get("status")}, sort_keys=True))
        return 0

    state_path = os.path.join(work_root, "state.json")
    pinned_paths = {}
    supp_paths = {}
    af_records = {}
    if os.path.exists(state_path) and not args.probe:
        with open(state_path, "r", encoding="utf-8") as fh:
            state = json.load(fh)
        acquisition = state["acquisition"]
        verification = state["verification"]
        af_records = state["af_records"]
        for target, rec in state["pinned"].items():
            raw_path = os.path.join(
                dl_dir, hashlib.sha256(rec["url"].encode()).hexdigest()
                + ".download")
            pinned_paths[target] = {"path": raw_path,
                                    "sha256": rec["sha256"],
                                    "bytes": rec["bytes"]}
        for target, rec in state["supp"].items():
            converted_path = os.path.join(dl_dir, target + ".converted.pdb")
            supp_paths[target] = {"path": converted_path,
                                  "sha256": rec["convert_sha256"],
                                  "bytes": rec["convert_bytes"]}
    else:
        for target, meta in sorted(coord_manifest.items()):
            if meta["acquisition_status"] != "acquired":
                detail = (probe_url(meta["source"])
                          if not args.skip_fetch
                          else "skipped (--skip-fetch)")
                verification["targets"][target] = {
                    "manifest_status": meta["acquisition_status"],
                    "manifest_bytes": meta.get("bytes"),
                    "manifest_sha256": meta.get("sha256"),
                    "fetched_bytes": None, "fetched_sha256": None,
                    "bytes_match": False, "sha_match": False,
                    "probe_detail": detail}
                acquisition["v24_unavailable"] += 1
                continue
            rec = fetch_cached(meta["source"], dl_dir,
                               expected_sha=meta["sha256"],
                               expected_bytes=meta["bytes"])
            if not rec["sha_match"] or not rec["bytes_match"]:
                raise SystemExit(
                    f"{target}: v24 pinned bytes failed verification: {rec}")
            verification["targets"][target] = {
                "manifest_status": meta["acquisition_status"],
                "manifest_bytes": meta["bytes"],
                "manifest_sha256": meta["sha256"],
                "fetched_bytes": rec["bytes"],
                "fetched_sha256": rec["sha256"],
                "bytes_match": rec["bytes_match"],
                "sha_match": rec["sha_match"],
                "probe_detail": None}
            acquisition["v24_acquired_verified"] += 1
            pinned_paths[target] = rec

    def acquire_supp(target: str, accession: str) -> None:
        """Metadata-first AFDB v6 acquisition (delivered-record semantics);
        attempts are idempotent and update the acquisition counters."""
        try:
            metadata = json.loads(http_get(
                AFDB_API.format(accession=accession)))
            # AFDB serves some accessions under migrated numeric
            # entity ids; prefer the canonical AF-{acc}-F* id, fall
            # back to whatever served entity the metadata lists
            # (mirrors the delivered acquisition records).
            candidates = [e for e in metadata
                          if e.get("entryId", "").startswith(
                              f"AF-{accession}-F")]
            chosen = None
            for entry in candidates:
                if chosen is None or int(
                        entry.get("latestVersion", 0) or 0) > int(
                            chosen.get("latestVersion", 0) or 0):
                    chosen = entry
            if chosen is None:
                chosen = metadata[0]
            url = chosen.get("cifUrl") or chosen.get("pdbUrl")
            rec = fetch_cached(url, dl_dir)
            if rec["path"] is None:
                raise RuntimeError(rec.get("detail", "fetch failed"))
            entry_id = chosen["entryId"]
            migrated = entry_id if (
                not entry_id.startswith(f"AF-{accession}-")
            ) else None
            af_records[target] = {
                "fmt": "cif" if url.endswith(".cif") else "pdb",
                "latest_version": int(
                    chosen.get("latestVersion") or 1),
                "migrated_identifier": migrated,
                "model_entity_id": entry_id,
                "requested_entity_id": entry_id,
                "sha256": rec["sha256"],
                "source_url": url}
            converted = cif_to_pdb_bytes(open(rec["path"], "rb").read())
            convert_sha = sha256_bytes(converted)
            converted_path = os.path.join(dl_dir,
                                          target + ".converted.pdb")
            with open(converted_path, "wb") as fh:
                fh.write(converted)
            supp_paths[target] = {"path": converted_path,
                                  "sha256": convert_sha,
                                  "bytes": len(converted)}
            acquisition["afdb_v6_acquired"] += 1
            if target in acquisition["afdb_v6_missing"]:
                acquisition["afdb_v6_missing"].remove(target)
        except Exception as exc:
            if target not in acquisition["afdb_v6_missing"]:
                acquisition["afdb_v6_missing"].append(target)
            af_records[target] = {"model_entity_id": target,
                                  "unavailable": str(exc)[:160]}
        with open(os.path.join(c5_out, f"{target}.json"), "w",
                  encoding="utf-8", newline="\n") as fh:
            json.dump(af_records[target], fh, indent=1, sort_keys=True)
            fh.write("\n")

    resume_state = os.path.exists(state_path) and not args.probe
    if resume_state:
        # Only re-attempt the targets the previous run failed to acquire;
        # kept attempts are idempotent but not re-counted.
        retry = list(acquisition["afdb_v6_missing"])
        for target in retry:
            acquire_supp(target, target.replace("AF-", "").replace("-F1", ""))
    else:
        for target, meta in coord_manifest.items():
            if meta["acquisition_status"] == "acquired":
                continue
            accession = target.replace("AF-", "").replace("-F1", "")
            acquisition["afdb_v6_attempted"] += 1
            acquire_supp(target, accession)
        # drop supplemental records the retry could not acquire so the
        # evaluator never sees a stale file for an unavailable target
        supp_paths = {t: r for t, r in supp_paths.items()
                      if t not in acquisition["afdb_v6_missing"]}
        af_records = {t: r for t, r in af_records.items()
                      if t not in acquisition["afdb_v6_missing"]}

    state = {
        "acquisition": acquisition, "verification": verification,
        "af_records": af_records,
        "pinned": {t: {"url": coord_manifest[t]["source"],
                       "sha256": r["sha256"], "bytes": r["bytes"]}
                   for t, r in pinned_paths.items()},
        "supp": {t: {"url": af_records[t].get("source_url"),
                     "convert_sha256": r["sha256"],
                     "convert_bytes": r["bytes"]}
                 for t, r in supp_paths.items()}}
    with open(state_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(state, fh, indent=1, sort_keys=True)
        fh.write("\n")

    gov_cache = args.work_root

    evals = {"old": {}, "base": {}, "new": {}}
    supp_evals = {"old": {}, "base": {}, "new": {}}
    if args.phase in ("pinned", "all"):
        for target, rec in sorted(pinned_paths.items()):
            evals["old"][target] = evaluate_cached(
                "old", old_backend, rec["path"], rec["sha256"], gov_cache)
            evals["base"][target] = evaluate_cached(
                "base", base_backend, rec["path"], rec["sha256"], gov_cache)
            evals["new"][target] = evaluate_cached(
                "new", new_backend, rec["path"], rec["sha256"], gov_cache)
            print(f"evaluated pinned {target}", flush=True)
    else:
        for target in sorted(pinned_paths):
            for engine, backend in (("old", old_backend),
                                    ("base", base_backend),
                                    ("new", new_backend)):
                evals[engine][target] = evaluate_cached(
                    engine, backend, pinned_paths[target]["path"],
                    pinned_paths[target]["sha256"], gov_cache)
    if args.phase in ("support", "all"):
        for target, rec in sorted(supp_paths.items()):
            supp_evals["old"][target] = evaluate_cached(
                "old", old_backend, rec["path"], rec["sha256"], gov_cache)
            supp_evals["base"][target] = evaluate_cached(
                "base", base_backend, rec["path"], rec["sha256"], gov_cache)
            supp_evals["new"][target] = evaluate_cached(
                "new", new_backend, rec["path"], rec["sha256"], gov_cache)
            print(f"evaluated supplementary {target}", flush=True)
    else:
        for target in sorted(supp_paths):
            for engine, backend in (("old", old_backend),
                                    ("base", base_backend),
                                    ("new", new_backend)):
                supp_evals[engine][target] = evaluate_cached(
                    engine, backend,
                    supp_paths[target]["path"],
                    supp_paths[target]["sha256"], gov_cache)

    if args.phase in ("pinned", "support"):
        missing = ("pinned" if (args.phase == "pinned" and
                                len(evals["new"]) != len(pinned_paths)) else
                   "support" if (args.phase == "support" and
                                 len(supp_evals["new"]) != len(supp_paths))
                   else "done")
        print(f"PHASE {args.phase}: {missing}", flush=True)
        if missing == "done":
            args.phase = "emit"
        else:
            return 0

    # ---- cohort files ------------------------------------------------------
    def cohort(label_desc, evs, engine_desc):
        targets = {}
        for target, ev in sorted(evs.items()):
            entry = dict(ev)
            targets[target] = entry
        return {"backend": engine_desc, "label": label_desc,
                "targets": targets}

    old_fl = cohort("old-c78", evals["old"],
                    f"git materialization of engine rev {old_rev[:12]}")
    base_fl = cohort(f"base-{base_rev[:12]}", evals["base"],
                     f"git materialization of engine rev {base_rev[:12]}")
    new_fl = cohort(f"new-candidate-{new_rev12}", evals["new"],
                    f"git materialization of engine rev {new_rev[:12]}")
    old_sf = cohort("old-c78-af-v6", supp_evals["old"],
                    f"git materialization of engine rev {old_rev[:12]}")
    base_sf = cohort(f"base-{base_rev[:12]}-af-v6", supp_evals["base"],
                     f"git materialization of engine rev {base_rev[:12]}")
    new_sf = cohort(f"new-candidate-{new_rev12}-af-v6",
                    supp_evals["new"],
                    f"git materialization of engine rev {new_rev[:12]}")
    for name, payload in (("v24-pinned-old-c78.json", old_fl),
                          (f"v24-pinned-base-{base_rev[:12]}.json", base_fl),
                          (f"v24-pinned-new-{new_rev12}.json", new_fl),
                          ("af-v6-supp-old-c78.json", old_sf),
                          (f"af-v6-supp-base-{base_rev[:12]}.json", base_sf),
                          (f"af-v6-supp-new-{new_rev12}.json", new_sf)):
        with open(os.path.join(c5_out, name), "w", encoding="utf-8",
                  newline="\n") as fh:
            json.dump(payload, fh, indent=1, sort_keys=True)
            fh.write("\n")
    with open(os.path.join(c5_out, "v24-coordinate-verification.json"),
              "w", encoding="utf-8", newline="\n") as fh:
        json.dump(verification, fh, indent=1, sort_keys=True)
        fh.write("\n")

    # ---- narrative ---------------------------------------------------------
    def pair_rows(evs_old, evs_new, sorted_targets):
        rows = []
        for target in sorted_targets:
            o, n = evs_old.get(target), evs_new.get(target)
            deltas = law_status_deltas(o["laws"], n["laws"])
            delta_cell = "none" if not deltas else "; ".join(
                f"{d['law']}: {d['old_status']} -> {d['new_status']}"
                for d in deltas)
            rows.append((target, o, n, deltas, delta_cell))
        return rows

    pinned_rows = pair_rows(evals["old"], evals["new"], pinned_paths)
    supp_rows = pair_rows(supp_evals["old"], supp_evals["new"], supp_paths)

    verdict_changes_pinned = [r for r in pinned_rows if r[1]["verdict"] != r[2]["verdict"]]
    verdict_changes_supp = [r for r in supp_rows if r[1]["verdict"] != r[2]["verdict"]]

    law_changes_pinned = [(t, deltas, engine_label)
                          for (t, _o, _n, deltas, _c) in pinned_rows
                          for deltas_ignored in [None]
                          for engine_label in [None]
                          if deltas]
    # flatten (target law-status-change rows)
    law_change_rows = []
    for (t, o, n, deltas, _c) in pinned_rows:
        for d in deltas:
            law_change_rows.append((t, d))
    for (t, o, n, deltas, _c) in supp_rows:
        for d in deltas:
            law_change_rows.append((t, d))

    lines = []
    lines.append("# C5 — v24 verdict-change audit (regenerated for the candidate)")
    lines.append("")
    lines.append("**Status:** diagnostic evidence for owner review. Not a threshold validation.")
    lines.append(f"**Regenerated by:** "
                 f"`docs/audit/2026-10-22_od_port_candidate/scripts/regenerate_c5_evidence.py`")
    lines.append(f"**Base engine (candidate branch base, pre-port):** `{base_rev}`")
    lines.append(f"**New engine (candidate review branch tip):** `{new_rev}`")
    lines.append(f"**Cross-lineage reference engine (od-scope lineage):** `{old_rev}`")
    lines.append(f"**Coordinate manifest pin (verified):** {V24_COORD_MANIFEST_PIN}")
    lines.append(f"**Run manifest pin (verified):** {V24_RUN_MANIFEST_PIN}")
    lines.append("**Primary comparison rule:** base-vs-new. The candidate "
                 "branch differs from its base ONLY by the reviewed "
                 "round-25 port, so the base-vs-new delta attributes each "
                 "verdict/law change to the reviewed OD-03/04/20 package.")
    lines.append("**Cross-lineage comparison rule:** the od-scope-lineage "
                 "engine (old rev) additionally carries the earlier P0/WP-SCI-1 "
                 "measurement-definition items (file-order residue pool, "
                 "conformer selection, coverage denominator, honest empty "
                 "states) that were NOT part of the reviewed round-25 delta; "
                 "that comparison documents lineage differences only.")
    lines.append("")
    lines.append("## Provenance (portable)")
    lines.append("")
    lines.append("- Engine checkouts are materialized byte-exactly from git "
                 "objects via `git archive`; every evaluation runs in an "
                 "isolated subprocess, so verdict deltas come only from "
                 "engine differences.")
    lines.append("- All v24 coordinate bytes are re-fetched from the pinned "
                 "URLs and verified against the pinned byte counts and "
                 "SHA-256 values (`v24-coordinate-verification.json`). No "
                 "coordinate bytes were substituted or fabricated.")
    lines.append("- Supplementary AlphaFold models are fetched metadata-first "
                 "from the AlphaFold DB API; mmCIF payloads are converted to "
                 "PDB records by a deterministic converter inside the script, "
                 "applied identically to every engine.")
    lines.append("- No machine-specific path is recorded in any artifact.")
    lines.append("")
    lines.append("## Acquisition summary")
    lines.append("")
    lines.append(f"- v24 manifest: **{acquisition['v24_manifest_total']} targets**.")
    lines.append(f"- v24 acquired and SHA-256 verified: "
                 f"**{acquisition['v24_acquired_verified']}/{acquisition['v24_manifest_total']}**.")
    lines.append(f"- v24 acquisition failures: "
                 f"**{acquisition['v24_unavailable']}/{acquisition['v24_manifest_total']}** "
                 f"(AlphaFold `model_v4` URLs; re-probed during regeneration).")
    lines.append(f"- Supplementary current-AFDB acquisition: "
                 f"**{acquisition['afdb_v6_acquired']}/{acquisition['afdb_v6_attempted']}**"
                 + (f" unavailable: {', '.join(acquisition['afdb_v6_missing'])}."
                    if acquisition['afdb_v6_missing'] else "."))
    lines.append("")

    # ---- primary attribution comparison: base vs new ----------------------
    base_rows = pair_rows(evals["base"], evals["new"], pinned_paths)
    base_supp_rows = pair_rows(supp_evals["base"], supp_evals["new"],
                               supp_paths)
    vc_base = [r for r in base_rows if r[1]["verdict"] != r[2]["verdict"]]
    vc_base_supp = [r for r in base_supp_rows
                    if r[1]["verdict"] != r[2]["verdict"]]
    law_row_changes_base = []
    for (t, o, n, deltas, _c) in base_rows:
        for d in deltas:
            law_row_changes_base.append((t, d))
    for (t, o, n, deltas, _c) in base_supp_rows:
        for d in deltas:
            law_row_changes_base.append((t, d))

    lines.append("## ROUND-25 ATTRIBUTION: base vs candidate (primary)")
    lines.append("")
    lines.append(f"### v24-pinned cohort ({len(pinned_paths)}/{acquisition['v24_manifest_total']})")
    lines.append("")
    lines.append("Base: " + verdict_bucket([r[1]["verdict"] for r in base_rows]) + ".")
    lines.append("New: " + verdict_bucket([r[2]["verdict"] for r in base_rows]) + ".")
    lines.append("")
    lines.append(f"**Verdict changes: {len(vc_base)}/{len(pinned_paths)}.**")
    lines.append("")
    lines.append("| Target | Base verdict | New verdict | Changed law rows |")
    lines.append("|---|---|---|---|")
    for (t, o, n, deltas, cell) in base_rows:
        if o["verdict"] != n["verdict"] or deltas:
            lines.append(f"| {t} | {o['verdict']} | {n['verdict']} | {cell} |")
    lines.append("")
    lines.append(f"### Supplementary predicted cohort ({len(supp_paths)}/{acquisition['afdb_v6_attempted']})")
    lines.append("")
    lines.append(f"These rows are supplementary (current AlphaFold DB, not "
                 f"v24-pinned) and must not be reported as v24-pinned evidence.")
    lines.append("")
    lines.append("Base: " + verdict_bucket([r[1]["verdict"] for r in base_supp_rows]) + ".")
    lines.append("New: " + verdict_bucket([r[2]["verdict"] for r in base_supp_rows]) + ".")
    lines.append("")
    lines.append(f"**Verdict changes: {len(vc_base_supp)}/{len(supp_paths)}.**")
    lines.append("")
    if base_supp_rows:
        lines.append("| Target | Base verdict | New verdict | Changed law rows |")
        lines.append("|---|---|---|---|")
        for (t, o, n, deltas, cell) in base_supp_rows:
            if o["verdict"] != n["verdict"] or deltas:
                lines.append(f"| {t} | {o['verdict']} | {n['verdict']} | {cell} |")
        lines.append("")
    total_base = len(vc_base) + len(vc_base_supp)
    lines.append(f"**{total_base} verdict changes across "
                 f"{len(pinned_paths) + len(supp_paths)} measured structures** "
                 f"({len(pinned_paths)} v24-pinned + {len(supp_paths)} supplementary), "
                 f"attributable to the reviewed round-25 port.")
    lines.append("")
    lines.append(f"{len(law_row_changes_base)} law rows changed status "
                 f"(base vs candidate, both cohorts; per-row table in "
                 f"`v24-pinned-base-*.json` / `af-v6-supp-base-*.json` against "
                 f"the `new-*` files):")
    lines.append("")
    lines.append("| Target | Law | Base status | New status | Observed base -> new | Sample base -> new |")
    lines.append("|---|---|---|---|---|---|")
    for (t, d) in law_row_changes_base:
        obs = f"{d['observed_old']} -> {d['observed_new']}"
        smp = f"{d['sample_old']} -> {d['sample_new']}"
        lines.append(f"| {t} | {d['law']} | {d['old_status']} | {d['new_status']} | {obs} | {smp} |")
    lines.append("")

    # ---- cross-lineage documentary comparison: c78 lineage vs candidate ---
    pinned_rows = pair_rows(evals["old"], evals["new"], pinned_paths)
    supp_rows = pair_rows(supp_evals["old"], supp_evals["new"], supp_paths)
    verdict_changes_pinned = [r for r in pinned_rows if r[1]["verdict"] != r[2]["verdict"]]
    verdict_changes_supp = [r for r in supp_rows if r[1]["verdict"] != r[2]["verdict"]]
    law_change_rows = []
    for (t, o, n, deltas, _c) in pinned_rows:
        for d in deltas:
            law_change_rows.append((t, d))
    for (t, o, n, deltas, _c) in supp_rows:
        for d in deltas:
            law_change_rows.append((t, d))

    lines.append("## CROSS-LINEAGE documentation (od-scope lineage vs candidate)")
    lines.append("")
    lines.append("The od-scope-lineage engine additionally carries the earlier "
                 "P0/WP-SCI-1 measurement-definition items; the candidate is "
                 "main + the reviewed port only. Differences below therefore "
                 "document lineage divergence, NOT the effect of the "
                 "reviewed round-25 package. In particular: `1PPE` differences "
                 "are the documented P0 item-4 symptom (deprecated (chain, seq, "
                 "icode) sort re-pairs 1PPE chain E around 183/184/184A/185 and "
                 "fabricates backbone gaps); `5T4O` differences come from "
                 "main-lineage LAW-150 empty-sample and LAW-195 "
                 "not-applicable semantics, both outside the ported scope.")
    lines.append("")
    lines.append("Base cross-check: the od-scope-lineage engine regenerates the "
                 "delivered `old-c78` artifacts (e.g. 1AKE atoms 3312, det "
                 "9/10; LAW-182 observed values identical to the frozen record).")
    lines.append("")
    lines.append(f"### v24-pinned cohort ({len(pinned_paths)}/{acquisition['v24_manifest_total']})")
    lines.append("")
    lines.append("Old (c78): " + verdict_bucket([r[1]["verdict"] for r in pinned_rows]) + ".")
    lines.append("New: " + verdict_bucket([r[2]["verdict"] for r in pinned_rows]) + ".")
    lines.append("")
    lines.append(f"**Verdict changes: {len(verdict_changes_pinned)}/{len(pinned_paths)}.**")
    lines.append("")
    lines.append("| Target | Old verdict | New verdict | Changed law rows |")
    lines.append("|---|---|---|---|")
    for (t, o, n, deltas, cell) in pinned_rows:
        if o["verdict"] != n["verdict"] or deltas:
            lines.append(f"| {t} | {o['verdict']} | {n['verdict']} | {cell} |")
    lines.append("")
    lines.append(f"### Supplementary predicted cohort ({len(supp_paths)}/{acquisition['afdb_v6_attempted']})")
    lines.append("")
    lines.append("Old (c78): " + verdict_bucket([r[1]["verdict"] for r in supp_rows]) + ".")
    lines.append("New: " + verdict_bucket([r[2]["verdict"] for r in supp_rows]) + ".")
    lines.append("")
    lines.append(f"**Verdict changes: {len(verdict_changes_supp)}/{len(supp_paths)}.**")
    lines.append("")
    if supp_rows:
        lines.append("| Target | Old verdict | New verdict | Old coverage | New coverage | Old atom coverage | New atom coverage |")
        lines.append("|---|---|---|---|---|---|---|")
        for (t, o, n, deltas, cell) in supp_rows:
            lines.append(
                f"| {t} | {o['verdict']} | {n['verdict']} | "
                f"{o['coverage']:.2f} | {n['coverage']:.2f} | "
                f"{(o.get('atom_coverage_return') if isinstance(o.get('atom_coverage_return'), (int, float)) else 'n/a')} | "
                f"{(n.get('atom_coverage_return') if isinstance(n.get('atom_coverage_return'), (int, float)) else 'n/a')} |")
        lines.append("")
    total_changes = len(verdict_changes_pinned) + len(verdict_changes_supp)
    lines.append(f"**{total_changes} cross-lineage verdict changes across "
                 f"{len(pinned_paths) + len(supp_paths)} measured structures** "
                 f"({len(pinned_paths)} v24-pinned + {len(supp_paths)} supplementary).")
    lines.append("")
    lines.append(f"{len(law_change_rows)} cross-lineage law rows changed status:")
    lines.append("")
    lines.append("| Target | Law | Old status | New status | Observed old -> new | Sample old -> new |")
    lines.append("|---|---|---|---|---|---|")
    for (t, d) in law_change_rows:
        obs = f"{d['observed_old']} -> {d['observed_new']}"
        smp = f"{d['sample_old']} -> {d['sample_new']}"
        lines.append(f"| {t} | {d['law']} | {d['old_status']} | {d['new_status']} | {obs} | {smp} |")
    lines.append("")
    if base_supp_rows:
        lines.append("## Predicted-cohort measurement scopes (supplementary, attributed comparison)")
        lines.append("")
        lines.append("| Target | Coverage base -> new | LAW-110 sample base -> new | LAW-130 sample base -> new | LAW-155 observed/sample base -> new | LAW-200 observed/sample base -> new |")
        lines.append("|---|---|---|---|---|---|")
        for (t, o, n, deltas, _c) in base_supp_rows:
            def cell(law_id, field):
                oo = o["laws"].get(law_id, {}).get(field)
                nn = n["laws"].get(law_id, {}).get(field)
                return f"{oo} -> {nn}" if (oo is not None or nn is not None) else "ABSENT"
            lines.append(
                f"| {t} | {o['coverage']:.2f} -> {n['coverage']:.2f} | "
                f"{cell('LAW-110','sample_size')} | {cell('LAW-130','sample_size')} | "
                f"{cell('LAW-155','observed')} / {cell('LAW-155','sample_size')} | "
                f"{cell('LAW-200','observed')} / {cell('LAW-200','sample_size')} |")
        lines.append("")
    lines.append("## Interpretation constraints (unchanged from the frozen record)")
    lines.append("")
    lines.append("- Diagnostic-only; no validation of v24 calibration, "
                 "thresholds, cohorts, splits or the historical verdict-change "
                 "table (OD-9, OD-14, OD-16).")
    lines.append("- `contained_law_failures` is an old-engine artifact; the "
                 "candidate intentionally has no containment policy, so its "
                 "rows always show `[]`.")
    lines.append("- `det_total` is derived per evaluation from rows whose "
                 "`method` label is `deterministic`, identically in both engines.")
    lines.append("- `LAW-110` sample semantics differ by design (old engine: "
                 "positional set; candidate: `pairs_checked`, the pairs "
                 "actually measured); this is verdict-neutral and documented "
                 "in the candidate's registry row.")
    lines.append("- `atom_coverage_return` is candidate-side only (the "
                 "candidate's run_full_audit third slot); experimental "
                 f"structures stay at 100.0.")
    lines.append("")
    narrative = "\n".join(lines)
    with open(os.path.join(out_root, f"C5-verdict-change-{new_rev12}.md"),
              "w", encoding="utf-8", newline="\n") as fh:
        fh.write(narrative)

    print("WROTE:", out_root)
    print("ATTRIBUTION verdict changes (base vs new): pinned",
          len(vc_base), "supp:", len(vc_base_supp))
    print("cross-lineage verdict changes (c78 vs new): pinned",
          len(verdict_changes_pinned), "supp:", len(verdict_changes_supp))
    print("law status changes: attributed",
          len(law_row_changes_base), "cross-lineage", len(law_change_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
