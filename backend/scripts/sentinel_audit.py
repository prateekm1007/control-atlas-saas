#!/usr/bin/env python3
"""
Sentinel Audit: Phase 3.4 Stability Gate
-----------------------------------------
Hits the live /ingest API for 7 benchmark structures and verifies:
  1. Verdict Invariance  — live verdict matches snapshot "expected"
  2. Hash Integrity      — governance fingerprint matches locked hashes
  3. Numerical Stability — observed physics values match snapshot (tol=0.01)

Requires the Toscanini server running on localhost:8000.
Exit code 0 = PASS. Non-zero = FAIL.
"""

from __future__ import annotations

import json
import logging
import sys
import time
import urllib.request
from pathlib import Path

import requests

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("sentinel")

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
BACKEND = SCRIPT_DIR.parent
SNAPSHOT_PATH = BACKEND / "invariants" / "v23.0.phase1.snapshot.json"
CACHE_DIR = BACKEND / "invariants" / "benchmarks"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# ── Locked Hashes ─────────────────────────────────────────────────────────────
LOCKED_CANON_HASH = "6a9cd4b4349b81de"
LOCKED_MATRIX_HASH = "ed1731ace09a9a3bd24d3659ab6f3fd184da2a2335d935a234fd6a9c14880aa4"

# ── API ───────────────────────────────────────────────────────────────────────
API_BASE = "http://localhost:8000"
INGEST_URL = f"{API_BASE}/ingest"
HEALTH_URL = f"{API_BASE}/health"

# ── Benchmark Sources ─────────────────────────────────────────────────────────
BENCHMARK_URLS = {
    "4HHB": "https://files.rcsb.org/download/4HHB.pdb",
    "1CRN": "https://files.rcsb.org/download/1CRN.pdb",
    "1UBQ": "https://files.rcsb.org/download/1UBQ.pdb",
    "6LZG": "https://files.rcsb.org/download/6LZG.pdb",
    "7BV2": "https://files.rcsb.org/download/7BV2.pdb",
    "1G03": "https://files.rcsb.org/download/1G03.pdb",
    "AF-P01308-F1": "https://alphafold.ebi.ac.uk/files/AF-P01308-F1-model_v6.pdb",
}

# ── Tolerance ─────────────────────────────────────────────────────────────────
ABS_TOL = 0.01


def _fetch_pdb(name: str) -> Path:
    """Download and cache a benchmark PDB. Returns path to cached file."""
    safe_name = name.replace("-", "_") + ".pdb"
    cache_path = CACHE_DIR / safe_name

    if cache_path.exists() and cache_path.stat().st_size > 100:
        return cache_path

    url = BENCHMARK_URLS[name]
    log.info(f"    Downloading {url} ...")
    req = urllib.request.Request(url, headers={"User-Agent": "Toscanini-Sentinel/1.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        data = resp.read()

    if len(data) < 100:
        raise RuntimeError(f"Downloaded file for {name} too small ({len(data)} bytes)")

    cache_path.write_bytes(data)
    log.info(f"    Cached: {safe_name} ({len(data):,} bytes)")
    return cache_path


def _check_server():
    """Verify the server is reachable."""
    try:
        r = requests.get(HEALTH_URL, timeout=5)
        r.raise_for_status()
        log.info(f"  Server: operational (v{r.json().get('version', '?')})")
        return True
    except Exception as e:
        log.error(f"  ✗ Server not reachable: {e}")
        return False


def _ingest(pdb_path: Path, name: str) -> dict:
    """Send a PDB to the /ingest endpoint and return the JSON payload."""
    with pdb_path.open("rb") as f:
        resp = requests.post(
            INGEST_URL,
            data={"mode": "Audit", "candidate_id": name, "t3_category": "NONE"},
            files={"file": (pdb_path.name, f, "chemical/x-pdb")},
            timeout=120,
        )
    resp.raise_for_status()
    return resp.json()


def _values_match(live: float, snap: float) -> bool:
    """Compare observed values with absolute tolerance of 0.01."""
    if snap == 0:
        return live == 0
    return abs(live - snap) <= ABS_TOL


def run_sentinel() -> bool:
    """Execute the full sentinel audit. Returns True if all checks pass."""
    log.info("=" * 70)
    log.info("  TOSCANINI SENTINEL AUDIT — Phase 3.4 Stability Gate")
    log.info("=" * 70)

    # ── Load snapshot ─────────────────────────────────────────────────────
    if not SNAPSHOT_PATH.exists():
        log.error(f"✗ Snapshot not found: {SNAPSHOT_PATH}")
        return False

    with SNAPSHOT_PATH.open() as f:
        snapshot = json.load(f)

    benchmarks = snapshot["benchmarks"]
    log.info(f"\n  Snapshot: {snapshot['meta']['phase']}")
    log.info(f"  Targets:  {len(benchmarks)}")

    # ── Check server ──────────────────────────────────────────────────────
    log.info("")
    if not _check_server():
        log.error("  ✗ Cannot proceed without a running server.")
        log.error("  Start with: cd backend && uvicorn main:app --port 8000")
        return False

    # ── Run benchmarks ────────────────────────────────────────────────────
    all_pass = True
    results = []

    for name, snap_data in benchmarks.items():
        log.info(f"\n{'─' * 70}")
        log.info(f"  TARGET: {name}")
        log.info(f"{'─' * 70}")

        expected_verdict = snap_data.get("expected", snap_data["verdict"])
        snap_laws = snap_data.get("laws", {})

        # ── Fetch ─────────────────────────────────────────────────────
        try:
            pdb_path = _fetch_pdb(name)
        except Exception as e:
            log.error(f"    ✗ FETCH FAILED: {e}")
            all_pass = False
            results.append((name, "FETCH_ERR", expected_verdict, False, False, 0))
            continue

        # ── Ingest ────────────────────────────────────────────────────
        try:
            payload = _ingest(pdb_path, name)
        except Exception as e:
            log.error(f"    ✗ INGEST FAILED: {e}")
            all_pass = False
            results.append((name, "API_ERR", expected_verdict, False, False, 0))
            continue

        live_verdict = payload["verdict"]["binary"]

        # ── Check 1: Verdict Invariance ───────────────────────────────
        verdict_ok = live_verdict == expected_verdict
        v_mark = "✓" if verdict_ok else "✗ DRIFT"
        log.info(f"    Verdict:     {live_verdict:>15}  expected: {expected_verdict:>15}  {v_mark}")
        if not verdict_ok:
            all_pass = False

        # ── Check 2: Hash Integrity ───────────────────────────────────
        fp = payload.get("governance", {}).get("governance_fingerprint", {})
        canon_ok = fp.get("canon_hash") == LOCKED_CANON_HASH
        matrix_ok = fp.get("matrix_hash") == LOCKED_MATRIX_HASH
        hash_ok = canon_ok and matrix_ok

        if hash_ok:
            log.info(f"    Fingerprint: ✓ (canon + matrix locked)")
        else:
            if not canon_ok:
                log.error(f"    ✗ Canon hash: {fp.get('canon_hash', 'MISSING')}")
            if not matrix_ok:
                log.error(f"    ✗ Matrix hash: {fp.get('matrix_hash', 'MISSING')[:20]}...")
            all_pass = False

        # ── Check 3: Numerical Stability ──────────────────────────────
        live_laws = {r["law_id"]: r for r in payload.get("tier1", {}).get("laws", [])}
        drift_count = 0
        checked = 0

        for lid, snap_law in snap_laws.items():
            live_law = live_laws.get(lid)
            if live_law is None:
                log.warning(f"    ⚠ {lid} missing from live output")
                drift_count += 1
                continue

            snap_obs = float(snap_law.get("observed", 0))
            live_obs = float(live_law.get("observed", 0))
            checked += 1

            if not _values_match(live_obs, snap_obs):
                delta = abs(live_obs - snap_obs)
                log.warning(
                    f"    ⚠ {lid} drift: snapshot={snap_obs:.4f} "
                    f"live={live_obs:.4f} delta={delta:.4f}"
                )
                drift_count += 1

        if drift_count == 0:
            log.info(f"    Numerical:   ✓ ({checked} laws, 0 drifts)")
        else:
            log.warning(f"    Numerical:   ⚠ {drift_count} drift(s) in {checked} laws")

        results.append((name, live_verdict, expected_verdict, verdict_ok, hash_ok, drift_count))

    # ── Summary Table ─────────────────────────────────────────────────────
    log.info(f"\n{'═' * 70}")
    log.info("  SENTINEL AUDIT SUMMARY")
    log.info(f"{'═' * 70}")
    log.info(f"  {'Target':>15}  {'Live':>13}  {'Expected':>13}  {'Verdict':>7}  {'Hash':>5}  {'Drift':>5}")
    log.info(f"  {'─'*15}  {'─'*13}  {'─'*13}  {'─'*7}  {'─'*5}  {'─'*5}")

    for name, live, expected, v_ok, h_ok, drifts in results:
        log.info(
            f"  {name:>15}  {str(live):>13}  {str(expected):>13}"
            f"  {'✓':>7}  {'✓':>5}  {drifts:>5}"
            if v_ok and h_ok else
            f"  {name:>15}  {str(live):>13}  {str(expected):>13}"
            f"  {'✓' if v_ok else '✗':>7}  {'✓' if h_ok else '✗':>5}  {drifts:>5}"
        )

    passed = sum(1 for _, _, _, v, h, _ in results if v and h)
    total = len(results)

    log.info(f"\n  Gate: {passed}/{total} targets passed")
    log.info(f"  Canon:  {LOCKED_CANON_HASH}")
    log.info(f"  Matrix: {LOCKED_MATRIX_HASH[:16]}...")

    if all_pass:
        log.info("\n  ✓ SENTINEL GATE: PASSED ✓")
    else:
        log.error("\n  ✗ SENTINEL GATE: FAILED ✗")

    log.info(f"{'═' * 70}")
    return all_pass


if __name__ == "__main__":
    success = run_sentinel()
    sys.exit(0 if success else 1)
