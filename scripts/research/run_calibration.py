"""
TOSCANINI Phase 2 Calibration Runner

Measures all 15 laws across the 50-structure calibration dataset.
Saves results to backend/calibration/calibration_results.json.

Usage: python3 run_calibration.py
Requires: Backend running on localhost:8000
"""
import requests
import json
import sys
import time
import os
from datetime import datetime, timezone

sys.path.insert(0, "backend")
from calibration.dataset_v23 import get_all_targets, CALIBRATION_DATASET

API_URL = "http://localhost:8000/ingest"
OUTPUT_PATH = "backend/calibration/calibration_results.json"

# Check backend is running
try:
    health = requests.get("http://localhost:8000/health", timeout=5)
    print(f"Backend: {health.json()}")
except Exception:
    print("ERROR: Backend not running on localhost:8000")
    sys.exit(1)

targets = get_all_targets()
print(f"\nRunning calibration on {len(targets)} structures...")
print(f"{'#':<4} {'PDB':<18} {'Category':<18} {'Verdict':<15} {'LAW-120':>8} {'LAW-160':>8} {'LAW-150':>8} {'Status'}")
print("=" * 105)

results = {
    "meta": {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "phase": "2-calibration",
        "total_targets": len(targets),
        "backend_version": health.json().get("version", "unknown"),
    },
    "measurements": {}
}

success_count = 0
error_count = 0

for idx, (pdb_id, category, resolution, desc) in enumerate(targets, 1):
    try:
        resp = requests.post(API_URL, data={
            "mode": "Benchmark",
            "candidate_id": pdb_id,
            "t3_category": "NONE"
        }, timeout=120)

        if resp.status_code != 200:
            print(f"{idx:<4} {pdb_id:<18} {category:<18} {'HTTP-' + str(resp.status_code):<15} {'':>8} {'':>8} {'':>8} ❌")
            error_count += 1
            continue

        d = resp.json()
        verdict = d["verdict"]["binary"]
        laws = {l["law_id"]: l for l in d["tier1"]["laws"]}

        law120 = laws.get("LAW-120", {}).get("observed", "—")
        law160 = laws.get("LAW-160", {}).get("observed", "—")
        law150 = laws.get("LAW-150", {}).get("observed", "—")

        results["measurements"][pdb_id] = {
            "category": category,
            "resolution": resolution,
            "desc": desc,
            "verdict": verdict,
            "coverage": d["verdict"]["coverage_pct"],
            "laws": {l["law_id"]: {
                "observed": l["observed"],
                "threshold": l["threshold"],
                "status": l["status"],
                "method": l["method"],
                "sample_size": l["sample_size"]
            } for l in d["tier1"]["laws"]}
        }

        flag = "✅" if verdict in ("PASS", "INDETERMINATE") else "⚠️"
        print(f"{idx:<4} {pdb_id:<18} {category:<18} {verdict:<15} {str(law120):>8} {str(law160):>8} {str(law150):>8} {flag}")
        success_count += 1

        # Brief pause to avoid hammering RCSB
        time.sleep(0.5)

    except Exception as e:
        print(f"{idx:<4} {pdb_id:<18} {category:<18} {'ERROR':<15} {'':>8} {'':>8} {'':>8} ❌ {str(e)[:30]}")
        error_count += 1

print("=" * 105)
print(f"\nCompleted: {success_count} success, {error_count} errors")

# Save results
os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
with open(OUTPUT_PATH, "w") as f:
    json.dump(results, f, indent=4)

print(f"Results saved to: {OUTPUT_PATH}")

# Quick LAW-120 summary by category
print("\n" + "=" * 60)
print("LAW-120 SUMMARY BY CATEGORY")
print("=" * 60)
for category in CALIBRATION_DATASET.keys():
    values = []
    for pdb_id, data in results["measurements"].items():
        if data["category"] == category:
            v = data["laws"].get("LAW-120", {}).get("observed")
            if v is not None:
                values.append(v)
    if values:
        import statistics
        print(f"  {category:<20} n={len(values):<3} mean={statistics.mean(values):6.2f}  min={min(values):6.2f}  max={max(values):6.2f}")
    else:
        print(f"  {category:<20} no data")
