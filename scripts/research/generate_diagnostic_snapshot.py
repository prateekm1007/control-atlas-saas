import requests
import json
import os
import sys
from datetime import datetime, timezone

TARGETS = [
    ("4HHB", "PASS"),
    ("1CRN", "PASS"),
    ("1UBQ", "PASS"),
    ("6LZG", "PASS"),
    ("7BV2", "PASS"),
    ("1G03", "PASS"),
    ("AF-P01308-F1", "INDETERMINATE")
]

API_URL = "http://localhost:8000/ingest"
snapshot = {
    "meta": {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "phase": "1.5-diagnostic-lock",
        "note": "Pre-calibration physics baseline. Contains expected VETOs due to uncalibrated thresholds (LAW-120/160)."
    },
    "benchmarks": {}
}

print(f"{'Target':<15} | {'Expected':<15} | {'Actual':<15} | {'Breaking Laws'}")
print("-" * 80)

for sid, expected in TARGETS:
    try:
        resp = requests.post(API_URL, data={"mode": "Benchmark", "candidate_id": sid, "t3_category": "NONE"})
        d = resp.json()
        actual = d["verdict"]["binary"]
        
        # Identify breakers (Deterministic laws that are not PASS)
        breakers = []
        for l in d["tier1"]["laws"]:
            if l["method"] == "deterministic" and l["status"] != "PASS":
                breakers.append(f"{l['law_id']}({l['observed']})")
        
        breaker_str = ", ".join(breakers) if breakers else "—"
        
        print(f"{sid:<15} | {expected:<15} | {actual:<15} | {breaker_str}")
        
        # Capture full physics state for this target
        snapshot["benchmarks"][sid] = {
            "verdict": actual,
            "expected": expected,
            "laws": {l["law_id"]: l for l in d["tier1"]["laws"]}
        }
        
    except Exception as e:
        print(f"{sid:<15} | ERROR: {str(e)}")

# Write invariant file
os.makedirs("backend/invariants", exist_ok=True)
filepath = "backend/invariants/v23.0.phase1.snapshot.json"
with open(filepath, "w") as f:
    json.dump(snapshot, f, indent=4)

print("-" * 80)
print(f"✅ DIAGNOSTIC SNAPSHOT LOCKED: {filepath}")
print("This file contains the evidence required for Phase 2 Calibration.")
