#!/usr/bin/env python3
"""
Entry 025 — Batch Compound Screening (High-Throughput)

Robust dynamic loader for Entry 020.
Supports BOTH:
- class-based Gatekeeper (UniversalGatekeeper)
- function-based Gatekeeper (validate)

No physics. No duplication. Single initialization.
"""

import argparse
import csv
import sys
import importlib.util
from pathlib import Path
from time import time

# --- Dynamic Import of Entry 020 ---
GK_PATH = Path(__file__).resolve().parents[2] / "entries/020_candidate_validation/universal_gatekeeper.py"

spec = importlib.util.spec_from_file_location("universal_gatekeeper", GK_PATH)
gk_module = importlib.util.module_from_spec(spec)
sys.modules["universal_gatekeeper"] = gk_module
spec.loader.exec_module(gk_module)

# --- Resolve Gatekeeper API ---
if hasattr(gk_module, "UniversalGatekeeper"):
    MODE = "class"
    GatekeeperClass = gk_module.UniversalGatekeeper
elif hasattr(gk_module, "validate"):
    MODE = "function"
    validate_fn = gk_module.validate
else:
    sys.stderr.write(
        "CRITICAL: Entry 020 exposes neither UniversalGatekeeper nor validate().\n"
    )
    sys.exit(1)
# --------------------------------

def parse_smi(smi_path):
    with open(smi_path, "r") as f:
        for idx, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smiles = parts[0]
            cid = parts[1] if len(parts) > 1 else f"cmp_{idx+1:06d}"
            yield cid, smiles

def main():
    ap = argparse.ArgumentParser(description="Batch compound screening via Universal Gatekeeper")
    ap.add_argument("--target", required=True)
    ap.add_argument("--smi", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    print(f"[*] Initializing Gatekeeper ({MODE}) for {args.target}...")

    if MODE == "class":
        gk = GatekeeperClass()

    out_fields = [
        "compound_id", "smiles", "status", "confidence",
        "reason", "volume", "exposure", "hydrophobic_pct"
    ]

    count = 0
    t0 = time()

    with open(args.out, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=out_fields)
        writer.writeheader()

        for cid, smiles in parse_smi(args.smi):
            try:
                if MODE == "class":
                    result = gk.validate(args.target, smiles)
                else:
                    result = validate_fn(args.target, smiles)

                status = result.get("status", "ERROR")
                reasons = "; ".join(result.get("reasons", []))
                metrics = result.get("metrics", {}) or {}

                writer.writerow({
                    "compound_id": cid,
                    "smiles": smiles,
                    "status": status,
                    "confidence": metrics.get("confidence", 0.0),
                    "reason": reasons,
                    "volume": metrics.get("volume", ""),
                    "exposure": metrics.get("exposure", ""),
                    "hydrophobic_pct": metrics.get("hydrophobic_pct", "")
                })
                count += 1

            except Exception as e:
                writer.writerow({
                    "compound_id": cid,
                    "smiles": smiles,
                    "status": "ERROR",
                    "confidence": 0.0,
                    "reason": str(e),
                    "volume": "",
                    "exposure": "",
                    "hydrophobic_pct": ""
                })

    dt = time() - t0
    rate = count / dt if dt else 0
    print(f"[+] Complete: {count} compounds in {dt:.2f}s ({rate:.1f} cmp/s)")
    print(f"[+] Output: {args.out}")

if __name__ == "__main__":
    main()
