#!/usr/bin/env python3
"""
Entry 026 — Multi-Target Batch Screening

Screens compound libraries against ALL VALIDATED targets.
Pure orchestration over UniversalGatekeeper (Entry 020).
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

if not hasattr(gk_module, "UniversalGatekeeper"):
    sys.stderr.write("CRITICAL: UniversalGatekeeper not found in Entry 020\n")
    sys.exit(1)

UniversalGatekeeper = gk_module.UniversalGatekeeper
# -----------------------------------

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

def get_validated_targets(gk, explicit=None):
    validated = {
        t for t, d in gk.catalog.items()
        if d.get("status") == "VALIDATED"
    }
    if explicit:
        return sorted(t for t in explicit if t in validated)
    return sorted(validated)

def main():
    ap = argparse.ArgumentParser(description="Multi-target batch screening")
    ap.add_argument("--smi", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--targets", default="ALL")
    ap.add_argument("--summary")
    args = ap.parse_args()

    print("[*] Initializing Universal Gatekeeper...")
    gk = UniversalGatekeeper()

    if args.targets.upper() == "ALL":
        targets = get_validated_targets(gk)
    else:
        explicit = [t.strip() for t in args.targets.split(",")]
        targets = get_validated_targets(gk, explicit)

    if not targets:
        sys.stderr.write("ERROR: No VALIDATED targets selected\n")
        sys.exit(1)

    compounds = list(parse_smi(args.smi))
    total = len(compounds) * len(targets)

    print(f"[*] {len(compounds)} compounds × {len(targets)} targets = {total} screens")

    fields = [
        "compound_id", "smiles", "target",
        "status", "confidence",
        "reason", "volume", "exposure", "hydrophobic_pct"
    ]

    best = {t: {"compound_id": None, "confidence": -1, "status": None} for t in targets}

    t0 = time()
    count = 0

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()

        for cid, smiles in compounds:
            for target in targets:
                try:
                    result = gk.validate(target, smiles, cid)
                    metrics = result.get("metrics", {}) or {}
                    conf = metrics.get("confidence", 0.0)

                    writer.writerow({
                        "compound_id": cid,
                        "smiles": smiles,
                        "target": target,
                        "status": result.get("status", "ERROR"),
                        "confidence": conf,
                        "reason": "; ".join(result.get("reasons", [])),
                        "volume": metrics.get("volume", ""),
                        "exposure": metrics.get("exposure", ""),
                        "hydrophobic_pct": metrics.get("hydrophobic_pct", "")
                    })

                    if conf > best[target]["confidence"]:
                        best[target] = {
                            "compound_id": cid,
                            "confidence": conf,
                            "status": result.get("status")
                        }

                    count += 1
                except Exception as e:
                    writer.writerow({
                        "compound_id": cid,
                        "smiles": smiles,
                        "target": target,
                        "status": "ERROR",
                        "confidence": 0.0,
                        "reason": str(e),
                        "volume": "",
                        "exposure": "",
                        "hydrophobic_pct": ""
                    })

    dt = time() - t0
    print(f"[+] Complete: {count} screens in {dt:.2f}s ({count/dt:.1f} screens/s)")
    print(f"[+] Output: {args.out}")

    if args.summary:
        with open(args.summary, "w", newline="") as sf:
            sw = csv.DictWriter(sf, fieldnames=["target", "top_compound", "confidence", "status"])
            sw.writeheader()
            for t in targets:
                sw.writerow({
                    "target": t,
                    "top_compound": best[t]["compound_id"],
                    "confidence": best[t]["confidence"],
                    "status": best[t]["status"]
                })
        print(f"[+] Summary: {args.summary}")

if __name__ == "__main__":
    main()
