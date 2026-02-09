#!/usr/bin/env python3
"""
Entry 026 — Multi-Target Batch Screening

Deterministic orchestration over:
- Entry 020 (Universal Gatekeeper)
- Entry 025 (.smi parsing + batch loop pattern)

Screens ONE compound library against MANY targets.
"""

import argparse
import csv
import sys
import importlib.util
from pathlib import Path
from time import time

# --- Dynamic import of Entry 020 ---
GK_PATH = Path(__file__).resolve().parents[2] / "entries/020_candidate_validation/universal_gatekeeper.py"

spec = importlib.util.spec_from_file_location("universal_gatekeeper", GK_PATH)
gk_module = importlib.util.module_from_spec(spec)
sys.modules["universal_gatekeeper"] = gk_module
spec.loader.exec_module(gk_module)

if not hasattr(gk_module, "UniversalGatekeeper"):
    sys.stderr.write("CRITICAL: UniversalGatekeeper class not found in Entry 020\n")
    sys.exit(1)

Gatekeeper = gk_module.UniversalGatekeeper


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


def load_targets(arg_targets, validated_only):
    gk = Gatekeeper()
    catalog = gk.catalog

    if validated_only:
        return [
            t for t, v in catalog.items()
            if v.get("status") == "VALIDATED"
        ]

    return arg_targets


def main():
    ap = argparse.ArgumentParser(description="Multi-target batch screening")
    ap.add_argument("--smi", required=True, help="Compound library (.smi)")
    ap.add_argument("--out", required=True, help="Output CSV (long-form)")
    ap.add_argument("--targets", nargs="+", help="Explicit list of targets")
    ap.add_argument("--validated-only", action="store_true",
                    help="Screen all VALIDATED targets in catalog")

    args = ap.parse_args()

    if not args.targets and not args.validated_only:
        sys.stderr.write("ERROR: provide --targets or --validated-only\n")
        sys.exit(1)

    targets = load_targets(args.targets, args.validated_only)

    if not targets:
        sys.stderr.write("ERROR: no targets selected\n")
        sys.exit(1)

    compounds = list(parse_smi(args.smi))

    fields = [
        "compound_id", "smiles", "target",
        "status", "confidence",
        "volume", "exposure", "hydrophobic_pct"
    ]

    t0 = time()
    total = 0

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()

        for target in targets:
            print(f"[*] Initializing Gatekeeper for {target}")
            gk = Gatekeeper()

            for cid, smiles in compounds:
                try:
                    result = gk.validate(target, smiles, cid)
                    metrics = result.get("metrics", {}) or {}

                    writer.writerow({
                        "compound_id": cid,
                        "smiles": smiles,
                        "target": target,
                        "status": result.get("status", "ERROR"),
                        "confidence": metrics.get("confidence", 0.0),
                        "volume": metrics.get("volume", ""),
                        "exposure": metrics.get("exposure", ""),
                        "hydrophobic_pct": metrics.get("hydrophobic_pct", "")
                    })
                    total += 1
                except Exception as e:
                    writer.writerow({
                        "compound_id": cid,
                        "smiles": smiles,
                        "target": target,
                        "status": "ERROR",
                        "confidence": 0.0,
                        "volume": "",
                        "exposure": "",
                        "hydrophobic_pct": ""
                    })

    dt = time() - t0
    rate = total / dt if dt else 0
    print(f"[+] Complete: {total} evaluations in {dt:.2f}s ({rate:.1f} eval/s)")
    print(f"[+] Output: {args.out}")


if __name__ == "__main__":
    main()
