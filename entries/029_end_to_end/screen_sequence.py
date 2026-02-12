#!/usr/bin/env python3
"""
Entry 029 — End-to-End Sequence Screening Pipeline
With Decision Provenance Export (Entry 030)
"""

import argparse
import csv
import json
import sys
from pathlib import Path
from time import time

BASE = Path(__file__).resolve().parents[2] / "entries"

sys.path.insert(0, str(BASE / "028_structure_prediction"))
sys.path.insert(0, str(BASE / "027_pocket_detection"))
sys.path.insert(0, str(BASE / "020_candidate_validation"))
sys.path.insert(0, str(BASE / "030_provenance"))

from structure_provider import StructureProvider
from pocket_detector import PocketDetector
from universal_gatekeeper import UniversalGatekeeper
from provenance import generate_provenance, generate_batch_provenance, save_provenance


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
    parser = argparse.ArgumentParser(description="End-to-End Sequence Screening")
    parser.add_argument("--sequence", required=True)
    parser.add_argument("--smi", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--max-pockets", type=int, default=3)
    parser.add_argument("--provenance-dir", help="Directory for provenance records")
    parser.add_argument("--json-trace", help="Optional trace output")
    args = parser.parse_args()

    t0 = time()

    # STAGE 1: Structure
    print("\n" + "="*60)
    print("STAGE 1: Structure Prediction")
    print("="*60)

    sp = StructureProvider()
    sr = sp.predict(args.sequence)

    if sr["status"] != "SUCCESS":
        print(f"[!] Failed: {sr['error']}")
        sys.exit(1)

    pdb_path = sr["pdb_path"]
    struct_conf = sr["confidence_global"]
    print(f"[+] Structure: {pdb_path}")
    print(f"[+] Confidence: {struct_conf:.1%}")

    # STAGE 2: Pockets
    print("\n" + "="*60)
    print("STAGE 2: Pocket Detection")
    print("="*60)

    pd = PocketDetector(max_pockets=args.max_pockets)
    pr = pd.detect(pdb_path, struct_conf)

    if pr["status"] != "SUCCESS":
        print(f"[!] Failed: {pr['error']}")
        sys.exit(1)

    pockets = [p for p in pr["pockets"] if p["status"] in ("VALIDATED", "CANDIDATE")]

    if not pockets:
        print("[!] No druggable pockets found.")
        sys.exit(0)

    print(f"[+] Screenable pockets: {len(pockets)}")
    for p in pockets:
        print(f"    - {p['pocket_id']}: {p['status']} (conf: {p['confidence']:.2f})")

    # STAGE 3: Screening
    print("\n" + "="*60)
    print("STAGE 3: Compound Screening")
    print("="*60)

    gk = UniversalGatekeeper()
    for p in pockets:
        tid = f"SEQ_{p['pocket_id']}"
        gk.catalog[tid] = {
            "status": "VALIDATED",
            "volume": p["volume"],
            "hydrophobic_pct": p["hydrophobic_pct"],
            "exposure": p["exposure"]
        }

    compounds = list(parse_smi(args.smi))
    print(f"[*] Screening {len(compounds)} compounds x {len(pockets)} pockets")

    fields = ["compound_id", "smiles", "pocket", "status", "confidence", "reason"]
    results = []
    provenance_records = []

    for cid, smi in compounds:
        for p in pockets:
            tid = f"SEQ_{p['pocket_id']}"
            try:
                r = gk.validate(tid, smi)
                m = r.get("metrics", {}) or {}
                conf = struct_conf * p["confidence"] * m.get("confidence", 0)
                
                decision = {
                    "status": r.get("status", "ERROR"),
                    "confidence": round(conf, 3),
                    "reasons": r.get("reasons", [])
                }
                
                results.append({
                    "compound_id": cid,
                    "smiles": smi,
                    "pocket": p["pocket_id"],
                    "status": decision["status"],
                    "confidence": decision["confidence"],
                    "reason": "; ".join(decision["reasons"])
                })
                
                # Generate provenance
                if args.provenance_dir:
                    prov = generate_provenance(
                        sequence=args.sequence,
                        compound_id=cid,
                        smiles=smi,
                        decision=decision,
                        structure_result=sr,
                        pocket_data=p,
                        checks_passed=[],
                        checks_failed=decision["reasons"]
                    )
                    provenance_records.append(prov)
                    save_provenance(prov, Path(args.provenance_dir), cid)
                    
            except Exception as e:
                results.append({
                    "compound_id": cid, "smiles": smi, "pocket": p["pocket_id"],
                    "status": "ERROR", "confidence": 0, "reason": str(e)
                })

    results.sort(key=lambda x: x["confidence"], reverse=True)

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    dt = time() - t0
    v = sum(1 for r in results if r["status"] == "VALID")
    c = sum(1 for r in results if r["status"] == "CANDIDATE")
    rj = sum(1 for r in results if r["status"] == "REJECT")

    print(f"\n[+] Done: {len(results)} screens in {dt:.2f}s")
    print(f"[+] Output: {args.out}")
    print(f"    VALID: {v} | CANDIDATE: {c} | REJECT: {rj}")

    # Save batch provenance
    if args.provenance_dir:
        batch_prov = generate_batch_provenance(args.sequence, sr, results, pockets)
        batch_path = Path(args.provenance_dir) / "batch_provenance.json"
        with open(batch_path, "w") as f:
            json.dump(batch_prov, f, indent=2)
        print(f"[+] Provenance: {args.provenance_dir}/")

    if args.json_trace:
        trace = {
            "total_time": round(dt, 2),
            "valid": v, "candidate": c, "reject": rj
        }
        with open(args.json_trace, "w") as f:
            json.dump(trace, f, indent=2)


if __name__ == "__main__":
    main()
