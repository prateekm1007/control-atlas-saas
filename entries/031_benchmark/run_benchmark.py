#!/usr/bin/env python3
"""
Entry 031 — DUD-E Benchmark Runner

Evaluates Control Atlas on DUD-E actives vs decoys.
Measures: Acceptance rate for actives, Rejection rate for decoys.
"""

import argparse
import csv
import json
import sys
from pathlib import Path
from time import time

BASE = Path(__file__).resolve().parents[2] / "entries"
sys.path.insert(0, str(BASE / "020_candidate_validation"))

from universal_gatekeeper import UniversalGatekeeper


def load_smiles(smi_path: Path, max_compounds: int = None) -> list:
    """Load SMILES from DUD-E format files."""
    compounds = []
    with open(smi_path, "r") as f:
        for i, line in enumerate(f):
            if max_compounds and i >= max_compounds:
                break
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                smiles, cid = parts[0], parts[1]
            else:
                smiles, cid = parts[0], f"cmp_{i}"
            compounds.append((cid, smiles))
    return compounds


def run_benchmark_target(
    target_name: str,
    target_dir: Path,
    gk: UniversalGatekeeper,
    max_actives: int = 100,
    max_decoys: int = 500
) -> dict:
    """Run benchmark on a single target."""
    
    # Load actives and decoys
    actives_file = target_dir / "actives_final.ism"
    decoys_file = target_dir / "decoys_final.ism"
    
    if not actives_file.exists() or not decoys_file.exists():
        return {"error": "Missing data files"}
    
    actives = load_smiles(actives_file, max_actives)
    decoys = load_smiles(decoys_file, max_decoys)
    
    print(f"    Actives: {len(actives)}, Decoys: {len(decoys)}")
    
    # Create a synthetic target in the catalog
    target_id = f"DUDE_{target_name.upper()}"
    gk.catalog[target_id] = {
        "status": "VALIDATED",
        "volume": 450.0,
        "hydrophobic_pct": 0.55,
        "exposure": 0.35
    }
    
    # Screen actives
    active_results = {"VALID": 0, "CANDIDATE": 0, "REJECT": 0, "ERROR": 0}
    for cid, smiles in actives:
        try:
            r = gk.validate(target_id, smiles)
            status = r.get("status", "ERROR")
            if status in active_results:
                active_results[status] += 1
            else:
                active_results["ERROR"] += 1
        except:
            active_results["ERROR"] += 1
    
    # Screen decoys
    decoy_results = {"VALID": 0, "CANDIDATE": 0, "REJECT": 0, "ERROR": 0}
    for cid, smiles in decoys:
        try:
            r = gk.validate(target_id, smiles)
            status = r.get("status", "ERROR")
            if status in decoy_results:
                decoy_results[status] += 1
            else:
                decoy_results["ERROR"] += 1
        except:
            decoy_results["ERROR"] += 1
    
    # Calculate metrics
    total_actives = len(actives)
    total_decoys = len(decoys)
    
    # True Positive Rate (actives accepted)
    actives_accepted = active_results["VALID"] + active_results["CANDIDATE"]
    tpr = actives_accepted / total_actives if total_actives > 0 else 0
    
    # True Negative Rate (decoys rejected)
    decoys_rejected = decoy_results["REJECT"]
    tnr = decoys_rejected / total_decoys if total_decoys > 0 else 0
    
    # Enrichment Factor (simplified)
    decoys_accepted = decoy_results["VALID"] + decoy_results["CANDIDATE"]
    total_accepted = actives_accepted + decoys_accepted
    if total_accepted > 0:
        ef = (actives_accepted / total_accepted) / (total_actives / (total_actives + total_decoys))
    else:
        ef = 0
    
    return {
        "target": target_name,
        "actives_total": total_actives,
        "decoys_total": total_decoys,
        "actives_accepted": actives_accepted,
        "actives_rejected": active_results["REJECT"],
        "decoys_accepted": decoys_accepted,
        "decoys_rejected": decoys_rejected,
        "tpr": round(tpr, 3),
        "tnr": round(tnr, 3),
        "enrichment_factor": round(ef, 2)
    }


def main():
    parser = argparse.ArgumentParser(description="DUD-E Benchmark")
    parser.add_argument("--data-dir", default="data", help="DUD-E data directory")
    parser.add_argument("--max-actives", type=int, default=100)
    parser.add_argument("--max-decoys", type=int, default=500)
    parser.add_argument("--out", default="results/benchmark_results.json")
    args = parser.parse_args()
    
    data_dir = Path(__file__).parent / args.data_dir
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    
    print("="*60)
    print("Control Atlas — DUD-E Benchmark")
    print("="*60)
    
    gk = UniversalGatekeeper()
    
    targets = ["egfr", "braf", "cdk2", "src", "vgfr2"]
    all_results = []
    
    t0 = time()
    
    for target in targets:
        target_dir = data_dir / target
        if not target_dir.exists():
            print(f"\n[!] {target}: Data not found (run download_dude.py first)")
            continue
        
        print(f"\n[*] Benchmarking {target.upper()}...")
        result = run_benchmark_target(
            target,
            target_dir,
            gk,
            max_actives=args.max_actives,
            max_decoys=args.max_decoys
        )
        
        if "error" not in result:
            all_results.append(result)
            print(f"    TPR (actives accepted): {result['tpr']:.1%}")
            print(f"    TNR (decoys rejected):  {result['tnr']:.1%}")
            print(f"    Enrichment Factor:      {result['enrichment_factor']:.1f}x")
    
    dt = time() - t0
    
    # Aggregate metrics
    if all_results:
        avg_tpr = sum(r["tpr"] for r in all_results) / len(all_results)
        avg_tnr = sum(r["tnr"] for r in all_results) / len(all_results)
        avg_ef = sum(r["enrichment_factor"] for r in all_results) / len(all_results)
        
        summary = {
            "benchmark": "DUD-E",
            "targets_evaluated": len(all_results),
            "total_time_s": round(dt, 2),
            "aggregate": {
                "avg_tpr": round(avg_tpr, 3),
                "avg_tnr": round(avg_tnr, 3),
                "avg_enrichment_factor": round(avg_ef, 2)
            },
            "per_target": all_results
        }
        
        print("\n" + "="*60)
        print("AGGREGATE RESULTS")
        print("="*60)
        print(f"  Targets:             {len(all_results)}")
        print(f"  Avg TPR (actives):   {avg_tpr:.1%}")
        print(f"  Avg TNR (decoys):    {avg_tnr:.1%}")
        print(f"  Avg Enrichment:      {avg_ef:.1f}x")
        print(f"  Time:                {dt:.1f}s")
        
        # Save results
        out_path = Path(__file__).parent / args.out
        with open(out_path, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"\n[+] Results saved to {out_path}")


if __name__ == "__main__":
    main()
