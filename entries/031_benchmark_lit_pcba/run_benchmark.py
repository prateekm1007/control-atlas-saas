#!/usr/bin/env python3
"""
Entry 031 — LIT-PCBA Benchmark Runner
Evaluates EF1%, TNR, and FPR.
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

# Pre-computed reference physics (from Entry 027)
TARGET_PHYSICS = {
    "IDH1": {
        "status": "VALIDATED",
        "volume": 420.0,
        "hydrophobic_pct": 0.60,
        "exposure": 0.40
    },
    "TP53": {
        "status": "VALIDATED",
        "volume": 380.0,
        "hydrophobic_pct": 0.55,
        "exposure": 0.45
    }
}

def load_smiles(path):
    data = []
    if not path.exists(): return []
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                data.append((parts[0], parts[1]))
    return data

def run_target(target_name, data_dir, gk):
    actives = load_smiles(data_dir / "actives.smi")
    inactives = load_smiles(data_dir / "inactives.smi")
    
    if not actives or not inactives:
        return {"error": "Missing data"}

    print(f"    Actives: {len(actives)} | Inactives: {len(inactives)}")
    
    # Inject physics
    gk.catalog[target_name] = TARGET_PHYSICS[target_name]
    
    results = []
    
    # Screen Actives
    for smiles, cid in actives:
        r = gk.validate(target_name, smiles)
        score = r["metrics"].get("confidence", 0.0)
        results.append({"type": "active", "score": score, "status": r["status"]})
        
    # Screen Inactives
    for smiles, cid in inactives:
        r = gk.validate(target_name, smiles)
        score = r["metrics"].get("confidence", 0.0)
        results.append({"type": "inactive", "score": score, "status": r["status"]})
        
    results.sort(key=lambda x: x["score"], reverse=True)
    
    # Metrics
    n_total = len(results)
    n_actives = len(actives)
    top_1_pct_n = max(1, int(n_total * 0.01))
    
    # EF1%
    top_1_pct_set = results[:top_1_pct_n]
    actives_found = sum(1 for r in top_1_pct_set if r["type"] == "active")
    precision_1pct = actives_found / top_1_pct_n
    random_rate = n_actives / n_total
    ef_1pct = precision_1pct / random_rate if random_rate > 0 else 0
    
    # TNR (True Negative Rate) & FPR (False Positive Rate)
    rejected_inactives = sum(1 for r in results if r["type"] == "inactive" and r["status"] == "REJECT")
    tnr = rejected_inactives / len(inactives)
    fpr = 1.0 - tnr
    
    return {
        "target": target_name,
        "ef_1pct": round(ef_1pct, 2),
        "tnr": round(tnr, 3),
        "fpr": round(fpr, 3),
        "actives_recovered": actives_found,
        "inactives_rejected": rejected_inactives
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="results/benchmark_results.json")
    args = parser.parse_args()
    
    data_root = Path(__file__).parent / "data"
    gk = UniversalGatekeeper()
    summary = []
    
    print("="*60)
    print("Control Atlas — LIT-PCBA Benchmark (Synthetic Verification)")
    print("="*60)
    
    for target in ["IDH1", "TP53"]:
        print(f"\n[*] Benchmarking {target}...")
        res = run_target(target, data_root / target, gk)
        if "error" not in res:
            print(f"    EF1%: {res['ef_1pct']}x")
            print(f"    TNR (Rejection): {res['tnr']:.1%}")
            print(f"    FPR (False Pos): {res['fpr']:.1%}")
            summary.append(res)
            
    with open(args.out, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n[+] Results saved to {args.out}")

if __name__ == "__main__":
    main()
