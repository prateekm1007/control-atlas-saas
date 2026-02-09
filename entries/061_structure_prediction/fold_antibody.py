#!/usr/bin/env python3
"""
Entry 061 — Antibody Folding
Uses ESMFold API (reused from Entry 028) to fold antibody candidates.
"""

import sys
import json
import argparse
from pathlib import Path

# Reuse Entry 028
BASE = Path(__file__).resolve().parents[2] / "entries"
sys.path.insert(0, str(BASE / "028_structure_prediction"))
from structure_provider import StructureProvider

def fold_candidates(candidates_file, output_dir):
    with open(candidates_file) as f:
        candidates = json.load(f)
        
    provider = StructureProvider()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results = []
    print(f"[*] Folding {len(candidates)} antibodies...")
    
    for cand in candidates:
        cid = cand["id"]
        # Create scFv: Heavy + Linker + Light
        linker = "GGGGSGGGGSGGGGS"
        full_seq = cand["heavy_chain"] + linker + cand["light_chain"]
        
        print(f" -> Folding {cid}...")
        res = provider.predict(full_seq)
        
        status = "FOLDED"
        if res["status"] != "SUCCESS":
            status = "FAILED"
        elif res["confidence_global"] < 0.70:
            status = "UNSTABLE"
            
        result_entry = {
            "id": cid,
            "status": status,
            "plddt": res.get("confidence_global", 0.0),
            "pdb_path": res.get("pdb_path"),
            "sequence_length": len(full_seq)
        }
        results.append(result_entry)
        
        # Save Metadata
        with open(output_dir / f"{cid}_meta.json", "w") as f:
            json.dump(result_entry, f, indent=2)
            
    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--out", default="folding_results.json")
    args = parser.parse_args()
    
    results = fold_candidates(args.candidates, "folded_antibodies")
    
    with open(args.out, "w") as f:
        json.dump(results, f, indent=2)
        
    folded = sum(1 for r in results if r["status"] == "FOLDED")
    print(f"\n[+] Complete. {folded}/{len(results)} antibodies folded successfully.")

if __name__ == "__main__":
    main()
