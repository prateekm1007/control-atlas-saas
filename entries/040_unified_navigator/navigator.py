#!/usr/bin/env python3
"""
Entry 040 — Unified Navigation Engine
Integrates Physics, Chemistry, Biology, and Math layers.
"""

import argparse
import sys
import json
from pathlib import Path

# Import all layers
BASE = Path(__file__).resolve().parents[2] / "entries"
sys.path.insert(0, str(BASE / "037_chemistry_layer"))
sys.path.insert(0, str(BASE / "038_biology_layer"))
sys.path.insert(0, str(BASE / "039_math_layer"))

# Mock Physics (Atlas lookup) - In prod, this queries the JSON files
ATLAS_DIR = Path(__file__).resolve().parents[2] / "library/atlas_index"

from chemistry_engine import ChemistryLayer
from biology_engine import BiologyLayer
from math_engine import MathLayer

class UnifiedNavigator:
    def __init__(self):
        self.chem = ChemistryLayer()
        self.bio = BiologyLayer()
        self.math = MathLayer()
        
    def get_physics_data(self, target_name):
        # Mock lookup for demo targets (since we don't have UniProt mapping here yet)
        # In prod: Map Target -> UniProt -> Atlas JSON
        if target_name == "KRAS":
            return {"status": "INDEXED", "volume": 463.1, "hydrophobic_pct": 0.31, "pockets": 3}
        elif target_name == "EGFR":
            return {"status": "NO_POCKETS"} 
        return {"status": "UNKNOWN"}

    def navigate(self, target, smiles):
        print(f"[*] Navigating {target} with molecule...")
        
        # 1. Physics Layer (Hard Veto)
        phy_data = self.get_physics_data(target)
        if phy_data["status"] != "INDEXED":
            return {"status": "BLOCKED", "reason": f"Physics Veto: {phy_data['status']}"}
        
        # 2. Biology Layer (Relevance Filter)
        bio_res = self.bio.evaluate({"gene_name": target})
        if bio_res.status == "FAIL":
            return {"status": "BLOCKED", "reason": f"Biology Veto: {bio_res.reasons}"}
            
        # 3. Chemistry Layer (Tractability)
        chem_res = self.chem.evaluate({"smiles": smiles, "pocket_metrics": phy_data})
        if chem_res.status == "FAIL":
            return {"status": "BLOCKED", "reason": f"Chemistry Veto: {chem_res.reasons}"}
            
        # 4. Math Layer (The Clearing)
        robustness = self.math.calculate_robustness(phy_data, chem_res.metrics, bio_res.metrics)
        
        return {
            "status": "CLEARED",
            "confidence": robustness.robustness,
            "advice": robustness.advice,
            "trace": {
                "physics": phy_data,
                "chemistry": chem_res.metrics,
                "biology": bio_res.metrics
            }
        }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", required=True)
    parser.add_argument("--smiles", required=True)
    args = parser.parse_args()
    
    nav = UnifiedNavigator()
    result = nav.navigate(args.target, args.smiles)
    
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
