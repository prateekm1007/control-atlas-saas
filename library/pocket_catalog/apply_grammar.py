#!/usr/bin/env python3
"""
Apply Pan-Target Grammar (Fixed Scoring)
"""
import os
import json
import numpy as np
from prody import fetchPDB, parsePDB

CATALOG = os.path.expanduser("~/control-atlas/library/pocket_catalog")
GRAMMAR = os.path.join(CATALOG, "pan_target_grammar.json")
PDB_CACHE = os.path.expanduser("~/control-atlas/library/pdb_cache")

TEST_TARGETS = [
    {"target": "JAK2", "mutation": "V617F", "pdb": "4IVA", "lining": [610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621]},
    {"target": "FLT3", "mutation": "ITD", "pdb": "4XUF", "lining": [830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 842]},
    {"target": "IDH2", "mutation": "R140Q", "pdb": "4JA8", "lining": [132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 170, 171, 172]}
]

def compute_physics(pdb, lining):
    try:
        f = fetchPDB(pdb, folder=PDB_CACHE)
        if not f: return None
        s = parsePDB(f)
        sel = s.select(f"resnum {' '.join(map(str, lining))}")
        if not sel: return None
        
        coords = sel.getCoords()
        from scipy.spatial import ConvexHull
        vol = ConvexHull(coords).volume if len(coords) >= 4 else 0
        
        return {"volume_A3": vol, "exposure": 25.0, "hydrophobic_pct": 30.0} # Simplified proxy for debug
    except:
        return None

def classify(physics, grammar):
    best = "Unclassified"
    min_dist = float("inf")
    
    vol = physics["volume_A3"]
    
    for name, data in grammar["pocket_classes"].items():
        mean_vol = data["physics"]["volume"]["mean"]
        std_vol = data["physics"]["volume"]["std"]
        
        # Z-score distance (volume dominant)
        dist = abs(vol - mean_vol) / (std_vol + 100)
        
        if dist < min_dist:
            min_dist = dist
            best = name
            
    conf = max(0, 100 - min_dist * 10)
    return best, round(conf, 1)

def main():
    print("=== GRAMMAR TRANSFER DEBUG ===")
    with open(GRAMMAR) as f: grammar = json.load(f)
    
    for t in TEST_TARGETS:
        name = f"{t['target']}_{t['mutation']}"
        phys = compute_physics(t["pdb"], t["lining"])
        
        if not phys:
            print(f"{name}: Physics Failed")
            continue
            
        cls, conf = classify(phys, grammar)
        print(f"{name}: {cls} (Conf: {conf}%) | Vol: {phys['volume_A3']:.1f}")
        
        # Save
        d = os.path.join(CATALOG, name)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "grammar_prediction.json"), "w") as f:
            json.dump({"class": cls, "conf": conf}, f)

if __name__ == "__main__":
    main()
