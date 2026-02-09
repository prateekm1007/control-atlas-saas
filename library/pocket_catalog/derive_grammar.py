#!/usr/bin/env python3
"""
Pan-Target Grammar Derivation (v1.1)
Robust to missing/failed targets
"""

import os
import json
import numpy as np

CATALOG = os.path.expanduser("~/control-atlas/library/pocket_catalog")
OUTPUT = os.path.expanduser("~/control-atlas/library/pocket_catalog/pan_target_grammar.json")

# Pocket class definitions (Updated based on successful physics)
POCKET_CLASSES = {
    "SwitchII": ["KRAS_G12C", "KRAS_G12D", "NRAS_Q61R"], # HRAS failed
    "KinaseHinge": ["BRAF_V600E", "EGFR_L858R", "ALK_F1174L", "RET_M918T"],
    "ActivationLoop": ["KIT_D816V", "PDGFRA_D842V", "FGFR3_S249C"],
    "Metabolic": ["IDH1_R132H"],
    "PHDomain": ["AKT1_E17K"],
    "Kinase": ["PIK3CA_H1047R"],
    "Juxtamembrane": ["MET_Y1003"],
    "FRBDomain": ["MTOR_S2215Y"]
}

def load_physics(target):
    path = os.path.join(CATALOG, target, "physics_metrics.json")
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return None

def derive_class_rules(targets):
    volumes = []
    exposures = []
    hydros = []
    
    for t in targets:
        m = load_physics(t)
        if m and m.get("status") == "computed":
            volumes.append(m["volume_A3"])
            exposures.append(m["exposure"])
            hydros.append(m["hydrophobic_pct"])
    
    if not volumes: return None
    
    return {
        "volume": {"mean": round(np.mean(volumes), 1), "std": round(np.std(volumes), 1)},
        "exposure": {"mean": round(np.mean(exposures), 2), "std": round(np.std(exposures), 2)},
        "hydrophobicity": {"mean": round(np.mean(hydros), 1), "std": round(np.std(hydros), 1)},
        "n_targets": len(volumes)
    }

def derive_chemistry_rules(stats):
    rules = {}
    vol = stats["volume"]["mean"]
    exp = stats["exposure"]["mean"]
    
    # Core Size
    if vol > 1500: rules["core_size"] = "large (3-4 rings)"
    elif vol > 1000: rules["core_size"] = "medium (2-3 rings)"
    else: rules["core_size"] = "small (1-2 rings)"
    
    # Polar Tolerance
    if exp > 30: rules["polar_tolerance"] = "high"
    elif exp > 20: rules["polar_tolerance"] = "medium"
    else: rules["polar_tolerance"] = "low"
    
    return rules

def main():
    print("=== PAN-TARGET GRAMMAR DERIVATION (v1.1) ===")
    grammar = {"pocket_classes": {}, "universal_rules": {}, "class_specific_rules": {}}
    
    all_volumes = []
    
    for class_name, targets in POCKET_CLASSES.items():
        stats = derive_class_rules(targets)
        if stats:
            grammar["pocket_classes"][class_name] = {"targets": targets, "physics": stats}
            grammar["class_specific_rules"][class_name] = derive_chemistry_rules(stats)
            print(f"{class_name}: Vol={stats['volume']['mean']}")
            
            # Collect universal stats
            for t in targets:
                m = load_physics(t)
                if m and m.get("status") == "computed":
                    all_volumes.append(m["volume_A3"])

    grammar["universal_rules"] = {
        "volume_range": [min(all_volumes), max(all_volumes)],
        "common_requirements": {"aromatic_core": True}
    }
    
    with open(OUTPUT, "w") as f:
        json.dump(grammar, f, indent=2)
    print(f"Grammar saved to: {OUTPUT}")

if __name__ == "__main__":
    main()
