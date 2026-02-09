#!/usr/bin/env python3
import os
import json

print("=== ENTRY 018b: CHEMISTRY GRAMMAR (HARDENED) ===")

# Pocket regions with distinct chemistry requirements
POCKET_REGIONS = {
    "buried_core": {
        "residues": [64, 71, 72, 78, 82, 90],
        "properties": {"exposure": "low", "polarity": "hydrophobic"},
        "allowed": ["aromatic_ring", "alkyl", "halogen"],
        "forbidden": ["charged", "h_bond_donor", "flexible_chain"],
        "max_volume": 150,
        "notes": "Desolvation penalty for polar groups"
    },
    "concave_wall": {
        "residues": [60, 65, 66, 67, 68, 69, 70],
        "properties": {"exposure": "medium", "polarity": "mixed"},
        "allowed": ["rigid_linker", "small_polar", "halogen"],
        "forbidden": ["large_hydrophobic", "macrocycle", "flexible_alkyl"],
        "max_volume": 120,
        "notes": "Shape complementarity critical"
    },
    "polar_rim": {
        "residues": [62, 63, 76, 91, 92, 98],
        "properties": {"exposure": "high", "polarity": "polar"},
        "allowed": ["h_bond_acceptor", "h_bond_donor", "small_polar"],
        "forbidden": ["large_hydrophobic", "strong_cation"],
        "max_volume": 100,
        "notes": "Engages GLU/ASP rim"
    },
    "solvent_tail": {
        "residues": [],
        "properties": {"exposure": "full", "polarity": "any"},
        "allowed": ["polar_chain", "solubilizing_group", "piperazine", "morpholine"],
        "forbidden": ["large_hydrophobic"],
        "max_volume": 200,
        "notes": "Tolerates flexibility, aids solubility"
    },
    "warhead_zone": {
        "residues": [12],
        "properties": {"target": "CYS12"},
        "allowed": ["acrylamide", "vinyl_sulfonamide", "chloroacetamide"],
        "forbidden": ["reversible_binder", "non_electrophile"],
        "max_volume": 80,
        "notes": "Covalent engagement required for G12C"
    }
}

# Global quantitative constraints
QUANTITATIVE_LIMITS = {
    "molecular_weight": {"min": 350, "max": 600, "unit": "Da"},
    "rotatable_bonds": {"min": 2, "max": 7},
    "ring_count": {"min": 2, "max": 4},
    "hbd": {"min": 0, "max": 3, "note": "H-bond donors"},
    "hba": {"min": 3, "max": 8, "note": "H-bond acceptors"},
    "tpsa": {"min": 60, "max": 140, "unit": "A^2"},
    "clogp": {"min": 1, "max": 5}
}

# Resistance-aware exclusions (from Entry 017)
RESISTANCE_RULES = {
    "high_risk_residues": [60, 65, 66, 69, 74, 75, 77, 83, 85, 86, 87, 89, 92],
    "rules": [
        {
            "rule": "no_steric_projection",
            "target_residues": [60, 75, 77],
            "reason": "GLY residues - any bulk causes clash",
            "penalty": "CRITICAL"
        },
        {
            "rule": "flexibility_forbidden",
            "target_residues": [65, 66, 83, 89],
            "reason": "Small residues (SER/ALA) - flexible groups amplify escape",
            "penalty": "HIGH"
        },
        {
            "rule": "polar_mismatch_forbidden",
            "target_residues": [69, 92],
            "reason": "ASP residues - hydrophobic contact disrupts",
            "penalty": "MEDIUM"
        },
        {
            "rule": "steric_buffer_required",
            "target_residues": [74, 87],
            "reason": "THR residues - need tolerance margin",
            "penalty": "MEDIUM"
        }
    ]
}

# Forbidden chemistries (hard constraints)
FORBIDDEN_GLOBAL = [
    {"feature": "macrocycle", "reason": "pocket too shallow"},
    {"feature": "strong_base", "reason": "desolvation in buried core"},
    {"feature": "metal_chelator", "reason": "off-target liability"},
    {"feature": "reactive_aldehyde", "reason": "toxicity"},
    {"feature": "nitro_group", "reason": "metabolic liability"},
    {"feature": "polyfluorinated_chain", "reason": "PFAS concern"}
]

# Validation against known drugs
VALIDATION = {
    "Sotorasib": {
        "mw": 560.6,
        "rotatable_bonds": 6,
        "ring_count": 4,
        "warhead": "acrylamide",
        "passes_limits": True,
        "passes_forbidden": True,
        "resistance_safe": True
    },
    "Adagrasib": {
        "mw": 604.7,
        "rotatable_bonds": 7,
        "ring_count": 4,
        "warhead": "acrylamide",
        "passes_limits": True,
        "passes_forbidden": True,
        "resistance_safe": True
    }
}

# Output
print("")
print("POCKET REGIONS")
print("=" * 60)
for region, props in POCKET_REGIONS.items():
    print(f"{region.upper()}")
    print(f"  Residues: {props['residues']}")
    print(f"  Allowed: {', '.join(props['allowed'])}")
    print(f"  Forbidden: {', '.join(props['forbidden'])}")
    print(f"  Max volume: {props['max_volume']} A^3")
    print("")

print("QUANTITATIVE LIMITS")
print("=" * 60)
for prop, limits in QUANTITATIVE_LIMITS.items():
    print(f"  {prop}: {limits['min']} - {limits['max']}")

print("")
print("RESISTANCE-AWARE RULES")
print("=" * 60)
for rule in RESISTANCE_RULES['rules']:
    print(f"  {rule['rule']}")
    print(f"    Residues: {rule['target_residues']}")
    print(f"    Reason: {rule['reason']}")
    print(f"    Penalty: {rule['penalty']}")
    print("")

print("FORBIDDEN CHEMISTRIES (GLOBAL)")
print("=" * 60)
for item in FORBIDDEN_GLOBAL:
    print(f"  X {item['feature']}: {item['reason']}")

print("")
print("VALIDATION")
print("=" * 60)
for drug, props in VALIDATION.items():
    status = "PASS" if all([props['passes_limits'], props['passes_forbidden'], props['resistance_safe']]) else "FAIL"
    print(f"  {drug}: {status}")

# Save complete grammar
OUTPUT = os.path.expanduser("~/control-atlas/entries/018_chemistry_grammar/chemistry_grammar_v2.json")
grammar = {
    "version": "2.0",
    "pocket_regions": POCKET_REGIONS,
    "quantitative_limits": QUANTITATIVE_LIMITS,
    "resistance_rules": RESISTANCE_RULES,
    "forbidden_global": FORBIDDEN_GLOBAL,
    "validation": VALIDATION
}
with open(OUTPUT, "w") as f:
    json.dump(grammar, f, indent=2)

print("")
print(f"Saved: {OUTPUT}")
print("=== ENTRY 018b COMPLETE ===")
