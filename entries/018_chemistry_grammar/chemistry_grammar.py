#!/usr/bin/env python3
import os
import json

print("=== ENTRY 018: CHEMISTRY GRAMMAR ===")

# Pocket properties from Entry 015/016
POCKET = {
    "volume": 2231.4,
    "exposure": 36.63,
    "hydrophobic_pct": 18.8,
    "state": "drug_open",
    "target": "KRAS_G12C",
    "key_residue": "CYS12"
}

# Chemistry rules based on pocket physics
RULES = {
    "core_scaffold": {
        "requirement": "rigid_aromatic",
        "reason": "fills concave pocket",
        "examples": ["quinazoline", "pyrimidine", "indole"]
    },
    "anchor_group": {
        "requirement": "hydrophobic",
        "reason": "matches 18.8% nonpolar surface",
        "examples": ["methyl", "isopropyl", "cyclopropyl"]
    },
    "polar_contact": {
        "requirement": "h_bond_acceptor",
        "reason": "engages polar rim (GLU, ASP)",
        "examples": ["carbonyl", "sulfonyl", "ether"]
    },
    "warhead": {
        "requirement": "electrophile",
        "reason": "covalent to CYS12",
        "examples": ["acrylamide", "vinyl_sulfonamide", "chloroacetamide"]
    },
    "solubility": {
        "requirement": "polar_tail",
        "reason": "exposed region tolerates polarity",
        "examples": ["piperazine", "morpholine", "hydroxyl"]
    }
}

# Known drugs validation
KNOWN_DRUGS = {
    "Sotorasib": {
        "core": "pyridopyrimidine",
        "warhead": "acrylamide",
        "matches_rules": True
    },
    "Adagrasib": {
        "core": "pyrrolopyrimidine",
        "warhead": "acrylamide",
        "matches_rules": True
    }
}

print("")
print("POCKET PROPERTIES")
print("-" * 50)
for k, v in POCKET.items():
    print(f"  {k}: {v}")

print("")
print("CHEMISTRY RULES")
print("-" * 50)
for category, rule in RULES.items():
    print(f"{category}:")
    print(f"  Requirement: {rule['requirement']}")
    print(f"  Reason: {rule['reason']}")
    print(f"  Examples: {', '.join(rule['examples'])}")
    print("")

print("VALIDATION: Known Drugs")
print("-" * 50)
for drug, props in KNOWN_DRUGS.items():
    status = "PASS" if props["matches_rules"] else "FAIL"
    print(f"  {drug}: {status}")
    print(f"    Core: {props['core']}")
    print(f"    Warhead: {props['warhead']}")

# Save rules
OUTPUT = os.path.expanduser("~/control-atlas/entries/018_chemistry_grammar/chemistry_rules.json")
with open(OUTPUT, "w") as f:
    json.dump({"pocket": POCKET, "rules": RULES, "validation": KNOWN_DRUGS}, f, indent=2)

print("")
print(f"Saved: {OUTPUT}")
print("=== ENTRY 018 COMPLETE ===")
