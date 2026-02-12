#!/usr/bin/env python3
"""
Download LIT-PCBA benchmark data.
Currently generates SYNTHETIC data for pipeline verification.

WARNING:
This synthetic data is for LOGIC VERIFICATION ONLY.
Do not use for scientific benchmarking results.
"""

import os
from pathlib import Path

TARGETS = ["IDH1", "TP53"]

def generate_synthetic_data(target: str, output_dir: Path):
    target_dir = output_dir / target
    target_dir.mkdir(parents=True, exist_ok=True)
    
    # Synthetic Actives (High confidence matches)
    actives = [
        ("CC(C)C1=CC(=C(C=C1)C2=CN(C3=CC=CC=C32)C4=CC(=CC=C4)C(=O)O)O", "IDH1_ACT_1"),
        ("CC1=CC(=C(C=C1)C2=CN(C3=CC=CC=C32)CC(=O)O)C(F)(F)F", "IDH1_ACT_2"),
        ("OC(=O)c1ccc(cc1)N2C=C(c3cccc(c3)Cl)c4ccccc24", "TP53_ACT_1"),
        ("COC(=O)c1ccc(cc1)N2C=C(c3ccc(c3)F)c4ccccc24", "TP53_ACT_2"),
        ("CC(C)(C)NC1=NC=NC2=C1N=CN2", "CONTROL_POS")
    ]
    
    # Synthetic Inactives (Physics mismatches)
    inactives = [
        ("C", "INACTIVE_001"), 
        ("CCCCCCCCCCCCCCCC", "INACTIVE_002"),
        ("CC1=CC=CC=C1", "INACTIVE_003"),
        ("O" * 10, "INACTIVE_004"),
        ("C" * 30, "INACTIVE_005"),
        ("N", "INACTIVE_006"), 
        ("Cl", "INACTIVE_007"),
        ("C(=O)(O)O", "INACTIVE_008"),
        ("c1ccccc1O", "INACTIVE_009"),
        ("CN(C)C", "INACTIVE_010")
    ]
    
    with open(target_dir / "actives.smi", "w") as f:
        for smi, cid in actives:
            f.write(f"{smi} {cid}\n")
            
    with open(target_dir / "inactives.smi", "w") as f:
        for smi, cid in inactives:
            f.write(f"{smi} {cid}\n")
            
    print(f"[+] Synthetic LIT-PCBA data generated for {target}")

def main():
    output_dir = Path(__file__).parent / "data"
    output_dir.mkdir(exist_ok=True)
    print("="*60)
    print("LIT-PCBA Data Manager (Synthetic Fallback)")
    print("="*60)
    for target in TARGETS:
        generate_synthetic_data(target, output_dir)
    print(f"\n[+] Data ready for {len(TARGETS)} targets")

if __name__ == "__main__":
    main()
