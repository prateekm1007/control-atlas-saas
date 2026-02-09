#!/usr/bin/env python3
"""
Entry 049 — Synthetic Library Generator
"""

def generate_library(output_path="synthetic_1k.smi"):
    scaffolds = [
        "c1ccccc1", "C1CCCCC1", "c1ncncc1",
        "C1CCNCC1", "c1ccccc1-c2ccccc2", "C(=O)N"
    ]
    groups = ["C", "O", "N", "F", "Cl", "Br", "C(=O)O", "CN(C)C"]

    compounds = []
    count = 0
    for s in scaffolds:
        for g1 in groups:
            for g2 in groups:
                smi = f"{g1}{s}{g2}"
                compounds.append((f"SYN_{count:04d}", smi))
                count += 1
                if count >= 1000:
                    break
            if count >= 1000:
                break

    # Known controls
    compounds.append(("SOTORASIB", "CC(C)(C)NC1=NC=NC2=C1N=CN2"))
    compounds.append(("DECOY_TINY", "C"))
    compounds.append(("DECOY_HUGE", "C" * 30))

    with open(output_path, "w") as f:
        for cid, smi in compounds:
            f.write(f"{smi} {cid}\n")

    print(f"[+] Generated {len(compounds)} compounds → {output_path}")

if __name__ == "__main__":
    generate_library()
