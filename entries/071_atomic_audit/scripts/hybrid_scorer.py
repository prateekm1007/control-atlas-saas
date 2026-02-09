#!/usr/bin/env python3
print(">>> hybrid_scorer.py loaded")

import json
import random
from pathlib import Path

ROOT = Path.home() / "control-atlas"
CANDIDATES_FILE = ROOT / "entries/070_affinity_maturation/candidates_v5_mutants.json"
RESULTS_FILE = ROOT / "entries/071_atomic_audit/hybrid_results/final_ranking.csv"

def simulate_rf2_pae(warhead):
    pae = 10.0
    if any(x in warhead for x in "RK"):
        pae -= 3.0
    if any(x in warhead for x in "YW"):
        pae -= 2.0
    return max(2.0, pae + random.uniform(-1.0, 1.0))

def simulate_rosetta_dG(warhead, pae):
    if pae > 8.0:
        return random.uniform(-10.0, -5.0)
    dG = -20.0
    if any(x in warhead for x in "RK"):
        dG -= 8.0
    if any(x in warhead for x in "YW"):
        dG -= 5.0
    return dG + random.uniform(-2.0, 2.0)

def main():
    print(">>> main() entered")

    if not CANDIDATES_FILE.exists():
        raise FileNotFoundError(CANDIDATES_FILE)

    cands = json.load(open(CANDIDATES_FILE))
    print(f">>> Loaded {len(cands)} candidates")

    results = []

    for c in cands:
        pae = simulate_rf2_pae(c["warhead"])
        dG  = simulate_rosetta_dG(c["warhead"], pae)
        hybrid = pae + dG / 5.0
        results.append((c["id"], c["warhead"], pae, dG, hybrid))

    results.sort(key=lambda x: x[-1])

    print("\nTOP 5 CANDIDATES")
    for r in results[:5]:
        print(r)

    RESULTS_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_FILE, "w") as f:
        f.write("id,warhead,pae,dG,hybrid\n")
        for r in results:
            f.write(",".join(map(str, r)) + "\n")

    print("\n>>> Hybrid audit complete")
    print(f">>> Winner: {results[0][0]}")

if __name__ == "__main__":
    main()
