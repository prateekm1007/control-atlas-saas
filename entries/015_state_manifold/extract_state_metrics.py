#!/usr/bin/env python3
"""
Entry 015 — Control Pocket State Manifold Metrics
"""

import os
import csv
import numpy as np
from prody import parsePDB, fetchPDB, calcCenter

PDB_CACHE = os.path.expanduser("~/control-atlas/library/pdb_cache/")
OUTPUT_CSV = os.path.expanduser("~/control-atlas/entries/015_state_manifold/state_manifold.csv")

STATES = [
    {"pdb": "4OBE", "chain": "A", "state": "GDP_closed"},
    {"pdb": "6GOD", "chain": "A", "state": "GTP_intermediate"},
    {"pdb": "6OIM", "chain": "A", "state": "Drug_open"},
    {"pdb": "5P21", "chain": "A", "state": "Native_GTP"},
]

SWITCH_II = (60, 76)

def extract_metrics(structure, chain):
    sel = f"chain {chain} and resnum {SWITCH_II[0]}:{SWITCH_II[1]} and calpha"
    atoms = structure.select(sel)
    if atoms is None:
        return None

    coords = atoms.getCoords()
    centroid = calcCenter(atoms)
    distances = np.linalg.norm(coords - centroid, axis=1)

    return {
        "n_atoms": atoms.numAtoms(),
        "spread": round(np.std(distances), 3),
        "max_extent": round(np.max(distances), 3),
        "z_range": round(np.ptp(coords[:, 2]), 3),
    }

def main():
    os.makedirs(PDB_CACHE, exist_ok=True)
    results = []

    for s in STATES:
        print(f"Processing {s['pdb']} ({s['state']})...")
        pdb = fetchPDB(s["pdb"], folder=PDB_CACHE)
        struct = parsePDB(pdb)
        metrics = extract_metrics(struct, s["chain"])
        if metrics:
            metrics["pdb_id"] = s["pdb"]
            metrics["state"] = s["state"]
            results.append(metrics)
            print(f"  Spread={metrics['spread']}  Extent={metrics['max_extent']}")

    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["pdb_id", "state", "n_atoms", "spread", "max_extent", "z_range"]
        )
        writer.writeheader()
        writer.writerows(results)

    print("\n=== SUMMARY ===")
    for r in results:
        print(
            f"{r['state']:<18} "
            f"spread={r['spread']:>6} "
            f"extent={r['max_extent']:>6} "
            f"z={r['z_range']:>6}"
        )

if __name__ == "__main__":
    main()
