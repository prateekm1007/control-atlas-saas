#!/usr/bin/env python3
"""
Entry 015b — Pocket-Frame Manifold (Physics Metrics)
"""

import os
import csv
import numpy as np
from prody import parsePDB, fetchPDB
from scipy.spatial import ConvexHull

PDB_CACHE = os.path.expanduser("~/control-atlas/library/pdb_cache/")
OUTPUT_CSV = os.path.expanduser("~/control-atlas/entries/015_state_manifold/pocket_physics.csv")

STATES = [
    {"pdb": "4OBE", "chain": "A", "state": "GDP_closed"},
    {"pdb": "6GOD", "chain": "A", "state": "GTP_intermediate"},
    {"pdb": "6OIM", "chain": "A", "state": "Drug_open"},
    {"pdb": "5P21", "chain": "A", "state": "Native_GTP"},
]

SWITCH_II = (60, 76)
HYDROPHOBIC = {'ALA','VAL','LEU','ILE','MET','PHE','TRP','PRO'}

def exposure_proxy(structure, chain):
    pocket = structure.select(f"chain {chain} and resnum {SWITCH_II[0]}:{SWITCH_II[1]}")
    protein = structure.select("protein")
    if pocket is None: return None, None

    pcoords = pocket.getCoords()
    allc = protein.getCoords()

    exposures = []
    for p in pcoords:
        d = np.linalg.norm(allc - p, axis=1)
        exposures.append(100 - min(np.sum((d > 0) & (d < 8.0)), 100))

    return round(np.sum(exposures),1), round(np.mean(exposures),2)

def pocket_volume(structure, chain):
    sel = structure.select(f"chain {chain} and resnum {SWITCH_II[0]}:{SWITCH_II[1]}")
    if sel is None or sel.numAtoms() < 4: return None
    return round(ConvexHull(sel.getCoords()).volume,1)

def hydrophobic_fraction(structure, chain):
    sel = structure.select(f"chain {chain} and resnum {SWITCH_II[0]}:{SWITCH_II[1]} and name CA")
    if sel is None: return None
    return round(100 * sum(r in HYDROPHOBIC for r in sel.getResnames()) / len(sel),1)

def main():
    os.makedirs(PDB_CACHE, exist_ok=True)
    results = []

    for s in STATES:
        print(f"\n{s['pdb']} ({s['state']})")
        pdb = fetchPDB(s["pdb"], folder=PDB_CACHE)
        st = parsePDB(pdb)

        e_tot, e_mean = exposure_proxy(st, s["chain"])
        vol = pocket_volume(st, s["chain"])
        hyd = hydrophobic_fraction(st, s["chain"])

        print(f"  Exposure mean: {e_mean}")
        print(f"  Volume: {vol}")
        print(f"  Hydrophobic %: {hyd}")

        results.append({
            "pdb_id": s["pdb"],
            "state": s["state"],
            "exposure_mean": e_mean,
            "volume": vol,
            "hydrophobic_pct": hyd
        })

    with open(OUTPUT_CSV,"w",newline="") as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys())
        w.writeheader()
        w.writerows(results)

    print("\nPOCKET PHYSICS SUMMARY")
    for r in results:
        print(f"{r['state']:<18} exp={r['exposure_mean']} vol={r['volume']} hyd={r['hydrophobic_pct']}")

if __name__ == "__main__":
    main()
