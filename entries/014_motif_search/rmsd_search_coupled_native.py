#!/usr/bin/env python3
"""
Entry #014e (Native Geometry)
Coupled Composite Motif Search using native HRAS motifs
"""

import os
import numpy as np
from prody import parsePDB, fetchPDB, superpose, calcRMSD

LIB = os.path.expanduser("~/control-atlas/library/motifs/")
PDB_CACHE = os.path.expanduser("~/control-atlas/library/pdb_cache/")

TARGETS = [
    ("5P21", "A"),  # HRAS native (positive control)
    ("2RGN", "A"),  # RhoA
    ("1UBQ", "A"),  # Ubiquitin (negative)
]

RMSD_MAX = 2.5
DIST_TOL = 2.5
WINDOW = 40

def centroid(x):
    return x.mean(axis=0)

def main():
    ref_I = parsePDB(os.path.join(
        LIB, "M002A_SwitchI_HRAS_WT_native_gtp.pdb"
    )).select("calpha").getCoords()

    ref_II = parsePDB(os.path.join(
        LIB, "M002B_SwitchII_HRAS_WT_native_gtp.pdb"
    )).select("calpha").getCoords()

    ref_dist = np.linalg.norm(centroid(ref_I) - centroid(ref_II))

    print(f"Reference I–II distance: {ref_dist:.2f} Å\n")
    print(f"{'TARGET':<8} | I_RMSD | II_RMSD | ΔDist | STATUS")
    print("-" * 50)

    for pdb_id, chain in TARGETS:
        pdb = fetchPDB(pdb_id, folder=PDB_CACHE)
        struct = parsePDB(pdb)
        ca = struct.select(f"chain {chain} and calpha").getCoords()

        found = False

        for j in range(len(ca) - len(ref_II)):
            wII = ca[j:j+len(ref_II)].copy()
            aII, _ = superpose(wII, ref_II)
            rII = calcRMSD(aII, ref_II)
            if rII > RMSD_MAX:
                continue

            for i in range(
                max(0, j - WINDOW),
                min(len(ca) - len(ref_I), j + WINDOW)
            ):
                wI = ca[i:i+len(ref_I)].copy()
                aI, _ = superpose(wI, ref_I)
                rI = calcRMSD(aI, ref_I)
                if rI > RMSD_MAX:
                    continue

                d = abs(
                    np.linalg.norm(centroid(aI) - centroid(aII)) - ref_dist
                )

                if d <= DIST_TOL:
                    print(f"{pdb_id:<8} | {rI:6.2f} | {rII:7.2f} | {d:5.2f} | MATCH")
                    found = True
                    break

            if found:
                break

        if not found:
            print(f"{pdb_id:<8} |   inf  |   inf   |  inf  | REJECT")

if __name__ == "__main__":
    main()
