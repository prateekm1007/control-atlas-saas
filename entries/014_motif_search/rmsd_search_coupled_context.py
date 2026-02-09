#!/usr/bin/env python3
"""
Entry #014f: Coupled Composite Geometry + Context Constraint
Adds orientation (principal-axis alignment) to reject false positives.
"""

import os
import numpy as np
from prody import parsePDB, fetchPDB, superpose, calcRMSD

LIB = os.path.expanduser("~/control-atlas/library/motifs/")
PDB_CACHE = os.path.expanduser("~/control-atlas/library/pdb_cache/")

TARGETS = [
    ("5P21", "A"),
    ("2RGN", "A"),
    ("1UBQ", "A"),
]

RMSD_MAX = 2.5
DIST_TOL = 2.5
AXIS_COS_MIN = 0.85
WINDOW = 40

def centroid(x):
    return x.mean(axis=0)

def principal_axis(coords):
    # PCA first principal component
    X = coords - coords.mean(axis=0)
    _, _, vh = np.linalg.svd(X, full_matrices=False)
    axis = vh[0]
    return axis / np.linalg.norm(axis)

def axis_cos(a, b):
    return abs(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))

def main():
    ref_I = parsePDB(os.path.join(LIB,"M002A_SwitchI_HRAS_WT_native_gtp.pdb")).select("calpha").getCoords()
    ref_II = parsePDB(os.path.join(LIB,"M002B_SwitchII_HRAS_WT_native_gtp.pdb")).select("calpha").getCoords()

    ref_dist = np.linalg.norm(centroid(ref_I) - centroid(ref_II))
    ref_axis_I = principal_axis(ref_I)
    ref_axis_II = principal_axis(ref_II)

    print(f"Ref I–II dist: {ref_dist:.2f} Å | Axis cos min: {AXIS_COS_MIN}\n")
    print(f"{'TARGET':<8} | I_RMSD | II_RMSD | ΔDist | cosI | cosII | STATUS")
    print("-"*78)

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

            axis_II = principal_axis(aII)
            cosII = axis_cos(axis_II, ref_axis_II)
            if cosII < AXIS_COS_MIN:
                continue

            for i in range(max(0, j-WINDOW), min(len(ca)-len(ref_I), j+WINDOW)):
                wI = ca[i:i+len(ref_I)].copy()
                aI, _ = superpose(wI, ref_I)
                rI = calcRMSD(aI, ref_I)
                if rI > RMSD_MAX:
                    continue

                axis_I = principal_axis(aI)
                cosI = axis_cos(axis_I, ref_axis_I)
                if cosI < AXIS_COS_MIN:
                    continue

                d = abs(np.linalg.norm(centroid(aI) - centroid(aII)) - ref_dist)
                if d <= DIST_TOL:
                    print(f"{pdb_id:<8} | {rI:6.2f} | {rII:7.2f} | {d:5.2f} | {cosI:4.2f} | {cosII:5.2f} | MATCH")
                    found = True
                    break
            if found:
                break

        if not found:
            print(f"{pdb_id:<8} |   inf  |   inf   |  inf  |  --  |  --   | REJECT")

if __name__ == "__main__":
    main()
