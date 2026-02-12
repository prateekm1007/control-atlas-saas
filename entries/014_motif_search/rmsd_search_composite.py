#!/usr/bin/env python3
"""
Entry #014d: Composite Motif Geometry Search
Probes:
- Switch I (9 residues)
- Switch II (17 residues)
Constraint:
- Spatial distance between probe centroids
"""

import os
from prody import parsePDB, fetchPDB, superpose, calcRMSD
import numpy as np

# Paths
MOTIF_PATH = os.path.expanduser("~/control-atlas/library/motifs/M001_Switch_GTPase.pdb")
PDB_CACHE = os.path.expanduser("~/control-atlas/library/pdb_cache/")

# Targets
CANDIDATES = [
    ("5P21", "A"),  # H-Ras (control)
    ("2RGN", "A"),  # RhoA (discovery)
    ("1UBQ", "A"),  # Ubiquitin (negative)
]

RMSD_A_THRESHOLD = 2.5  # Switch II
RMSD_B_THRESHOLD = 2.5  # Switch I
DIST_TOLERANCE = 2.0   # Å

def sliding_window_best(target_coords, probe_coords):
    best_rmsd = float("inf")
    best_centroid = None

    for i in range(len(target_coords) - len(probe_coords) + 1):
        window = target_coords[i:i+len(probe_coords)].copy()
        aligned, _ = superpose(window, probe_coords)
        rmsd = calcRMSD(aligned, probe_coords)

        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_centroid = aligned.mean(axis=0)

    return best_rmsd, best_centroid

def main():
    print("Loading reference motif (M001)")
    motif = parsePDB(MOTIF_PATH)
    ca = motif.select("calpha")

    # Split probes
    probe_B = ca[:9].getCoords()      # Switch I
    probe_A = ca[-17:].getCoords()    # Switch II

    # Reference centroid distance
    ref_dist = np.linalg.norm(probe_A.mean(axis=0) - probe_B.mean(axis=0))
    print(f"Reference A–B distance: {ref_dist:.2f} Å\n")

    print(f"{'TARGET':<8} | A_RMSD | B_RMSD | ΔDist | STATUS")
    print("-" * 50)

    for pdb_id, chain in CANDIDATES:
        try:
            pdb = fetchPDB(pdb_id, folder=PDB_CACHE)
            struct = parsePDB(pdb)
            target = struct.select(f"chain {chain} and calpha")
            coords = target.getCoords()

            rmsd_A, cen_A = sliding_window_best(coords, probe_A)
            rmsd_B, cen_B = sliding_window_best(coords, probe_B)

            if cen_A is None or cen_B is None:
                status = "REJECT"
                delta = float("inf")
            else:
                dist = np.linalg.norm(cen_A - cen_B)
                delta = abs(dist - ref_dist)

                status = (
                    "MATCH"
                    if rmsd_A <= RMSD_A_THRESHOLD
                    and rmsd_B <= RMSD_B_THRESHOLD
                    and delta <= DIST_TOLERANCE
                    else "REJECT"
                )

            print(f"{pdb_id:<8} | {rmsd_A:6.2f} | {rmsd_B:6.2f} | {delta:5.2f} | {status}")

        except Exception as e:
            print(f"{pdb_id:<8} | ERROR: {e}")

if __name__ == "__main__":
    main()
