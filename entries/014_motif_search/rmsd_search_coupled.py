#!/usr/bin/env python3
"""
Entry #014e: Coupled Composite Motif Search
Switch I + Switch II with enforced coupling and spatial constraint.
"""

import os
import numpy as np
from prody import parsePDB, fetchPDB, superpose, calcRMSD

# Paths
MOTIF_PATH = os.path.expanduser("~/control-atlas/library/motifs/M001_Switch_GTPase.pdb")
PDB_CACHE = os.path.expanduser("~/control-atlas/library/pdb_cache/")

# Targets
CANDIDATES = [
    ("5P21", "A"),  # H-Ras (positive control)
    ("2RGN", "A"),  # RhoA
    ("1UBQ", "A"),  # Ubiquitin (negative control)
]

# Thresholds
RMSD_I_MAX = 2.5
RMSD_II_MAX = 2.5
DIST_TOL = 2.0        # Å
WINDOW_LINK = 40      # Max residue separation between Switch I and II

def centroid(coords):
    return coords.mean(axis=0)

def main():
    print("Loading reference motif (M001)")
    motif = parsePDB(MOTIF_PATH)
    ca = motif.select("calpha")

    probe_I = ca[:9].getCoords()      # Switch I
    probe_II = ca[-17:].getCoords()   # Switch II

    ref_dist = np.linalg.norm(centroid(probe_I) - centroid(probe_II))
    print(f"Reference I–II distance: {ref_dist:.2f} Å\n")

    print(f"{'TARGET':<8} | I_RMSD | II_RMSD | ΔDist | STATUS")
    print("-" * 50)

    for pdb_id, chain in CANDIDATES:
        try:
            pdb = fetchPDB(pdb_id, folder=PDB_CACHE)
            struct = parsePDB(pdb)
            ca_t = struct.select(f"chain {chain} and calpha")
            coords = ca_t.getCoords()

            best = None

            # Slide Switch II
            for j in range(len(coords) - len(probe_II)):
                win_II = coords[j:j+len(probe_II)].copy()
                aln_II, _ = superpose(win_II, probe_II)
                rmsd_II = calcRMSD(aln_II, probe_II)
                if rmsd_II > RMSD_II_MAX:
                    continue

                # Search Switch I nearby
                for i in range(
                    max(0, j-WINDOW_LINK),
                    min(len(coords)-len(probe_I), j+WINDOW_LINK)
                ):
                    win_I = coords[i:i+len(probe_I)].copy()
                    aln_I, _ = superpose(win_I, probe_I)
                    rmsd_I = calcRMSD(aln_I, probe_I)
                    if rmsd_I > RMSD_I_MAX:
                        continue

                    dist = np.linalg.norm(centroid(aln_I) - centroid(aln_II))
                    delta = abs(dist - ref_dist)
                    if delta > DIST_TOL:
                        continue

                    best = (rmsd_I, rmsd_II, delta)
                    break

                if best:
                    break

            if best:
                rmsd_I, rmsd_II, delta = best
                status = "MATCH"
            else:
                rmsd_I = rmsd_II = delta = float("inf")
                status = "REJECT"

            print(f"{pdb_id:<8} | {rmsd_I:6.2f} | {rmsd_II:7.2f} | {delta:5.2f} | {status}")

        except Exception as e:
            print(f"{pdb_id:<8} | ERROR: {e}")

if __name__ == "__main__":
    main()
