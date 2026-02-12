#!/usr/bin/env python3
"""
Entry #014c: Geometry-Only Motif Search (Switch II Focus) [FINAL]
Target: Continuous Switch II (Residues 60–76)
Method: Sliding Window + Kabsch Superposition + calcRMSD
"""

import os
from prody import parsePDB, superpose, calcRMSD, fetchPDB

# Paths
MOTIF_PATH = os.path.expanduser("~/control-atlas/library/motifs/M001_Switch_GTPase.pdb")
PDB_CACHE = os.path.expanduser("~/control-atlas/library/pdb_cache/")

# Candidates
CANDIDATES = [
    ("5P21", "A"),  # H-Ras (control)
    ("2RGN", "A"),  # RhoA (discovery target)
    ("1YZN", "A"),  # Ran
    ("1HUR", "A"),  # Arf1
    ("1UBQ", "A"),  # Ubiquitin (negative control)
]

RMSD_THRESHOLD = 2.5

def sliding_window_search(target_ca, motif_coords):
    target_coords = target_ca.getCoords()
    motif_len = motif_coords.shape[0]

    best_rmsd = float("inf")
    best_start = -1

    for i in range(len(target_coords) - motif_len + 1):
        window = target_coords[i:i+motif_len].copy()  # CRITICAL COPY
        aligned_window, _ = superpose(window, motif_coords)
        rmsd = calcRMSD(aligned_window, motif_coords)

        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_start = i

    return best_rmsd, best_start

def main():
    print(f"Loading motif: {MOTIF_PATH}")
    motif = parsePDB(MOTIF_PATH)
    motif_ca = motif.select("calpha")[-17:]  # Switch II only
    motif_coords = motif_ca.getCoords()

    print(f"Probe size: {len(motif_coords)} residues (Switch II)\n")

    print(f"{'TARGET':<8} | {'RMSD (Å)':<10} | STATUS | LOC")
    print("-" * 45)

    for pdb_id, chain in CANDIDATES:
        try:
            pdb_file = fetchPDB(pdb_id, folder=PDB_CACHE)
            structure = parsePDB(pdb_file)
            target_ca = structure.select(f"chain {chain} and calpha")

            rmsd, start = sliding_window_search(target_ca, motif_coords)
            status = "MATCH" if rmsd <= RMSD_THRESHOLD else "REJECT"
            resnum = target_ca.getResnums()[start] if start != -1 else "N/A"

            print(f"{pdb_id:<8} | {rmsd:<10.2f} | {status:<6} | Res {resnum}")

        except Exception as e:
            print(f"{pdb_id:<8} | ERROR: {e}")

if __name__ == "__main__":
    main()
