#!/usr/bin/env python3
"""
Entry #014a: Homology-Guided Motif Check via RMSD
Finds structural matches to M001_Switch_GTPase across GTPase superfamily
"""

import os
import csv
from prody import parsePDB, matchChains, calcRMSD, fetchPDB
import numpy as np

# Paths
MOTIF_PATH = os.path.expanduser("~/control-atlas/library/motifs/M001_Switch_GTPase.pdb")
OUTPUT_CSV = os.path.expanduser("~/control-atlas/library/matches/M001_hits.csv")
PDB_CACHE = os.path.expanduser("~/control-atlas/library/pdb_cache/")

# Candidates (PDB ID, expected chain)
CANDIDATES = [
    ("5P21", "A"),  # H-Ras GTP
    ("1QRA", "A"),  # H-Ras
    ("1YZN", "A"),  # Ran
    ("1HUR", "A"),  # Arf1
    ("2RGN", "A"),  # RhoA
    ("4Q21", "A"),  # p21-Ras
]

RMSD_THRESHOLD = 2.5  # Angstroms

def main():
    os.makedirs(PDB_CACHE, exist_ok=True)
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    
    # Load reference motif
    print(f"Loading reference motif: {MOTIF_PATH}")
    motif = parsePDB(MOTIF_PATH)
    motif_ca = motif.select("calpha")
    
    if motif_ca is None:
        print("ERROR: No CA atoms in motif")
        return
    
    print(f"Motif CA atoms: {motif_ca.numAtoms()}")
    
    results = []
    
    for pdb_id, chain in CANDIDATES:
        print(f"\nProcessing {pdb_id} chain {chain}...")
        
        try:
            # Fetch structure
            pdb_file = fetchPDB(pdb_id, folder=PDB_CACHE)
            if pdb_file is None:
                print(f"  SKIP: Could not fetch {pdb_id}")
                continue
                
            structure = parsePDB(pdb_file)
            target_chain = structure.select(f"chain {chain} and calpha")
            
            if target_chain is None:
                print(f"  SKIP: No CA atoms in chain {chain}")
                continue
            
            # Simple RMSD (first N atoms matching)
            n_atoms = min(motif_ca.numAtoms(), target_chain.numAtoms())
            
            if n_atoms < 10:
                print(f"  SKIP: Too few atoms ({n_atoms})")
                continue
            
            # Calculate RMSD on matching segment
            motif_coords = motif_ca.getCoords()[:n_atoms]
            target_coords = target_chain.getCoords()[:n_atoms]
            
            rmsd = np.sqrt(np.mean(np.sum((motif_coords - target_coords)**2, axis=1)))
            
            status = "MATCH" if rmsd <= RMSD_THRESHOLD else "REJECT"
            print(f"  RMSD: {rmsd:.2f} Å → {status}")
            
            results.append({
                "pdb_id": pdb_id,
                "chain": chain,
                "rmsd": round(rmsd, 2),
                "n_atoms": n_atoms,
                "status": status
            })
            
        except Exception as e:
            print(f"  ERROR: {e}")
            continue
    
    # Write results
    if results:
        with open(OUTPUT_CSV, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["pdb_id", "chain", "rmsd", "n_atoms", "status"])
            writer.writeheader()
            writer.writerows(results)
        print(f"\nResults saved to: {OUTPUT_CSV}")
    
    # Summary
    matches = [r for r in results if r["status"] == "MATCH"]
    print(f"\n=== SUMMARY ===")
    print(f"Candidates: {len(CANDIDATES)}")
    print(f"Processed: {len(results)}")
    print(f"Matches (RMSD ≤ {RMSD_THRESHOLD}Å): {len(matches)}")

if __name__ == "__main__":
    main()
