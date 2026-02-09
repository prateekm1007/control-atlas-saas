#!/usr/bin/env python3
"""
Extract Split Motifs (Geometry-First, Resolution-Aware)
Accepts missing residues if geometry is resolved.
"""
import os
from prody import parsePDB, fetchPDB, writePDB

LIBRARY = os.path.expanduser("~/control-atlas/library/motifs/")
PDB_CACHE = os.path.expanduser("~/control-atlas/library/pdb_cache/")

EXTRACTIONS = [
    {"pdb": "6OIM", "chain": "A", "range": (30, 38), "out": "M001A_SwitchI_KRAS_G12C_drug_bound.pdb"},
    {"pdb": "6OIM", "chain": "A", "range": (60, 76), "out": "M001B_SwitchII_KRAS_G12C_drug_bound.pdb"},
    {"pdb": "5P21", "chain": "A", "range": (30, 38), "out": "M002A_SwitchI_HRAS_WT_native_gtp.pdb"},
    {"pdb": "5P21", "chain": "A", "range": (60, 76), "out": "M002B_SwitchII_HRAS_WT_native_gtp.pdb"},
]

def main():
    os.makedirs(LIBRARY, exist_ok=True)
    os.makedirs(PDB_CACHE, exist_ok=True)

    for ext in EXTRACTIONS:
        print(f"\nExtracting {ext['out']}")

        pdb_file = fetchPDB(ext["pdb"], folder=PDB_CACHE)
        structure = parsePDB(pdb_file)

        start, end = ext["range"]
        sel = structure.select(f"chain {ext['chain']} and resnum {start}:{end}")
        if sel is None:
            print("  ERROR: selection failed")
            continue

        ca = sel.select("calpha")
        if ca is None or ca.numAtoms() < 6:
            print("  ERROR: insufficient CA atoms")
            continue

        out_path = os.path.join(LIBRARY, ext["out"])
        writePDB(out_path, sel)
        print(f"  ✓ Saved motif with {ca.numAtoms()} CA atoms")

if __name__ == "__main__":
    main()
