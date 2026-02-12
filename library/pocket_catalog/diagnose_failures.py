#!/usr/bin/env python3
"""
Diagnose PDB Structure Issues
Inspects chains and residue ranges for failed targets.
"""

import os
from prody import fetchPDB, parsePDB

PDB_CACHE = os.path.expanduser("~/control-atlas/library/pdb_cache")
FAILED_TARGETS = [
    {"name": "JAK2_V617F", "pdb": "4IVA", "expected_res": 617},
    {"name": "EGFR_T790M", "pdb": "5Y3A", "expected_res": 790},
    {"name": "MET_Y1003",  "pdb": "3DKC", "expected_res": 1003}
]

def analyze_pdb(target):
    print(f"--- Analyzing {target['name']} ({target['pdb']}) ---")
    try:
        pdb = fetchPDB(target['pdb'], folder=PDB_CACHE)
        if not pdb:
            print("  [FAIL] Download failed")
            return

        structure = parsePDB(pdb)
        hier = structure.getHierView()
        
        for chain in hier:
            chid = chain.getChid()
            resnums = chain.getResnums()
            if len(resnums) == 0: continue
            
            min_res = min(resnums)
            max_res = max(resnums)
            count = len(resnums)
            
            # Check if expected residue is in this chain
            has_target = target['expected_res'] in resnums
            marker = "  <-- CONTAINS TARGET" if has_target else ""
            
            print(f"  Chain {chid}: Residues {min_res} to {max_res} ({count} atoms){marker}")
            
            # If target exists, print local context to verify sequence
            if has_target:
                try:
                    r = chain.getResidue(target['expected_res'])
                    print(f"    Target Residue: {r.getResname()} {r.getResnum()}")
                except:
                    pass

    except Exception as e:
        print(f"  [ERROR] {e}")
    print("")

if __name__ == "__main__":
    for t in FAILED_TARGETS:
        analyze_pdb(t)
