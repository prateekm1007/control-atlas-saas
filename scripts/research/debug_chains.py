from prody import parsePDB, fetchPDB
import os

# WP-10: repo-relative cache (was the hardcoded ~/control-atlas/... path).
PDB_CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "library", "pdb_cache", "")
os.makedirs(PDB_CACHE, exist_ok=True)
pdb_file = fetchPDB("6OIM", folder=PDB_CACHE)
structure = parsePDB(pdb_file)

print(f"Total atoms: {structure.numAtoms()}")

# Get all chains
all_ca = structure.select("name CA")
if all_ca:
    chains = set(all_ca.getChids())
    print(f"Available chains: {chains}")
    
    for chain in chains:
        sel = structure.select(f"chain {chain} and name CA")
        resnums = sorted(set(sel.getResnums()))
        print(f"Chain {chain}: residues {min(resnums)}-{max(resnums)} ({len(resnums)} total)")
else:
    print("No CA atoms found")
