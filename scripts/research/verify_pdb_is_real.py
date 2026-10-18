from pathlib import Path
import sys

pdb = Path("outputs/v16_P4/final_model.pdb")

if not pdb.exists():
    print("⛔ PDB missing")
    sys.exit(1)

with open(pdb) as f:
    lines = f.readlines()

atoms = [l for l in lines if l.startswith("ATOM")]

if not atoms:
    print("⛔ PDB has no ATOM records — placeholder detected")
    sys.exit(1)

print("✅ PDB contains atomic coordinates")
