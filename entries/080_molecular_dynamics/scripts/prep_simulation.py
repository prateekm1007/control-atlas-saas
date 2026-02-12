from pdbfixer import PDBFixer
from openmm.app import *
from openmm import *
from openmm.unit import *

INPUT_PDB = "../inputs/lead_structure.pdb"
OUTPUT_PDB = "../inputs/system_solvated.pdb"

print("🔧 [Entry 080] Preparing system for MD...")

fixer = PDBFixer(filename=INPUT_PDB)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.4)

forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
modeller = Modeller(fixer.topology, fixer.positions)

print("💧 Adding solvent...")
modeller.addSolvent(forcefield, boxSize=Vec3(7.0,7.0,7.0)*nanometer, ionicStrength=0.15*molar)
with open(OUTPUT_PDB, 'w') as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

print("✅ Solvated system written to:", OUTPUT_PDB)
