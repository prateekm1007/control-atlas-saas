#!/usr/bin/env python3
import sys
from openmm.app import *
from openmm import *
from openmm.unit import *
from pdbfixer import PDBFixer

INPUT_PDB = "../inputs/lead_structure.pdb" # The one we fixed chain IDs for

def count_atoms(topology, stage):
    atoms = list(topology.atoms())
    print(f"[{stage}] Atom Count: {len(atoms)} | Chains: {topology.getNumChains()}")
    # Check chain lengths
    for i, chain in enumerate(topology.chains()):
        print(f"  -> Chain {chain.id}: {len(list(chain.residues()))} residues")

def main():
    print("🔍 [Entry 080] Debugging MD Preparation Pipeline...")
    
    # 1. LOAD
    print("1. Loading PDBFixer...")
    fixer = PDBFixer(filename=INPUT_PDB)
    count_atoms(fixer.topology, "Initial Load")

    # 2. FIND MISSING
    fixer.findMissingResidues()
    print(f"   Missing Residues: {len(fixer.missingResidues)}")
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    count_atoms(fixer.topology, "Post-Standardize")

    fixer.findMissingAtoms()
    print(f"   Missing Atoms: {len(fixer.missingAtoms)}")
    
    # 3. ADD ATOMS
    fixer.addMissingAtoms()
    count_atoms(fixer.topology, "Post-AddAtoms")

    # 4. ADD HYDROGENS
    fixer.addMissingHydrogens(7.4)
    count_atoms(fixer.topology, "Post-Hydrogens")

    # 5. MODELLER TRANSFER
    print("5. Transferring to Modeller...")
    modeller = Modeller(fixer.topology, fixer.positions)
    count_atoms(modeller.topology, "Modeller Init")

    # 6. FORCEFIELD LOAD
    print("6. Loading Forcefield (Amber14)...")
    try:
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        # This step often deletes atoms if they don't match templates
        # We can't easily check 'system' atom count without creating it, 
        # but Modeller.addSolvent calls createSystem internally or uses templates.
        print("   Forcefield loaded.")
    except Exception as e:
        print(f"❌ Forcefield Error: {e}")
        return

    # 7. SOLVATION
    print("7. Adding Solvent...")
    try:
        modeller.addSolvent(forcefield, padding=1.0*nanometers, ionicStrength=0.15*molar)
        count_atoms(modeller.topology, "Post-Solvation")
    except Exception as e:
        print(f"❌ Solvation Error: {e}")
        # If this fails, it's likely because the topology doesn't match the forcefield
        # Print residue names to check for weirdness
        print("   Dumping Residue Names for inspection...")
        for res in modeller.topology.residues():
            if res.name not in forcefield._templates:
                # This check is internal/complex, simplified check:
                pass
        return

if __name__ == "__main__":
    main()
