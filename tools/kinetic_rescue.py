"""
Kinetic Rescue: Energy Minimization for Clash Resolution
Tests whether static clashes resolve under physics simulation.
"""

import sys
import json
from pathlib import Path

def run_rescue(cif_path, output_path=None):
    """
    Run energy minimization and re-audit structure.
    
    Returns dict with before/after metrics.
    """
    try:
        import openmm.app as app
        import openmm as mm
        from openmm import unit
        from pdbfixer import PDBFixer
    except ImportError:
        print("ERROR: Install openmm and pdbfixer first:")
        print("  pip install openmm pdbfixer")
        sys.exit(1)
    
    from Bio.PDB import MMCIFParser, PDBIO
    from scipy.spatial.distance import cdist
    import numpy as np
    import tempfile
    
    print(f"="*60)
    print(f"KINETIC RESCUE: {cif_path}")
    print(f"="*60)
    
    # --- PRE-MINIMIZATION AUDIT ---
    print("\n[1/4] Pre-minimization audit...")
    
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("complex", cif_path)[0]
    
    def get_heavy_atoms(chain):
        return [a for a in chain.get_atoms() 
                if a.element != 'H' 
                and a.get_parent().id[0] == ' ']
    
    target_atoms = get_heavy_atoms(structure['A'])
    binder_atoms = get_heavy_atoms(structure['B'])
    
    t_coords = np.array([a.coord for a in target_atoms])
    b_coords = np.array([a.coord for a in binder_atoms])
    
    pre_distances = cdist(t_coords, b_coords)
    pre_min_dist = float(np.min(pre_distances))
    pre_clashes_25 = int(np.sum(pre_distances < 2.5))
    pre_clashes_19 = int(np.sum(pre_distances < 1.9))
    pre_rho = int(np.sum(pre_distances < 4.5))
    
    print(f"   Min distance:      {pre_min_dist:.2f} Å")
    print(f"   Clashes (< 2.5Å):  {pre_clashes_25}")
    print(f"   Clashes (< 1.9Å):  {pre_clashes_19}")
    print(f"   Contact density:   {pre_rho}")
    
    # --- CONVERT TO PDB FOR OPENMM ---
    print("\n[2/4] Preparing for minimization...")
    
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        tmp_pdb = tmp.name
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(tmp_pdb)
    
    # --- ENERGY MINIMIZATION ---
    print("\n[3/4] Running energy minimization (100 steps)...")
    
    try:
        fixer = PDBFixer(tmp_pdb)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)
        
        forcefield = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
        system = forcefield.createSystem(
            fixer.topology, 
            nonbondedMethod=app.CutoffNonPeriodic,
            constraints=app.HBonds
        )
        
        integrator = mm.LangevinIntegrator(
            310*unit.kelvin, 
            1/unit.picosecond, 
            0.002*unit.picoseconds
        )
        
        simulation = app.Simulation(fixer.topology, system, integrator)
        simulation.context.setPositions(fixer.positions)
        
        # Get initial energy
        state = simulation.context.getState(getEnergy=True)
        initial_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"   Initial energy: {initial_energy:.1f} kJ/mol")
        
        # Minimize
        simulation.minimizeEnergy(maxIterations=100)
        
        # Get final energy
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"   Final energy:   {final_energy:.1f} kJ/mol")
        print(f"   Energy change:  {final_energy - initial_energy:.1f} kJ/mol")
        
        # Save minimized structure
        if output_path:
            with open(output_path, 'w') as f:
                app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
            print(f"   Saved to: {output_path}")
        
        minimization_success = True
        
    except Exception as e:
        print(f"   ERROR during minimization: {e}")
        minimization_success = False
    
    # --- POST-MINIMIZATION AUDIT ---
    print("\n[4/4] Post-minimization audit...")
    
    if minimization_success and output_path and Path(output_path).exists():
        from Bio.PDB import PDBParser
        parser2 = PDBParser(QUIET=True)
        min_structure = parser2.get_structure("minimized", output_path)[0]
        
        # Find chains (may have different IDs after OpenMM)
        chains = list(min_structure.get_chains())
        if len(chains) >= 2:
            chain_a = chains[0]
            chain_b = chains[1]
            
            post_target = get_heavy_atoms(chain_a)
            post_binder = get_heavy_atoms(chain_b)
            
            post_t_coords = np.array([a.coord for a in post_target])
            post_b_coords = np.array([a.coord for a in post_binder])
            
            post_distances = cdist(post_t_coords, post_b_coords)
            post_min_dist = float(np.min(post_distances))
            post_clashes_25 = int(np.sum(post_distances < 2.5))
            post_clashes_19 = int(np.sum(post_distances < 1.9))
            post_rho = int(np.sum(post_distances < 4.5))
            
            print(f"   Min distance:      {post_min_dist:.2f} Å")
            print(f"   Clashes (< 2.5Å):  {post_clashes_25}")
            print(f"   Clashes (< 1.9Å):  {post_clashes_19}")
            print(f"   Contact density:   {post_rho}")
        else:
            print("   ERROR: Could not identify chains in minimized structure")
            post_min_dist = pre_min_dist
            post_clashes_25 = pre_clashes_25
            post_clashes_19 = pre_clashes_19
            post_rho = pre_rho
    else:
        post_min_dist = pre_min_dist
        post_clashes_25 = pre_clashes_25
        post_clashes_19 = pre_clashes_19
        post_rho = pre_rho
    
    # --- VERDICT ---
    print("\n" + "="*60)
    print("VERDICT")
    print("="*60)
    
    print(f"\n{'Metric':<25} {'Before':<12} {'After':<12} {'Change'}")
    print("-" * 55)
    print(f"{'Min distance (Å)':<25} {pre_min_dist:<12.2f} {post_min_dist:<12.2f} {post_min_dist - pre_min_dist:+.2f}")
    print(f"{'Clashes (< 2.5Å)':<25} {pre_clashes_25:<12} {post_clashes_25:<12} {post_clashes_25 - pre_clashes_25:+d}")
    print(f"{'Clashes (< 1.9Å)':<25} {pre_clashes_19:<12} {post_clashes_19:<12} {post_clashes_19 - pre_clashes_19:+d}")
    print(f"{'Contact density':<25} {pre_rho:<12} {post_rho:<12} {post_rho - pre_rho:+d}")
    
    print("")
    
    if post_clashes_25 == 0:
        print("✅ VERDICT: SOVEREIGN_PASS (2.5Å threshold)")
        print("   All clashes resolved. Lead validated.")
        verdict = "SOVEREIGN_PASS"
    elif post_clashes_19 == 0:
        print("⚠️  VERDICT: MARGINAL_PASS (1.9Å threshold)")
        print("   No severe clashes, but minor overlaps remain.")
        verdict = "MARGINAL_PASS"
    else:
        print("❌ VERDICT: PHYSICAL_VETO")
        print("   Clashes persist after minimization.")
        verdict = "PHYSICAL_VETO"
    
    # Return results for programmatic use
    return {
        "pre": {
            "min_distance": pre_min_dist,
            "clashes_25": pre_clashes_25,
            "clashes_19": pre_clashes_19,
            "rho": pre_rho
        },
        "post": {
            "min_distance": post_min_dist,
            "clashes_25": post_clashes_25,
            "clashes_19": post_clashes_19,
            "rho": post_rho
        },
        "verdict": verdict
    }


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python kinetic_rescue.py <structure.cif> [output.pdb]")
        sys.exit(1)
    
    cif_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else "minimized_structure.pdb"
    
    results = run_rescue(cif_path, output_path)
    
    # Save results
    results_file = Path(cif_path).parent / "rescue_results.json"
    with open(results_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {results_file}")
