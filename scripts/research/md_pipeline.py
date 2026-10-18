import openmm as mm
import openmm.app as app
from openmm import unit
from pdbfixer import PDBFixer

class KineticFalsifier:
    def __init__(self, pdb_path):
        self.pdb_path = pdb_path
        
    def run_micro_burst(self, output_dcd="burst.dcd", steps=250000):
        """Executes a 500ps forensic burst in implicit solvent (OBC2)."""
        fixer = PDBFixer(self.pdb_path)
        fixer.findMissingResidues(); fixer.findMissingAtoms(); fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.4)
        
        forcefield = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
        system = forcefield.createSystem(fixer.topology, nonbondedMethod=app.CutoffNonPeriodic, 
                                        constraints=app.HBonds, hydrogenMass=4*unit.amu)
        
        integrator = mm.LangevinMiddleIntegrator(310*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
        simulation = app.Simulation(fixer.topology, system, integrator)
        simulation.context.setPositions(fixer.positions)
        
        print("⚡ Minimizing Energy...")
        simulation.minimizeEnergy()
        
        print(f"⚡ Running {steps} steps of Kinetic Falsification...")
        simulation.reporters.append(app.DCDReporter(output_dcd, 1000))
        simulation.step(steps)
        
        state = simulation.context.getState(getEnergy=True)
        return state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
