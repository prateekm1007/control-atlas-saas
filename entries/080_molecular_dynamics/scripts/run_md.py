#!/usr/bin/env python3
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

INPUT_PDB = "../inputs/system_solvated.pdb"
OUTPUT_PREFIX = "../outputs/entry080"

print("🚀 [Entry 080] Starting Production MD...")

pdb = PDBFile(INPUT_PDB)
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometer,
    constraints=HBonds
)

integrator = LangevinIntegrator(
    300*kelvin,
    1/picosecond,
    0.002*picoseconds
)

platform = Platform.getPlatformByName("CPU")
simulation = Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

print("🔧 Minimizing energy...")
simulation.minimizeEnergy()

simulation.reporters.append(
    StateDataReporter(
        stdout, 1000,
        step=True,
        potentialEnergy=True,
        temperature=True,
        progress=True,
        remainingTime=True,
        speed=True,
        totalSteps=5_000_000
    )
)

simulation.reporters.append(DCDReporter(f"{OUTPUT_PREFIX}.dcd", 5000))
simulation.reporters.append(
    StateDataReporter(
        f"{OUTPUT_PREFIX}_md_log.csv", 1000,
        step=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True
    )
)

print("▶️ Running 10 ns MD...")
simulation.step(5_000_000)

state = simulation.context.getState(getPositions=True)
with open(f"{OUTPUT_PREFIX}_final.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, state.getPositions(), f)

print("✅ MD complete.")
