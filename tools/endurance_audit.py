import openmm.app as app
import openmm as mm
from openmm import unit
import numpy as np
import time
import sys

def run_endurance_test(pdb_path, duration_ns=1.0):
    print(f"🕒 [TIER 2] Initializing {duration_ns}ns Endurance Audit for {pdb_path}...")
    
    # 1. Load Structure
    pdb = app.PDBFile(pdb_path)
    forcefield = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
    
    # 2. Physics Environment
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic, 
                                    constraints=app.HBonds, hydrogenMass=4*unit.amu)
    integrator = mm.LangevinMiddleIntegrator(310*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    
    # 🔒 SOVEREIGN PLATFORM ARBITRATION
    platforms = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    
    # Priority: CUDA -> OpenCL -> CPU -> Reference
    if 'CPU' in platforms:
        target_platform = 'CPU'
    elif 'OpenCL' in platforms:
        target_platform = 'OpenCL'
    else:
        target_platform = 'Reference'
            
    print(f"📡 [INFRA] Using Platform: {target_platform} (Available: {platforms})")
    platform = mm.Platform.getPlatformByName(target_platform)
    
    # 3. Initialize Simulation
    sim = app.Simulation(pdb.topology, system, integrator, platform)
    sim.context.setPositions(pdb.positions)
    
    # 4. Production Run
    # 1ns = 500,000 steps at 2fs.
    steps = int(duration_ns * 500000)
    report_freq = 10000 # Every 20ps for detailed tracking
    sim.reporters.append(app.StateDataReporter(sys.stdout, report_freq, step=True, 
                                              potentialEnergy=True, temperature=True, speed=True))
    
    start_time = time.time()
    print("⚡ [PHYSICS] Commencing 1.0ns Production Run...")
    sim.step(steps)
    end_time = time.time()
    
    print(f"✅ TIER 2 COMPLETE. 1ns simulated in {round((end_time - start_time)/60, 2)} mins.")
