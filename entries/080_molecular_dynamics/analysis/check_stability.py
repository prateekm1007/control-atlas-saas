import MDAnalysis as mda
from MDAnalysis.analysis import rms
import pandas as pd
import matplotlib.pyplot as plt

# CONFIGURATION
LOG_FILE = "../outputs/md_log.csv"
TRAJ_FILE = "../outputs/production_traj.dcd"
TOPOLOGY = "../inputs/prepared_system.pdb" # We will download this from Kaggle too

def analyze():
    print("🧪 [Entry 080] Analyzing Physical Evidence...")
    
    # 1. Check Energy Stability (from CSV)
    df = pd.read_csv(LOG_FILE)
    plt.figure(figsize=(10,5))
    plt.plot(df['#"Step"'], df['Potential Energy (kJ/mole)'])
    plt.title("Potential Energy vs Time (Thermal Stability)")
    plt.savefig("energy_plot.png")
    print("  -> Energy plot generated.")

    # 2. Check Structural Stability (RMSD)
    # This requires the DCD and PDB (Topology)
    # (Commented out until files exist)
    # u = mda.Universe(TOPOLOGY, TRAJ_FILE)
    # R = rms.RMSD(u, select="backbone").run()
    # plt.figure()
    # plt.plot(R.results.rmsd[:,2])
    # plt.title("Backbone RMSD (Structure Integrity)")
    # plt.savefig("rmsd_plot.png")

if __name__ == "__main__":
    analyze()
