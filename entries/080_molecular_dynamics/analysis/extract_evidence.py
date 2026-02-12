import MDAnalysis as mda
from MDAnalysis.analysis import rms
import pandas as pd
import matplotlib.pyplot as plt
import os

TOPOLOGY = '../outputs/starting_state.pdb'
TRAJECTORY = '../outputs/production_traj.dcd'

def run_harvest():
    print("🧪 [Entry 081] Extracting Structural Stability (RMSD)...")
    
    if not os.path.exists(TRAJECTORY):
        print(f"❌ Error: Trajectory not found at {TRAJECTORY}")
        return

    # Load Universe
    u = mda.Universe(TOPOLOGY, TRAJECTORY)
    print(f"  -> System loaded: {u.atoms.n_atoms} atoms.")
    
    # 1. ALIGN & CALCULATE RMSD
    # We select 'backbone' for alignment and 'protein' for the complex RMSD
    R = rms.RMSD(u, u, select="backbone", groupselections=["protein"])
    R.run()
    
    # 2. SAVE DATA
    rmsd_df = pd.DataFrame(R.results.rmsd, columns=['Frame', 'Time', 'Backbone', 'Complex'])
    rmsd_df.to_csv("rmsd_data.csv", index=False)
    
    # 3. PLOT
    plt.figure(figsize=(10, 5))
    plt.plot(rmsd_df['Time']/1000, rmsd_df['Complex'], color='crimson')
    plt.xlabel("Time (ns)")
    plt.ylabel("RMSD (Å)")
    plt.title("VAR_YRK Persistence Trace: 10ns Production")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.savefig("rmsd_persistence.png")
    
    avg_rmsd = rmsd_df['Complex'].mean()
    print(f"\n✅ DATA HARVESTED.")
    print(f"📊 Average RMSD: {avg_rmsd:.2f} Å")
    print(f"🖼️  Plot saved as: rmsd_persistence.png")

if __name__ == "__main__":
    run_harvest()
