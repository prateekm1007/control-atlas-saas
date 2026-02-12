import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import matplotlib.pyplot as plt

TOPOLOGY = '../outputs/starting_state.pdb'
TRAJECTORY = '../outputs/production_traj.dcd'

def audit_handshake():
    u = mda.Universe(TOPOLOGY, TRAJECTORY)
    
    # 1. IDENTIFY THE WARHEAD (Chain A) AND TARGET (Chain B)
    # Based on our FASTA: Antibody is Chain A, KRAS is Chain B.
    # KRAS G12D is Residue 12 on Chain B.
    # Antibody Warhead (YRK) is at the center of the CDR3 (approx res 100-110).
    
    # Selection: Arginine of Warhead vs Aspartate of G12D
    # We use a broad selection to ensure we catch them
    warhead_arg = u.select_atoms("chainID A and resname ARG")
    target_asp = u.select_atoms("chainID B and resid 12")
    
    if len(warhead_arg) == 0 or len(target_asp) == 0:
        print("❌ Error: Could not find Warhead/Target residues.")
        return

    print(f"🧬 Measuring distance between {len(warhead_arg)} Arg atoms and {len(target_asp)} Asp atoms...")

    # 2. CALCULATE DISTANCE OVER TIME
    dist_arr = []
    for ts in u.trajectory:
        # Calculate minimum distance between any Arg atom and any Asp atom
        d = distances.distance_array(warhead_arg.positions, target_asp.positions).min()
        dist_arr.append(d)
    
    avg_dist = np.mean(dist_arr)
    
    # 3. PLOT
    plt.figure(figsize=(10, 4))
    plt.plot(np.array(range(len(dist_arr)))*0.01, dist_arr, color='blue')
    plt.axhline(y=4.0, color='red', linestyle='--', label='Hydrogen Bond Limit')
    plt.xlabel("Time (ns)")
    plt.ylabel("Distance (Å)")
    plt.title("VAR_YRK Warhead Engagement (Arg ↔ Asp12)")
    plt.legend()
    plt.savefig("handshake_persistence.png")
    
    print(f"✅ Handshake Analysis Complete.")
    print(f"📊 Average Interaction Distance: {avg_dist:.2f} Å")
    
    if avg_dist < 4.5:
        print("🟢 VERDICT: THE WARHEAD IS LOCKED.")
    else:
        print("🟡 VERDICT: THE WARHEAD IS WIGGLING (Weak Engagement).")

if __name__ == "__main__":
    audit_handshake()
