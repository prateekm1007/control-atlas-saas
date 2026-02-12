import MDAnalysis as mda
import numpy as np

u = mda.Universe("../outputs/starting_state.pdb", "../outputs/production_traj.dcd")
ligand = u.select_atoms("resname ARG")
target = u.select_atoms("resname ASP")

u.trajectory[0] # JUMP TO START
dists = np.linalg.norm(ligand.positions[:, None, :] - target.positions[None, :, :], axis=-1)
min_dist = dists.min()

print(f"\n🔎 FRAME-0 MIN HANDSHAKE DISTANCE: {min_dist:.2f} Å")
if min_dist > 40.0:
    print("🟨 VERDICT: MODELING HALLUCINATION. The complex never started in the pocket.")
else:
    print("🟥 VERDICT: PHYSICAL REPULSION. It started close but the forcefield pushed it away.")
