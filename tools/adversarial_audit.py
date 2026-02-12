from Bio.PDB import MMCIFParser
from scipy.spatial.distance import cdist
import numpy as np

def run_adversarial_audit(path, clash_limit=2.0):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("audit", path)[0]
    
    # Identify chains
    chains = sorted(list(structure.get_chains()), key=lambda c: len(list(c.get_atoms())), reverse=True)
    target, binder = chains[0], chains[1]
    
    t_coords = np.array([a.coord for a in target.get_atoms() if a.element != 'H'])
    b_coords = np.array([a.coord for a in binder.get_atoms() if a.element != 'H'])
    
    dists = cdist(t_coords, b_coords)
    min_dist = np.min(dists)
    clashes = np.sum(dists < clash_limit)
    rho = np.sum(dists < 4.5)
    
    print(f"📊 ADVERSARIAL RESULT for {path}:")
    print(f"   Min Distance: {min_dist:.2f} Å")
    print(f"   Clashes (<{clash_limit}Å): {clashes}")
    print(f"   Contact Density (Rho): {rho}")
    
    if min_dist < clash_limit:
        print("❌ STATUS: CLASH_VETO")
    else:
        print("✅ STATUS: SOVEREIGN_PASS")

if __name__ == "__main__":
    import sys
    run_adversarial_audit(sys.argv[1])
