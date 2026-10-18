import sys
from backend_chai import ChaiBackend
from tools.native_audit import SovereignJudge

def run_cycle(target_fasta, binder_seq, name):
    print(f"🚀 INITIALIZING CYCLE: {name}")
    
    # 1. Inference (Witness A)
    engine = ChaiBackend()
    structures = engine.run_constrained_docking(target_fasta, seed=42)
    
    # 2. Audit (The Executioner)
    judge = SovereignJudge()
    for cif in structures:
        passed, dist, reason = judge.audit(str(cif))
        if passed:
            print(f"💎 SUCCESS: {name} Cleared at {dist}A.")
        else:
            print(f"❌ VETO: {name} Failed - {reason}")

if __name__ == "__main__":
    # Example command line usage
    pass
