import sys
import os
sys.path.append(os.getcwd())

from tools.asym_generator import AsymBipodGenerator
from backend_chai import ChaiBackend

# PD-L1 IgV Domain (Full sequence for accurate folding)
PDL1_TARGET = "FTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRARLLKDQLSLGNAALQITDVKLQDAGVYRCMISYGGADYKRITVKVNAPY"

def main():
    gen = AsymBipodGenerator()
    engine = ChaiBackend()
    
    warheads = ["YWPTG", "WPTGY", "RYWPTG", "YWPTGR"]
    
    print("\n🚀 WSL PRE-COMPUTE: PREPARING KAGLLE BATCH")
    print("="*50)

    for w in warheads:
        seq = gen.generate_candidate(w)
        name = f"ASYM_{w}"
        fasta = engine.prepare_fasta(PDL1_TARGET, seq, name)
        print(f"✅ GENERATED: {fasta}")

if __name__ == "__main__":
    main()
