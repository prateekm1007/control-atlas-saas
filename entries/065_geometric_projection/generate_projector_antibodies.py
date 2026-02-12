#!/usr/bin/env python3
import json
import random
import argparse

# CONFIGURATION
WARHEADS = ["YRY", "RWR", "YKY", "RFR"]
STEM_LEFT_POOL = ["PVR", "YVR", "AVR", "TVR"]
STEM_RIGHT_POOL = ["DVP", "DIP", "DYP"]
ANCHOR_SUFFIX = "FDY" 
FORBIDDEN = ["NG", "NS", "DG"]

def generate_projector_cdr3():
    warhead = random.choice(WARHEADS)
    stem_L = random.choice(STEM_LEFT_POOL)
    stem_R = random.choice(STEM_RIGHT_POOL)
    
    candidate = stem_L + warhead + stem_R + ANCHOR_SUFFIX
    
    if any(m in candidate for m in FORBIDDEN): return generate_projector_cdr3()
    if "C" in candidate: return generate_projector_cdr3()
    
    return candidate, warhead

def main():
    print(f"🧬 [Entry 065] Generating Geometric Projector Antibodies (v5)...")
    vh_prefix = "EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYAMSWVRQAPGKGLEWVAVISYDGSTYYSADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCSR"
    vh_suffix = "WGQGTLVTVSS"
    vl_seq = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"

    candidates = []
    for i in range(5):
        cid = f"CAND_v5_{i:04d}"
        cdr3, warhead = generate_projector_cdr3()
        candidates.append({
            "id": cid,
            "heavy_chain": vh_prefix + cdr3 + vh_suffix,
            "light_chain": vl_seq,
            "cdr3": cdr3,
            "warhead": warhead
        })
        print(f"  [+] {cid}: {cdr3} (Rigid Stem)")

    with open("candidates_v5_projector.json", "w") as f:
        json.dump(candidates, f, indent=2)
    print(f"\n✅ Generated 5 Projector Candidates → candidates_v5_projector.json")

if __name__ == "__main__":
    main()
