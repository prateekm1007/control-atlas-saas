#!/usr/bin/env python3
import json
import random
import argparse

# CONFIGURATION
EPITOPE_FILE = "kras_epitopes.json"
ANCHOR_SUFFIX = "FDY" # Keep the v3 structural anchor
FORBIDDEN = ["NG", "NS", "DG"]

def generate_motif_cdr3(motifs, target_len=13):
    """
    Injects a Switch-II killer motif into the loop.
    Structure: [Stem] - [Warhead] - [Stem] - [Anchor]
    """
    warhead = random.choice(motifs)
    
    # Calculate remaining space
    remaining = target_len - len(ANCHOR_SUFFIX) - len(warhead)
    if remaining < 2: remaining = 2 # Min stem
    
    # Flanking residues (Flexible/Neutral)
    flank_pool = "GSTSA" 
    
    while True:
        left_flank = "".join(random.choice(flank_pool) for _ in range(remaining // 2))
        right_flank = "".join(random.choice(flank_pool) for _ in range(remaining - len(left_flank)))
        
        # Assemble
        candidate = left_flank + warhead + right_flank + ANCHOR_SUFFIX
        
        # Safety Checks
        if any(m in candidate for m in FORBIDDEN): continue
        if "C" in candidate: continue
        
        return candidate, warhead

def main():
    with open(EPITOPE_FILE) as f:
        data = json.load(f)
    
    # Extract Motifs
    motifs = data["KRAS"]["switch_II"]["complementary_motifs"]
    print(f"🧬 [Entry 064] Loaded Warheads: {motifs}")

    # Framework (Trastuzumab)
    vh_prefix = "EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYAMSWVRQAPGKGLEWVAVISYDGSTYYSADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCSR"
    vh_suffix = "WGQGTLVTVSS"
    vl_seq = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"

    candidates = []
    for i in range(5):
        cid = f"CAND_v4_{i:04d}"
        cdr3, motif_used = generate_motif_cdr3(motifs)
        
        candidates.append({
            "id": cid,
            "heavy_chain": vh_prefix + cdr3 + vh_suffix,
            "light_chain": vl_seq,
            "cdr3": cdr3,
            "motif": motif_used
        })
        print(f"  [+] {cid}: {cdr3} (Warhead: {motif_used})")

    with open("candidates_v4_motif.json", "w") as f:
        json.dump(candidates, f, indent=2)
    print(f"\n✅ Generated 5 Motif-Enhanced Candidates → candidates_v4_motif.json")

if __name__ == "__main__":
    main()
