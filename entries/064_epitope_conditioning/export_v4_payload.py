#!/usr/bin/env python3
import json

CANDIDATES_FILE = "candidates_v4_motif.json"
TARGET_SEQ = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQ" # KRAS G12D

def main():
    with open(CANDIDATES_FILE) as f:
        candidates = json.load(f)
    
    print("\n" + "="*20 + " COPY BELOW THIS LINE " + "="*20)
    print("JOBS_PAYLOAD = [")
    for c in candidates:
        entry = {
            "id": f"{c['id']}_vs_KRAS",
            "antibody_seq": c["heavy_chain"] + c["light_chain"],
            "antigen_seq": TARGET_SEQ
        }
        print(f"    {json.dumps(entry)},")
    print("]")
    print("="*20 + " COPY ABOVE THIS LINE " + "="*20 + "\n")

if __name__ == "__main__":
    main()
