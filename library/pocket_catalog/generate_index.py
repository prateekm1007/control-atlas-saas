#!/usr/bin/env python3
"""
Generate catalog index from ingested pockets
"""

import os
import json

CATALOG = os.path.expanduser("~/control-atlas/library/pocket_catalog")
OUTPUT = os.path.join(CATALOG, "CATALOG_INDEX.json")

def generate_index():
    index = {
        "version": "1.0",
        "targets": []
    }
    
    for entry in os.listdir(CATALOG):
        path = os.path.join(CATALOG, entry, "pocket_frame.json")
        if os.path.exists(path):
            with open(path, 'r') as f:
                frame = json.load(f)
                index["targets"].append({
                    "id": f"{frame['identity']['target']}_{frame['identity']['mutation']}",
                    "domain": frame["pocket"]["classification"],
                    "state": frame["pocket"]["reference_state"],
                    "lining_count": len(frame["frame"]["lining_residues"])
                })
    
    with open(OUTPUT, 'w') as f:
        json.dump(index, f, indent=2)
    
    print(f"Index generated: {len(index['targets'])} targets")
    print(f"Saved to {OUTPUT}")

if __name__ == "__main__":
    generate_index()
