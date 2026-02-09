import os
import json
import pandas as pd

DATA_DIR = "entries/034_kinome_atlas/data/"
MANIFEST_PATH = "ledger/kinome_manifest.json"

def build_manifest():
    targets = []
    if not os.path.exists(DATA_DIR):
        print(f"❌ Error: {DATA_DIR} not found.")
        return

    for file in os.listdir(DATA_DIR):
        if file.endswith((".cif", ".pdb")):
            targets.append({
                "id": file.split('.')[0],
                "path": os.path.join(DATA_DIR, file),
                "type": "kinase",
                "priority": 1
            })
    
    with open(MANIFEST_PATH, "w") as f:
        json.dump(targets, f, indent=2)
    print(f"✅ Manifest built with {len(targets)} targets at {MANIFEST_PATH}")

if __name__ == "__main__":
    build_manifest()
