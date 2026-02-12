#!/usr/bin/env python3
"""
Entry 042 — Knowledge Updater
Orchestrates data refresh for the Atlas.
"""

import sys
import json
import glob
from pathlib import Path
from datetime import datetime

# Import Sources
from sources import afdb_source, depmap_source

# Import Ingestion Logic
BASE = Path(__file__).resolve().parents[2] / "entries"
sys.path.insert(0, str(BASE / "033_afdb_ingress")) # Re-use Ingress logic
from ingest_proteome import process_protein

ATLAS_DIR = Path(__file__).resolve().parents[2] / "library/atlas_index"

def update_atlas():
    print("="*60)
    print("Control Atlas — Knowledge Updater")
    print("="*60)
    
    # 1. Scan Atlas for Outdated Structures
    print("\n[*] Checking Structure Versions (AlphaFold DB)...")
    files = glob.glob(str(ATLAS_DIR / "*.json"))
    
    updates_needed = []
    
    for f in files:
        with open(f, "r") as json_file:
            entry = json.load(json_file)
            
        uid = entry.get("uniprot_id")
        current_ver = entry.get("model_version", "v4") # Default to v4 if missing
        
        latest_ver = afdb_source.get_latest_version(uid)
        
        if latest_ver and latest_ver > current_ver:
            print(f"    UPDATE FOUND: {uid} ({current_ver} -> {latest_ver})")
            updates_needed.append(uid)
        else:
            # print(f"    Current: {uid} ({current_ver})")
            pass
            
    # 2. Execute Updates (Re-Ingest)
    if updates_needed:
        print(f"\n[*] Applying {len(updates_needed)} Structure Updates...")
        for uid in updates_needed:
            print(f" -> Re-ingesting {uid}...")
            # We re-use the Entry 033 logic to fetch, filter, detect, and index
            # This ensures physics gates are applied to the NEW structure
            res = process_protein(uid, Path("data")) # Temp data dir
            print(f"    Result: {res['status']}")
    else:
        print("\n[+] All structures are up to date.")

    # 3. Biology Updates (Mocked Demo)
    print("\n[*] Checking Biology Context (DepMap)...")
    # In a real run, we'd loop through all targets. Demo: EGFR.
    new_bio = depmap_source.check_update("EGFR", "23Q4")
    if new_bio:
        print(f"    UPDATE FOUND: EGFR Essentiality (23Q4 -> {new_bio['version']})")
        # Here we would update the Biology Knowledge Graph JSON
    
    print("\n[+] Update Cycle Complete.")

if __name__ == "__main__":
    update_atlas()
