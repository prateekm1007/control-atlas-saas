#!/usr/bin/env python3
"""
Batch ingestion of COSMIC driver pockets
"""

import os
import csv
import subprocess

CATALOG = os.path.expanduser("~/control-atlas/library/pocket_catalog")
DRIVERS = os.path.join(CATALOG, "cosmic_drivers_v1.csv")
INGEST = os.path.join(CATALOG, "ingest_target.py")

def batch_ingest():
    print("=== BATCH POCKET INGESTION ===")
    
    with open(DRIVERS, 'r') as f:
        reader = csv.DictReader(f)
        targets = list(reader)
    
    print(f"Processing {len(targets)} targets...")
    
    success = 0
    failed = 0
    
    for t in targets:
        cmd = [
            "python3", INGEST,
            "--target", t["target"],
            "--mutation", t["mutation"],
            "--domain", t["domain"],
            "--uniprot", t["uniprot"],
            "--pdbs", t["pdbs"],
            "--lining", t["lining"]
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            print(f"  [OK] {t['target']}_{t['mutation']}")
            success += 1
        except subprocess.CalledProcessError as e:
            print(f"  [FAIL] {t['target']}_{t['mutation']}: {e}")
            failed += 1
    
    print(f"\n=== COMPLETE ===")
    print(f"Success: {success}")
    print(f"Failed: {failed}")

if __name__ == "__main__":
    batch_ingest()
