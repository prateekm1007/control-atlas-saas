#!/usr/bin/env python3
"""
Entry 052 — Kaggle Persistence Adapter
Syncs local results to a Kaggle Dataset to prevent data loss.
"""

import os
import sys
from pathlib import Path

def sync_results(local_dir, dataset_id):
    """
    Uploads contents of local_dir to a Kaggle Dataset.
    Requires kaggle.json token.
    """
    print(f"[*] Syncing {local_dir} to {dataset_id}...")
    
    # Check for API token
    if not Path("~/.kaggle/kaggle.json").expanduser().exists():
        print("[!] Kaggle API token not found.")
        return

    # Using Kaggle CLI (assumed installed in environment)
    # 1. Init metadata if missing
    if not (Path(local_dir) / "dataset-metadata.json").exists():
        os.system(f"kaggle datasets init -p {local_dir}")
        
    # 2. Upload (version update)
    cmd = f"kaggle datasets version -p {local_dir} -m 'Auto-sync from RFAA run'"
    os.system(cmd)
    print("[+] Sync initiated.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python kaggle_sync.py <dir> <dataset_slug>")
    else:
        sync_results(sys.argv[1], sys.argv[2])
