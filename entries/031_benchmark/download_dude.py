#!/usr/bin/env python3
"""
Download DUD-E benchmark data for selected targets.
"""

import os
import sys
import urllib.request
import tarfile
from pathlib import Path

DUDE_BASE_URL = "http://dude.docking.org/targets"

TARGETS = ["egfr", "braf", "cdk2", "src", "vgfr2"]

def download_target(target: str, output_dir: Path):
    """Download and extract DUD-E data for a target."""
    output_dir = Path(output_dir)
    target_dir = output_dir / target
    
    if target_dir.exists():
        print(f"[*] {target} already downloaded")
        return True
    
    url = f"{DUDE_BASE_URL}/{target}/{target}.tar.gz"
    tar_path = output_dir / f"{target}.tar.gz"
    
    print(f"[*] Downloading {target}...")
    try:
        urllib.request.urlretrieve(url, tar_path)
    except Exception as e:
        print(f"[!] Failed to download {target}: {e}")
        return False
    
    print(f"[*] Extracting {target}...")
    try:
        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extractall(output_dir)
        os.remove(tar_path)
    except Exception as e:
        print(f"[!] Failed to extract {target}: {e}")
        return False
    
    print(f"[+] {target} ready")
    return True


def main():
    output_dir = Path(__file__).parent / "data"
    output_dir.mkdir(exist_ok=True)
    
    print("="*60)
    print("DUD-E Benchmark Data Download")
    print("="*60)
    
    success = 0
    for target in TARGETS:
        if download_target(target, output_dir):
            success += 1
    
    print(f"\n[+] Downloaded {success}/{len(TARGETS)} targets")
    

if __name__ == "__main__":
    main()
