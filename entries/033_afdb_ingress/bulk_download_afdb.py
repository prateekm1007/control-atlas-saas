#!/usr/bin/env python3
""" Bulk download AFDB CIFs via FTP (2025 production method) """
import argparse
import subprocess
from pathlib import Path

FTP_BASE = "ftp://ftp.ebi.ac.uk/pub/databases/alphafold/latest"

def download_targets(uniprot_ids: list, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    
    lftp_commands = [
        "open " + FTP_BASE,
        "cd UP000005640_9606_HUMAN",  # Human proteome
        "lcd " + str(output_dir),
    ]
    
    mget_files = []
    for uid in uniprot_ids:
        mget_files.append(f"AF-{uid.upper()}-F1-model_v4.cif.gz")
    
    lftp_commands.append("mget " + " ".join(mget_files))
    lftp_commands.append("bye")
    
    cmd = ["lftp", "-c", "\n".join(lftp_commands)]
    print(f"[*] Downloading {len(uniprot_ids)} CIFs via FTP...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print("[!] FTP failed:", result.stderr)
    else:
        print("[+] Download complete")
    
    # Decompress
    print("[*] Decompressing...")
    for gz in output_dir.glob("*.cif.gz"):
        subprocess.run(["gunzip", str(gz)], check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ids", required=True, help="File with UniProt IDs")
    parser.add_argument("--out-dir", default="./data")
    args = parser.parse_args()
    
    with open(args.ids) as f:
        ids = [line.strip().upper() for line in f if line.strip()]
    
    download_targets(ids, Path(args.out_dir))
