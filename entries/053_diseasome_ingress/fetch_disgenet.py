#!/usr/bin/env python3
"""
Entry 053 — DisGeNET Ingress
Fetches Curated Gene-Disease Associations.
"""

import os
import requests
import gzip
import shutil
from pathlib import Path

# URL for Curated GDAs (Free for academic use)
DISGENET_URL = "https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz"

def fetch_data(output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    target_file = output_dir / "curated_gene_disease_associations.tsv"
    gz_file = target_file.with_suffix(".tsv.gz")
    
    # Mock Data Creation (since we can't download 100MB in this shell)
    # In production, uncomment the download logic below
    
    print("[*] Creating Mock DisGeNET Data (Demo Mode)...")
    with open(target_file, "w") as f:
        f.write("geneId\tgeneSymbol\tdiseaseId\tdiseaseName\tscore\tEI\n")
        # KRAS - Lung Cancer
        f.write("3845\tKRAS\tC0007131\tLung Neoplasms\t0.9\t1.0\n")
        # TP53 - Li-Fraumeni
        f.write("7157\tTP53\tC0023418\tLi-Fraumeni Syndrome\t1.0\t1.0\n")
        # EGFR - Glioblastoma
        f.write("1956\tEGFR\tC0017636\tGlioblastoma\t0.8\t1.0\n")
        # BRAF - Melanoma
        f.write("673\tBRAF\tC0025202\tMelanoma\t0.9\t1.0\n")
        # GAPDH - Alzheimer's (Low relevance link for test)
        f.write("2597\tGAPDH\tC0002395\tAlzheimer Disease\t0.3\t0.5\n")
    
    print(f"[+] Data ready at {target_file}")
    return target_file

    # --- Production Download Logic ---
    # print(f"[*] Downloading DisGeNET from {DISGENET_URL}...")
    # try:
    #     r = requests.get(DISGENET_URL, stream=True)
    #     if r.status_code == 200:
    #         with open(gz_file, 'wb') as f:
    #             f.write(r.raw.read())
    #         with gzip.open(gz_file, 'rb') as f_in, open(target_file, 'wb') as f_out:
    #             shutil.copyfileobj(f_in, f_out)
    #         print("[+] Download complete.")
    #         return target_file
    #     else:
    #         print(f"[!] Download failed: {r.status_code}")
    #         return None
    # except Exception as e:
    #     print(f"[!] Error: {e}")
    #     return None

if __name__ == "__main__":
    fetch_data(Path("data"))
