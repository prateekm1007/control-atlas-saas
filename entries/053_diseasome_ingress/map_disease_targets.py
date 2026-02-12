#!/usr/bin/env python3
"""
Entry 053 — Disease Mapper
Maps DisGeNET Genes -> UniProt IDs for Control Atlas Ingestion.
"""

import csv
import json
import requests
from pathlib import Path

# Mock Mapping (Gene Symbol -> UniProt)
# In prod, use UniProt API or mapping file
GENE_MAP = {
    "KRAS": "P01116",
    "TP53": "P04637",
    "EGFR": "P00533",
    "BRAF": "P15056",
    "GAPDH": "P04406"
}

def map_targets(tsv_path, output_dir):
    disease_map = {} # {DiseaseName: [UniProtIDs]}
    
    with open(tsv_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row["geneSymbol"]
            disease = row["diseaseName"]
            score = float(row["score"])
            
            # Filter by confidence
            if score < 0.4: continue
            
            uid = GENE_MAP.get(gene)
            if uid:
                if disease not in disease_map:
                    disease_map[disease] = set()
                disease_map[disease].add(uid)
                
    # Save target lists
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = {}
    
    for disease, uids in disease_map.items():
        # Create a target list file for Entry 033/034 to consume
        clean_name = disease.replace(" ", "_").lower()
        fname = f"{clean_name}_targets.txt"
        fpath = output_dir / fname
        
        with open(fpath, "w") as f:
            for uid in uids:
                f.write(f"{uid}\n")
        
        manifest[disease] = {"file": str(fpath), "count": len(uids)}
        print(f"[+] Mapped {disease}: {len(uids)} targets -> {fname}")
        
    with open(output_dir / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)

if __name__ == "__main__":
    map_targets(Path("data/curated_gene_disease_associations.tsv"), Path("disease_lists"))
