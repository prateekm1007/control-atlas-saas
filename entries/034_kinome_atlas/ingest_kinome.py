#!/usr/bin/env python3
"""
Entry 034 — Human Kinome Atlas Ingestion
FIXED: Uses gemmi to extract pLDDT (authoritative, no string parsing).
"""

import argparse
import json
import sys
from pathlib import Path
from time import time

BASE = Path(__file__).resolve().parents[2] / "entries"
sys.path.insert(0, str(BASE / "033_afdb_ingress"))
sys.path.insert(0, str(BASE / "027_pocket_detection"))

from fetch_afdb import fetch_structure
from pocket_detector import PocketDetector

ATLAS_DIR = Path(__file__).resolve().parents[2] / "library/atlas_index"
ATLAS_DIR.mkdir(parents=True, exist_ok=True)

def get_plddt_from_cif(cif_path):
    """Extract global pLDDT from AlphaFold CIF using gemmi (authoritative)."""
    import gemmi
    
    structure = gemmi.read_structure(str(cif_path))
    scores = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    b = atom.b_iso
                    if 0 < b <= 100:
                        scores.append(b)
    
    return sum(scores) / len(scores) if scores else 0.0

def process_kinase(uniprot_id, data_dir, min_plddt=70.0):
    paths = fetch_structure(uniprot_id, data_dir)
    if not paths:
        return {"id": uniprot_id, "status": "DOWNLOAD_FAILED"}
    
    cif_path = paths["cif"]
    pdb_path = paths["pdb"]
    
    raw_plddt = get_plddt_from_cif(cif_path)
    
    if raw_plddt < min_plddt:
        return {"id": uniprot_id, "status": "LOW_CONFIDENCE", "plddt": round(raw_plddt, 1)}
    
    conf = raw_plddt / 100.0
    
    detector = PocketDetector(max_pockets=5)
    result = detector.detect(pdb_path, structure_confidence=conf)
    
    if result["status"] != "SUCCESS":
        return {"id": uniprot_id, "status": "DETECTION_ERROR"}
    
    pockets = [p for p in result["pockets"] if p["status"] in ("VALIDATED", "CANDIDATE")]
    
    if not pockets:
        return {"id": uniprot_id, "status": "NO_POCKETS", "plddt": round(raw_plddt, 1)}
    
    entry = {
        "uniprot_id": uniprot_id,
        "protein_family": "kinase",
        "plddt_global": round(raw_plddt, 1),
        "plddt_source": "CIF_gemmi",
        "pockets": pockets,
        "source": "AlphaFoldDB_v6"
    }
    
    with open(ATLAS_DIR / f"{uniprot_id}.json", "w") as f:
        json.dump(entry, f, indent=2)
    
    return {"id": uniprot_id, "status": "INDEXED", "pockets": len(pockets), "plddt": round(raw_plddt, 1)}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--list", required=True)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()
    
    data_dir = Path(__file__).parent / "data"
    data_dir.mkdir(exist_ok=True)
    
    with open(args.list) as f:
        kinases = [line.strip() for line in f if line.strip()]
    
    if args.limit:
        kinases = kinases[:args.limit]
    
    print(f"[*] Processing {len(kinases)} kinases (gemmi pLDDT)...")
    print("="*60)
    
    t0 = time()
    stats = {"INDEXED": 0, "NO_POCKETS": 0, "LOW_CONFIDENCE": 0, "DOWNLOAD_FAILED": 0, "DETECTION_ERROR": 0}
    
    for i, uid in enumerate(kinases, 1):
        print(f"[{i}/{len(kinases)}] {uid}...", end=" ", flush=True)
        res = process_kinase(uid, data_dir)
        status = res["status"]
        stats[status] = stats.get(status, 0) + 1
        
        if status == "INDEXED":
            print(f"✓ {res['pockets']} pockets (pLDDT {res['plddt']})")
        elif status == "NO_POCKETS":
            print(f"○ No pockets (pLDDT {res['plddt']})")
        elif status == "LOW_CONFIDENCE":
            print(f"✗ Low conf (pLDDT {res.get('plddt', 0)})")
        else:
            print(f"✗ {status}")
    
    dt = time() - t0
    print("="*60)
    print(f"[+] Complete in {dt:.1f}s (gemmi pLDDT extraction)")
    print(f"    INDEXED:        {stats['INDEXED']}")
    print(f"    NO_POCKETS:     {stats['NO_POCKETS']}")
    print(f"    LOW_CONFIDENCE: {stats['LOW_CONFIDENCE']}")
    print(f"    FAILED:         {stats['DOWNLOAD_FAILED'] + stats['DETECTION_ERROR']}")

if __name__ == "__main__":
    main()
