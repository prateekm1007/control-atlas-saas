#!/usr/bin/env python3
"""
Entry 049 — Atlas Batch Screening (Production Hardened)
Supports Resume-on-Failure via Checkpointing.
"""

import sys
import json
import glob
import time
import os
from pathlib import Path
from dataclasses import asdict
from multiprocessing import Pool, cpu_count

BASE = Path(__file__).resolve().parents[2] / "entries"
sys.path.insert(0, str(BASE / "040_unified_navigator"))
sys.path.insert(0, str(BASE / "043_proof_engine"))
sys.path.insert(0, str(BASE / "044_knowledge_graph"))
sys.path.insert(0, str(BASE / "037_chemistry_layer"))

from navigator import UnifiedNavigator
from proof_generator import ProofEngine
from graph_builder import KnowledgeGraph

ATLAS_DIR = Path(__file__).resolve().parents[2] / "library/atlas_index"
KG_PATH = Path(__file__).resolve().parents[2] / "entries/044_knowledge_graph/atlas_kg.json"
CHECKPOINT_DIR = Path("checkpoints")

def load_library(smi_path):
    data = []
    with open(smi_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                data.append((parts[1], parts[0]))
    return data

def process_target_batch(args):
    """Worker function for parallel processing."""
    target_files, library = args
    results = []
    nav = UnifiedNavigator()
    prover = ProofEngine()
    
    CHECKPOINT_DIR.mkdir(exist_ok=True)
    
    for target_file in target_files:
        t_name = Path(target_file).stem
        ckpt_file = CHECKPOINT_DIR / f"{t_name}.done"
        
        # RESUME LOGIC: Skip if done
        if ckpt_file.exists():
            continue
            
        try:
            with open(target_file) as f:
                target_data = json.load(f)
            
            target_id = target_data["uniprot_id"]
            pockets = target_data.get("pockets", [])
            if not pockets: 
                # Mark done even if no pockets to avoid re-processing
                ckpt_file.touch()
                continue
            
            for cmp_id, smiles in library:
                for p in pockets:
                    # v2.0 TODO: Insert RFAA Complex Validation Here
                    
                    chem_res = nav.chem.evaluate({"smiles": smiles, "pocket_metrics": p})
                    nav_result = {
                        "compound_id": cmp_id,
                        "status": "BLOCKED" if chem_res.status == "FAIL" else "CLEARED",
                        "trace": {"physics": p, "chemistry": chem_res.metrics}
                    }
                    proof = prover.generate_proof(nav_result)
                    results.append({
                        "target": target_id,
                        "compound": cmp_id,
                        "proof": asdict(proof)
                    })
            
            # Mark target as completed
            ckpt_file.touch()
            
        except: continue
        
    return results

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--lib", required=True)
    args = parser.parse_args()
    
    library = load_library(args.lib)
    targets = glob.glob(str(ATLAS_DIR / "*.json"))
    
    print(f"[*] Screening {len(library)} compounds against {len(targets)} targets...")
    print(f"[*] Checkpoint Dir: {CHECKPOINT_DIR}")
    
    chunk_size = len(targets) // cpu_count() + 1
    chunks = [targets[i:i + chunk_size] for i in range(0, len(targets), chunk_size)]
    
    t0 = time.time()
    
    with Pool(processes=cpu_count()) as pool:
        batch_results = pool.map(process_target_batch, [(chunk, library) for chunk in chunks])
    
    kg = KnowledgeGraph()
    total_proofs = 0
    
    print("[*] Aggregating results into Knowledge Graph...")
    for batch in batch_results:
        for res in batch:
            kg.add_proof(res["target"], res["compound"], res["proof"])
            total_proofs += 1
            
    dt = time.time() - t0
    print(f"\n[+] Completed {total_proofs} interactions in {dt:.1f}s")
    
    kg.export_json(KG_PATH)
    print(f"[+] Knowledge Graph saved to {KG_PATH}")

if __name__ == "__main__":
    main()
