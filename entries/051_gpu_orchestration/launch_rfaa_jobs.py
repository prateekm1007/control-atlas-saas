#!/usr/bin/env python3
"""
Entry 051 — GPU Job Launcher
Identifies survivors needing RFAA validation and dispatches batch jobs.
"""

import sys
import json
import argparse
from pathlib import Path

BASE = Path(__file__).resolve().parents[2] / "entries"
sys.path.insert(0, str(BASE / "044_knowledge_graph"))
from graph_builder import KnowledgeGraph

# Output directory for RFAA batch lists
BATCH_DIR = Path("rfaa_batches")

def get_pending_hypotheses(kg_path):
    """Find hypotheses that passed pre-filters but lack RFAA evidence."""
    with open(kg_path, 'r') as f:
        data = json.load(f)
    
    # Simple parsing of graph JSON
    # In prod, use NetworkX traversal
    pending = []
    
    # Iterate nodes to find Hypotheses
    for node in data["nodes"]:
        if node.get("type") == "Hypothesis":
            # Check if it has a result (e.g. from Entry 049 screen)
            # We want ones that are "CLEARED" by logic but not "VALIDATED" by RFAA
            # For this demo, we assume any 'CLEARED' hypothesis needs RFAA
            if node.get("result") == "PROVEN TRUE" or node.get("result") == "CLEARED":
                pending.append(node["id"]) # ID format: Target_Compound
                
    return pending

def create_batches(pending_list, batch_size=50):
    BATCH_DIR.mkdir(exist_ok=True)
    
    chunks = [pending_list[i:i + batch_size] for i in range(0, len(pending_list), batch_size)]
    
    print(f"[*] Found {len(pending_list)} pending hypotheses.")
    print(f"[*] Creating {len(chunks)} batch files in {BATCH_DIR}...")
    
    for i, chunk in enumerate(chunks):
        batch_file = BATCH_DIR / f"batch_{i:04d}.json"
        with open(batch_file, "w") as f:
            json.dump(chunk, f, indent=2)
            
    print("[+] Batch generation complete.")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kg", required=True, help="Path to Knowledge Graph JSON")
    parser.add_argument("--size", type=int, default=50, help="Batch size")
    args = parser.parse_args()
    
    pending = get_pending_hypotheses(args.kg)
    if pending:
        create_batches(pending, args.size)
    else:
        print("[*] No pending hypotheses found. (Run Entry 049 first?)")

if __name__ == "__main__":
    main()
