#!/usr/bin/env python3
import sys
import os
import csv
import argparse
from datetime import datetime

sys.path.insert(0, os.path.expanduser("~/control-atlas/entries/020_candidate_validation"))
from validate_candidate import validate

def filter_library(input_file, output_dir):
    passed = []
    rejected = []
    errors = []
    
    print("=== LIBRARY FILTERING ===")
    print("Input:", input_file)
    
    with open(input_file, "r") as f:
        lines = f.readlines()
    
    total = len(lines)
    print("Processing", total, "compounds...")
    
    for i, line in enumerate(lines):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
            
        parts = line.split(",")
        smiles = parts[0].strip()
        mol_id = parts[1].strip() if len(parts) > 1 else "MOL_" + str(i+1)
        
        result = validate(smiles)
        
        if result["status"] == "VALID":
            passed.append({"id": mol_id, "smiles": smiles, "mw": result.get("metrics", {}).get("mw", "N/A")})
            print("  [PASS]", mol_id)
        elif result["status"] == "REJECT":
            notes = result.get("notes", [])
            rejected.append({"id": mol_id, "smiles": smiles, "reasons": "; ".join(notes)})
            print("  [REJECT]", mol_id, ":", notes[0] if notes else "Unknown")
        else:
            errors.append({"id": mol_id, "smiles": smiles, "error": "Parse error"})
            print("  [ERROR]", mol_id)
    
    os.makedirs(output_dir, exist_ok=True)
    
    with open(os.path.join(output_dir, "passed_candidates.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "smiles", "mw"])
        w.writeheader()
        w.writerows(passed)
    
    with open(os.path.join(output_dir, "rejected_candidates.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "smiles", "reasons"])
        w.writeheader()
        w.writerows(rejected)
    
    with open(os.path.join(output_dir, "filtering_stats.txt"), "w") as f:
        f.write("Library Filtering Report" + chr(10))
        f.write("Total: " + str(total) + chr(10))
        f.write("Passed: " + str(len(passed)) + chr(10))
        f.write("Rejected: " + str(len(rejected)) + chr(10))
        f.write("Errors: " + str(len(errors)) + chr(10))
    
    print("")
    print("=== COMPLETE ===")
    print("Passed:", len(passed))
    print("Rejected:", len(rejected))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default=".")
    args = parser.parse_args()
    filter_library(args.input, args.output)
