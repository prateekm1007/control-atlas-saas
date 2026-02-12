#!/usr/bin/env python3
import json
import numpy as np
import os
from Bio import PDB
from Bio.SeqUtils.ProtParam import ProteinAnalysis

CANDIDATES_FILE = "../060_generative_discovery/candidates.json"
FOLDING_FILE = "../061_structure_prediction/folding_results.json"
OUTPUT_FILE = "sieve_results_062.json"

MAX_GRAVY = 0.5
FORBIDDEN_MOTIFS = ["NG", "NS", "DG"]

class AntibodySieve:
    def __init__(self):
        self.parser = PDB.PDBParser(QUIET=True)

    def extract_plddt(self, pdb_path):
        if not os.path.exists(pdb_path):
            return np.array([])
        structure = self.parser.get_structure("scfv", pdb_path)
        scores = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if 'CA' in residue:
                        scores.append(residue['CA'].bfactor)
        scores = np.array(scores)
        if len(scores) > 0 and np.max(scores) <= 1.0:
            scores *= 100.0
        return scores

    def check_chemistry(self, full_seq, cdr3_seq):
        gravy = ProteinAnalysis(full_seq).gravy()
        liabilities = [m for m in FORBIDDEN_MOTIFS if m in cdr3_seq]
        free_cys = "C" in cdr3_seq
        return {
            "pass": (len(liabilities) == 0) and (not free_cys) and (gravy < MAX_GRAVY),
            "gravy": round(gravy, 3),
            "liabilities": liabilities
        }

def main():
    print("🧬 [Entry 062d] Context-Aware Antibody Sieve")

    with open(CANDIDATES_FILE) as f:
        candidates = {c["id"]: c for c in json.load(f)}

    with open(FOLDING_FILE) as f:
        folds = json.load(f)

    sieve = AntibodySieve()
    results = []

    print(f"\n{'ID':<15} | {'Chem':<6} | {'Decision'}")
    print("-" * 45)

    for entry in folds:
        cid = entry["id"]
        if cid not in candidates:
            continue

        cand = candidates[cid]
        full_seq = cand["heavy_chain"] + cand["light_chain"]

        chem = sieve.check_chemistry(full_seq, cand["cdr3"])

        if not chem["pass"]:
            decision = "REJECT_CHEMISTRY"
        else:
            decision = "CLEARED_FOR_RFAA (Contextual)"

        print(f"{cid:<15} | {str(chem['pass']):<6} | {decision}")

        results.append({
            "id": cid,
            "decision": decision,
            "metrics": chem
        })

    with open(OUTPUT_FILE, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n[+] Results written to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
