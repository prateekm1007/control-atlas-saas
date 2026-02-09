#!/usr/bin/env python3
"""
Entry 063 — Antibody–Antigen Interaction Builder
Prepares RFAA jobs for antibodies that passed Entry 062.
"""
import json
import argparse
import os
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sieve", default="../062_constraint_filter/sieve_results_062.json")
    parser.add_argument("--candidates", default="../060_generative_discovery/candidates.json")
    parser.add_argument("--antigen", required=True, help="FASTA file for target antigen")
    parser.add_argument("--out", default="interaction_jobs_063.json")
    args = parser.parse_args()

    if not os.path.exists(args.sieve):
        print(f"❌ Sieve file not found: {args.sieve}")
        return

    with open(args.sieve) as f:
        sieve = json.load(f)

    with open(args.candidates) as f:
        candidates = {c["id"]: c for c in json.load(f)}

    jobs = []
    print(f"🧬 [Entry 063] Staging Interaction Jobs...")
    
    for entry in sieve:
        if "CLEARED" not in entry["decision"]:
            continue

        cid = entry["id"]
        c = candidates[cid]

        # Structure the Job
        jobs.append({
            "job_id": f"{cid}_vs_KRAS",
            "antibody": {
                "id": cid,
                "sequence": c["heavy_chain"] + c["light_chain"], # scFv
                "cdr3": c["cdr3"]
            },
            "antigen": {
                "name": "KRAS_G12D",
                "fasta_path": os.path.abspath(args.antigen)
            },
            "status": "READY_FOR_RFAA"
        })
        print(f"  [+] Staged: {cid} vs KRAS_G12D")

    Path(args.out).write_text(json.dumps(jobs, indent=2))
    print(f"\n✅ Manifest written to {args.out} ({len(jobs)} jobs)")

if __name__ == "__main__":
    main()
