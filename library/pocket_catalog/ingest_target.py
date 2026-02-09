#!/usr/bin/env python3
"""
Step 1 — Pocket Catalog Ingestion (Frame Scaffold Only)
"""

import os
import json
import argparse

CATALOG = os.path.expanduser("~/control-atlas/library/pocket_catalog")

def ingest(target, mutation, domain, uniprot, pdbs, lining):
    name = f"{target}_{mutation}"
    base = os.path.join(CATALOG, name)

    os.makedirs(os.path.join(base, "states"), exist_ok=True)

    pocket_frame = {
        "schema_version": "1.0",
        "identity": {
            "target": target,
            "mutation": mutation,
            "domain": domain,
            "uniprot_id": uniprot
        },
        "pocket": {
            "pocket_id": f"{target}_{domain}",
            "classification": domain,
            "reference_state": "drug_bound",
            "source_pdbs": pdbs.split(",")
        },
        "frame": {
            "origin": [0.0, 0.0, 0.0],
            "axes": [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0]
            ],
            "lining_residues": [int(r) for r in lining.split(",")]
        },
        "quality": {
            "physics_consistency": 0.0,
            "notes": "Scaffold only — awaiting physics computation"
        }
    }

    with open(os.path.join(base, "pocket_frame.json"), "w") as f:
        json.dump(pocket_frame, f, indent=2)

    with open(os.path.join(base, "physics_metrics.json"), "w") as f:
        json.dump({}, f, indent=2)

    with open(os.path.join(base, "lining_residues.csv"), "w") as f:
        f.write("resnum\n")
        for r in pocket_frame["frame"]["lining_residues"]:
            f.write(f"{r}\n")

    print(f"[OK] Pocket scaffold created at {base}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--mutation", required=True)
    ap.add_argument("--domain", required=True)
    ap.add_argument("--uniprot", required=True)
    ap.add_argument("--pdbs", required=True)
    ap.add_argument("--lining", required=True)
    args = ap.parse_args()

    ingest(
        args.target,
        args.mutation,
        args.domain,
        args.uniprot,
        args.pdbs,
        args.lining
    )
