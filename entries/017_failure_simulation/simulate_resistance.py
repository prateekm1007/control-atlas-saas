#!/usr/bin/env python3
"""
Entry 017 — Failure Simulation (Resistance Prediction)
Predicts pocket collapse risk via steric perturbation.
"""

import os
import csv
from prody import fetchPDB, parsePDB

print("=== ENTRY 017 START ===")

CACHE = os.path.expanduser("~/control-atlas/library/pdb_cache/")
OUTPUT = os.path.expanduser("~/control-atlas/entries/017_failure_simulation/resistance_predictions.csv")

os.makedirs(CACHE, exist_ok=True)

# Residue volumes (Å³)
VOLUME = {
    'GLY':60,'ALA':88,'SER':89,'CYS':108,'ASP':111,'PRO':112,'ASN':114,'THR':116,
    'GLU':138,'VAL':140,'GLN':143,'HIS':153,'MET':162,'ILE':166,'LEU':166,
    'LYS':168,'ARG':173,'PHE':189,'TYR':193,'TRP':227
}

# Known resistance sites
KNOWN = {68:'R68S', 96:'Y96D', 99:'Q99L'}

# Load reference
pdb = fetchPDB("6OIM", folder=CACHE)
structure = parsePDB(pdb)

print("Loaded", structure.numAtoms(), "atoms")
print("")
print("POCKET RESIDUE RISK ANALYSIS")
print("-" * 60)

results = []

for resnum in range(60, 100):
    sel = structure.select(f"chain A and resnum {resnum} and name CA")
    if sel is None:
        continue

    resname = sel.getResnames()[0]
    base_vol = VOLUME.get(resname, 100)
    trp_delta = VOLUME['TRP'] - base_vol

    if trp_delta > 100:
        risk = "HIGH"
    elif trp_delta > 50:
        risk = "MED"
    else:
        risk = "LOW"

    known = KNOWN.get(resnum, "")
    flag = " <-- KNOWN RESISTANCE" if known else ""

    print(f"Res {resnum:3d} {resname:3s} | Vol {base_vol:3d} | ΔTrp {trp_delta:3d} | Risk {risk}{flag}")

    results.append({
        "resnum": resnum,
        "residue": resname,
        "volume": base_vol,
        "trp_delta": trp_delta,
        "risk": risk,
        "known": known
    })

# Write CSV
with open(OUTPUT, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=["resnum","residue","volume","trp_delta","risk","known"]
    )
    writer.writeheader()
    writer.writerows(results)

print("")
print("=== SUMMARY ===")
high = [r for r in results if r["risk"] == "HIGH"]
print("High-risk positions:", len(high))
for r in high:
    label = r["known"] if r["known"] else "Novel"
    print(f"  → Res {r['resnum']} {r['residue']} ({label})")

print("")
print("Saved:", OUTPUT)
print("=== ENTRY 017 DONE ===")
