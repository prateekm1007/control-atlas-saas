#!/usr/bin/env python3
"""
Pocket physics computation — ligand-aware centroid model (FINAL)
"""

import os
import json
import numpy as np
from prody import fetchPDB, parsePDB, confProDy

# Silence ProDy
confProDy(verbosity='none')

CATALOG = os.path.expanduser("~/control-atlas/library/pocket_catalog")
PDB_CACHE = os.path.expanduser("~/control-atlas/library/pdb_cache")

HYDROPHOBIC = {"ALA","VAL","LEU","ILE","MET","PHE","TRP","PRO","TYR"}

def get_centroid(structure, lining_residues):
    # 1. Try ligand centroid (preferred)
    lig = structure.select("hetero and not water")
    if lig is not None and lig.numAtoms() >= 5:
        return lig.getCoords().mean(axis=0), "ligand"

    # 2. Try CA centroid
    res_str = " ".join(str(r) for r in lining_residues)
    ca = structure.select(f"protein and name CA and resnum {res_str}")
    if ca is not None and ca.numAtoms() >= 3:
        return ca.getCoords().mean(axis=0), "ca_residues"

    return None, "none"

def compute_pocket_physics(pdb_id, lining_residues, radius=8.0):
    try:
        pdb_file = fetchPDB(pdb_id, folder=PDB_CACHE)
        if not pdb_file: return {"status": "failed", "error": "PDB download failed"}
        structure = parsePDB(pdb_file)

        centroid, source = get_centroid(structure, lining_residues)
        if centroid is None:
            return {"status": "failed", "error": "No centroid found"}

        protein = structure.select("protein")
        coords = protein.getCoords()

        dists = np.linalg.norm(coords - centroid, axis=1)
        idx = np.where(dists <= radius)[0]

        if len(idx) < 50:
            return {"status": "failed", "error": f"Insufficient pocket atoms ({len(idx)})"}

        pocket_atoms = protein.select(f"index {' '.join(map(str, idx))}")
        pocket_coords = pocket_atoms.getCoords()

        from scipy.spatial import ConvexHull
        hull = ConvexHull(pocket_coords)

        exposure_scores = []
        sample_coords = pocket_coords[::5] if len(pocket_coords) > 500 else pocket_coords
        for pc in sample_coords:
            n = np.sum(np.linalg.norm(coords - pc, axis=1) < 6.0)
            exposure_scores.append(max(0, 100 - n * 2))
        exposure = float(np.mean(exposure_scores))

        resnames = set(pocket_atoms.getResnames())
        hydro_count = sum(1 for r in resnames if r in HYDROPHOBIC)
        hydro_frac = (hydro_count / len(resnames)) * 100 if len(resnames) > 0 else 0

        return {
            "volume_A3": round(hull.volume, 1),
            "exposure": round(exposure, 2),
            "hydrophobic_pct": round(hydro_frac, 1),
            "atom_count": pocket_atoms.numAtoms(),
            "centroid_source": source,
            "status": "computed"
        }
    except Exception as e:
        return {"status": "failed", "error": str(e)}

def process_catalog():
    print("=== COMPUTING POCKET PHYSICS (LIGAND-AWARE) ===")
    success, failed = 0, 0

    for entry in sorted(os.listdir(CATALOG)):
        frame_path = os.path.join(CATALOG, entry, "pocket_frame.json")
        out_path = os.path.join(CATALOG, entry, "physics_metrics.json")

        if not os.path.exists(frame_path): continue

        with open(frame_path) as f: frame = json.load(f)
        pdbs = frame["pocket"]["source_pdbs"]
        residues = frame["frame"]["lining_residues"]

        if not pdbs: continue
        pdb_id = pdbs[0]
        
        print(f"Processing {entry} ({pdb_id})...")
        metrics = compute_pocket_physics(pdb_id, residues)

        if metrics["status"] == "computed":
            with open(out_path, "w") as f: json.dump(metrics, f, indent=2)
            print(f"  [OK] Vol={metrics['volume_A3']}  Hydro={metrics['hydrophobic_pct']}%  ({metrics['centroid_source']})")
            success += 1
        else:
            print(f"  [FAIL] {metrics['error']}")
            failed += 1

    print(f"\nSuccess: {success}, Failed: {failed}")

if __name__ == "__main__":
    process_catalog()
