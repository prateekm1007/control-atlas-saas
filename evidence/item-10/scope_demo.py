#!/usr/bin/env python3
"""P0 item 10 evidence: LAW-130/155/200 atom scope, quantified.

Verified from the code (commit d5fa87a, backend/tos/engine/
tier1_measurements.py):

  - LAW-130 (Clashscore):     all_coords = [a.pos for a in structure.atoms]
  - LAW-155 (Voxel Occupancy): same all_coords
  - LAW-200 (Packing Quality): same all_coords
  - LAW-100/110/120/125/135/145/150/160: iterate `core`
    (all residues for experimental; pLDDT >= 70 residues for predicted)

So for PREDICTED structures the three global-geometry heuristics score
EVERY atom regardless of model confidence, while every other law
restricts to the confident core. This script demonstrates that the
choice changes verdicts, with a synthetic AlphaFold-like structure:

  - a compact 20-residue helical-ish core, pLDDT 90
  - a 10-residue low-confidence tail (pLDDT 40) stretched far away
    (unstructured, as low-confidence AF regions genuinely are)

It then computes LAW-155 and LAW-200 both ways (all atoms vs core only)
and prints the verdict each way.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

import numpy as np

from tos.ingestion.processor import IngestionProcessor
from tos.engine.tier1_measurements import Tier1Measurements


def build_structure():
    """20-residue compact core (pLDDT 90) + 10-residue distant tail (pLDDT 40)."""
    lines = ["HEADER    SYNTHETIC AF-SCOPE DEMO", "REMARK   1 ALPHAFOLD-LIKE"]
    serial = 1
    # core: residues 1..20, compact blob around origin
    rng = np.random.default_rng(42)
    for seq in range(1, 21):
        for name in ("N", "CA", "C", "O"):
            pos = rng.normal(0.0, 4.0, 3)  # compact
            lines.append(
                f"ATOM  {serial:5d} {name:>4s} GLY A{seq:4d}"
                f"    {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}"
                f"{1.0:6.2f}{90.0:6.2f}          C  ")
            serial += 1
    # tail: residues 21..30, pLDDT 40, stretched along +x, far from core
    for i, seq in enumerate(range(21, 31)):
        for name in ("N", "CA", "C", "O"):
            x = 60.0 + i * 12.0  # extended, low-density
            lines.append(
                f"ATOM  {serial:5d} {name:>4s} GLY A{seq:4d}"
                f"    {x:8.3f}{0.0:8.3f}{0.0:8.3f}"
                f"{1.0:6.2f}{40.0:6.2f}          C  ")
            serial += 1
    lines.append("END")
    return ("\n".join(lines) + "\n").encode("ascii")


def main():
    content = build_structure()
    structure = IngestionProcessor.run(content, "demo.pdb", "Audit", "predicted")
    residues = Tier1Measurements._extract(structure)
    core = [r for r in residues if r["_conf"] >= 70]
    print(f"parsed residues: {len(residues)}  core (pLDDT>=70): {len(core)}")

    results, coverage, _ = Tier1Measurements.run_full_audit(structure)
    print(f"\nLAW-105 coverage (core/source): {results['LAW-105']['observed']}%")

    all_atoms = np.array([a.pos for a in structure.atoms])
    core_atoms = np.array([a.pos for a in structure.atoms if a.b_iso >= 70])

    def voxel_metric(coords):
        voxel_size = 2.0
        bbox_min = np.min(coords, axis=0)
        idx = set()
        for pos in coords:
            idx.add(tuple(int(np.floor((p - b) / voxel_size))
                          for p, b in zip(pos, bbox_min)))
        return len(idx)

    def packing_metric(coords):
        bbox = np.maximum(np.max(coords, axis=0) - np.min(coords, axis=0), 1.0)
        return float(bbox[0] * bbox[1] * bbox[2]) / len(coords)

    print("\nLAW-155 Voxel Occupancy (threshold >= 2.0 voxels/residue):")
    v_all = voxel_metric(all_atoms) / len(residues)
    v_core = voxel_metric(core_atoms) / len(core)
    print(f"  all atoms : {v_all:7.2f}  -> {'PASS' if v_all >= 2.0 else 'VETO'}   (current behavior)")
    print(f"  core only : {v_core:7.2f}  -> {'PASS' if v_core >= 2.0 else 'VETO'}   (if core-filtered)")

    print("\nLAW-200 Packing Quality (threshold <= 300 A^3/atom):")
    p_all = packing_metric(all_atoms)
    p_core = packing_metric(core_atoms)
    print(f"  all atoms : {p_all:9.2f}  -> {'PASS' if p_all <= 300 else 'VETO'}   (current behavior)")
    print(f"  core only : {p_core:9.2f}  -> {'PASS' if p_core <= 300 else 'VETO'}   (if core-filtered)")

    print("\nVerdict-relevant divergence:")
    print(f"  LAW-155: {'DIVERGES' if (v_all >= 2.0) != (v_core >= 2.0) else 'same'}"
          f"   LAW-200: {'DIVERGES' if (p_all <= 300) != (p_core <= 300) else 'same'}")
    print("\nOwner decision required: keep all-atoms scope (current) or apply the")
    print("core filter consistently. See evidence/item-10/DECISION_REQUIRED.md.")


if __name__ == "__main__":
    main()
