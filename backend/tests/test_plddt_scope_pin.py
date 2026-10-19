"""
P0 item 10: pin the CURRENT pLDDT scope of LAW-130/155/200.

Verified: for predicted structures these three laws measure ALL atoms
while every other law restricts to the pLDDT >= 70 core. Whether that
is correct is an OPEN OWNER DECISION (evidence/item-10/
DECISION_REQUIRED.md) — a verdict-changing scientific-semantics choice.

These tests pin the current behavior so it cannot drift silently in
either direction without a deliberate, reviewed change. If the owner
selects core-filtering, these tests are the ones to update IN THE SAME
PR as the law change, with the verdict-change table attached.
"""
import numpy as np

from tos.ingestion.processor import IngestionProcessor
from tos.engine.tier1_measurements import Tier1Measurements


def _line(serial, seq, name, x, y, z, b):
    return (f"ATOM  {serial:5d} {name:>4s} GLY A{seq:4d}"
            f"    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"{1.0:6.2f}{b:6.2f}          C  ")


def _predicted_structure():
    """Compact high-confidence core + distant low-confidence tail."""
    lines = ["HEADER    SCOPE PIN FIXTURE", "REMARK   1 ALPHAFOLD-LIKE"]
    serial = 1
    rng = np.random.default_rng(7)
    for seq in range(1, 11):  # core: pLDDT 90, compact
        for name in ("N", "CA", "C", "O"):
            p = rng.normal(0.0, 3.0, 3)
            lines.append(_line(serial, seq, name, *p, 90.0))
            serial += 1
    for i, seq in enumerate(range(11, 17)):  # tail: pLDDT 40, extended
        for name in ("N", "CA", "C", "O"):
            lines.append(_line(serial, seq, name, 80.0 + i * 15.0, 0.0, 0.0, 40.0))
            serial += 1
    lines.append("END")
    return ("\n".join(lines) + "\n").encode("ascii")


def test_law_155_200_sample_counts_all_atoms_for_predicted():
    structure = IngestionProcessor.run(
        _predicted_structure(), "pin.pdb", "Audit", "predicted")
    results, _, _ = Tier1Measurements.run_full_audit(structure)
    # All 64 atoms (16 residues x 4) are scored even though only the
    # 10-residue core has pLDDT >= 70. This is the pinned CURRENT scope.
    assert results["LAW-155"]["sample"] == len(structure.atoms)
    assert results["LAW-200"]["sample"] == len(structure.atoms)


def test_law_105_coverage_still_uses_core():
    # The core filter itself is unchanged: coverage reflects core/source.
    structure = IngestionProcessor.run(
        _predicted_structure(), "pin.pdb", "Audit", "predicted")
    results, _, _ = Tier1Measurements.run_full_audit(structure)
    assert results["LAW-105"]["observed"] == round(10 / 16 * 100, 2)


def test_experimental_structures_unaffected_by_scope_question():
    # For experimental structures core == all residues, so the scope
    # distinction is moot there.
    content = _predicted_structure().replace(b"ALPHAFOLD-LIKE", b"X-RAY")
    content = content.replace(b"REMARK   1 X-RAY", b"EXPDTA    X-RAY DIFFRACTION")
    structure = IngestionProcessor.run(content, "pin.pdb", "Audit", "experimental")
    results, _, _ = Tier1Measurements.run_full_audit(structure)
    assert results["LAW-105"]["observed"] == 100.0
    assert results["LAW-155"]["sample"] == len(structure.atoms)
