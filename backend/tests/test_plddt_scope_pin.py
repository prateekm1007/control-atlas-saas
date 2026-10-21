"""
OD-03 / OD-04: pin the pLDDT >= 70 core scope and the atom-coverage gate.

The owner clarified that predicted structures measure the pLDDT >= 70 core
for LAW-110/130/155/200, while experimental structures keep all atoms.
The global confidence gate is atom-based: fewer than 70 percent of atoms
in high-confidence residues is INDETERMINATE, never PASS. LAW-105 keeps
its residue-coverage meaning and reports the atom basis separately.

These tests are the deliberate scope-change pin. They must fail on the
old all-atom implementation and pass only with the reviewed OD-03/04
patch and its verdict-change evidence.
"""
from pathlib import Path

import numpy as np

from tos.engine.adjudicator import AdjudicationInput, adjudicate_laws
from tos.engine.tier1_measurements import Tier1Measurements
from tos.governance.station_sop import LAW_CANON
from tos.ingestion.processor import IngestionProcessor


def _line(serial, seq, name, x, y, z, b):
    return (f"ATOM  {serial:5d} {name:>4s} GLY A{seq:4d}"
            f"    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"{1.0:6.2f}{b:6.2f}          C  ")


def _predicted_structure():
    """Compact high-confidence core + distant low-confidence tail."""
    lines = ["HEADER    SCOPE PIN FIXTURE", "REMARK   1 ALPHAFOLD-LIKE"]
    serial = 1
    rng = np.random.default_rng(7)
    for seq in range(1, 11):  # core: pLDDT 90, compact, 4 atoms each
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


def _uneven_atom_structure():
    """3 core residues with 1 atom each + 1 low-confidence residue with 10.

    Residue coverage is 3/4 = 75 percent (passes LAW-105), but atom
    coverage is 3/13 = 23.08 percent (must block a PASS).
    """
    lines = ["HEADER    UNEVEN ATOM COVERAGE FIXTURE"]
    serial = 1
    for seq in range(1, 4):
        lines.append(_line(serial, seq, "CA", 0.0, 0.0, seq * 3.0, 90.0))
        serial += 1
    for i, name in enumerate(
            ("N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ", "OXT")):
        lines.append(_line(serial, 4, name, 50.0 + i, 0.0, 0.0, 40.0))
        serial += 1
    lines.append("END")
    return ("\n".join(lines) + "\n").encode("ascii")


def test_law_110_130_155_200_use_core_for_predicted():
    structure = IngestionProcessor.run(
        _predicted_structure(), "pin.pdb", "Audit", "predicted")
    results, _, _ = Tier1Measurements.run_full_audit(structure)

    core_atoms = 10 * 4
    assert results["LAW-110"]["sample"] == 9
    assert results["LAW-130"]["sample"] == core_atoms
    assert results["LAW-155"]["sample"] == core_atoms
    assert results["LAW-200"]["sample"] == core_atoms


def test_law_105_keeps_residue_coverage_and_reports_atom_basis():
    structure = IngestionProcessor.run(
        _predicted_structure(), "pin.pdb", "Audit", "predicted")
    results, coverage, _ = Tier1Measurements.run_full_audit(structure)

    # In this fixture every residue has 4 atoms, so both bases are 62.5.
    assert results["LAW-105"]["observed"] == round(10 / 16 * 100, 2)
    assert results["LAW-105"]["residue_coverage"] == 62.5
    assert results["LAW-105"]["atom_coverage"] == 62.5
    assert coverage == 62.5


def test_atom_coverage_below_70_blocks_pass_even_when_residue_coverage_passes():
    structure = IngestionProcessor.run(
        _uneven_atom_structure(), "uneven.pdb", "Audit", "predicted")
    results, coverage, _ = Tier1Measurements.run_full_audit(structure)

    assert round(coverage, 2) == round(3 / 13 * 100, 2)
    assert results["LAW-105"]["observed"] == 75.0
    assert results["LAW-105"]["residue_coverage"] == 75.0
    assert results["LAW-105"]["atom_coverage"] == round(3 / 13 * 100, 2)
    assert results["LAW-105"]["status"] == "PASS"

    all_pass = {
        lid: {"observed": 0, "deviation": "0.0", "sample": 1, "status": "PASS"}
        for lid in LAW_CANON
    }
    adjudicated = adjudicate_laws(AdjudicationInput(
        raw_measurements=all_pass,
        coverage=coverage,
        is_experimental=False,
        method_type="predicted",
        architecture="NONE",
    ))
    assert adjudicated.verdict == "INDETERMINATE"


def test_no_legacy_plddt_50_boundary_in_reporting():
    artifacts = Path(__file__).resolve().parents[1] / "tos" / "forensic_artifacts"
    for name in ("pdf_generator.py", "confidence_viz.py",
                 "structure_viz.py", "residue_heatmap.py"):
        text = (artifacts / name).read_text(encoding="utf-8")
        assert "threshold=50" not in text, f"{name} still has a 50 threshold"
        assert ">= 50" not in text, f"{name} still has a 50 boundary"
        assert "50-70" not in text, f"{name} still labels a 50-70 band"
        assert "Medium" not in text, f"{name} still has a 50-70 medium band"
        assert "[0, 50" not in text, f"{name} still has a 50 band edge"
        assert "for th in [50" not in text, f"{name} still draws a 50 line"


def test_experimental_structures_keep_all_atoms():
    content = _predicted_structure().replace(b"ALPHAFOLD-LIKE", b"X-RAY")
    content = content.replace(b"REMARK   1 X-RAY", b"EXPDTA    X-RAY DIFFRACTION")
    structure = IngestionProcessor.run(content, "pin.pdb", "Audit", "experimental")
    results, _, _ = Tier1Measurements.run_full_audit(structure)

    assert results["LAW-105"]["observed"] == 100.0
    assert results["LAW-110"]["sample"] == 15
    assert results["LAW-130"]["sample"] == len(structure.atoms)
    assert results["LAW-155"]["sample"] == len(structure.atoms)
    assert results["LAW-200"]["sample"] == len(structure.atoms)
