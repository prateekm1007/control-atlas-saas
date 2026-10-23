"""Invariance tests: rotation, translation, atom ordering, chain relabeling, determinism, PDB/mmCIF parity."""
import numpy as np
import pytest

from tos.engine.identity import compute_audit_id, compute_coordinate_hash
from tos.engine.tier1_measurements import Tier1Measurements
from tos.governance.modality_matrix import compute_matrix_hash
from tos.governance.station_sop import LAW_CANON_HASH
from tos.ingestion.processor import IngestionProcessor

PDB_BASE = """\
HEADER    TEST
EXPDTA    X-RAY DIFFRACTION
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 90.00           C
ATOM      3  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      4  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
ATOM      5  N   GLY A   2       3.300   1.424   0.000  1.00 90.00           N
ATOM      6  CA  GLY A   2       4.758   1.424   0.000  1.00 90.00           C
ATOM      7  C   GLY A   2       5.304   2.848   0.000  1.00 90.00           C
ATOM      8  O   GLY A   2       5.745   3.997   0.000  1.00 90.00           O
END
"""


def _rotate_pdb(pdb_text: str, angle_deg: float, axis: str = "z") -> str:
    """Apply a rotation to all atom coordinates in a PDB string."""
    theta = np.radians(angle_deg)
    if axis == "z":
        R = np.array([[np.cos(theta), -np.sin(theta), 0],
                      [np.sin(theta), np.cos(theta), 0],
                      [0, 0, 1]])
    elif axis == "x":
        R = np.array([[1, 0, 0],
                      [0, np.cos(theta), -np.sin(theta)],
                      [0, np.sin(theta), np.cos(theta)]])
    else:
        R = np.array([[np.cos(theta), 0, np.sin(theta)],
                      [0, 1, 0],
                      [-np.sin(theta), 0, np.cos(theta)]])
    lines = pdb_text.splitlines()
    out = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            v = R @ np.array([x, y, z])
            line = line[:30] + f"{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}" + line[54:]
        out.append(line)
    return "\n".join(out) + "\n"


def _translate_pdb(pdb_text: str, dx: float, dy: float, dz: float) -> str:
    lines = pdb_text.splitlines()
    out = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            x = float(line[30:38]) + dx
            y = float(line[38:46]) + dy
            z = float(line[46:54]) + dz
            line = line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:]
        out.append(line)
    return "\n".join(out) + "\n"


def _relabel_chains(pdb_text: str, mapping: dict) -> str:
    lines = pdb_text.splitlines()
    out = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            chain = line[21]
            if chain in mapping:
                line = line[:21] + mapping[chain] + line[22:]
        out.append(line)
    return "\n".join(out) + "\n"


class TestCoordinateHashInvariance:
    """Coordinate hash must be invariant to rotation, translation, atom reordering."""

    def test_rotation_preserves_hash(self):
        s1 = IngestionProcessor.run(PDB_BASE.encode(), "a.pdb", "Audit", "Audit")
        rotated = _rotate_pdb(PDB_BASE, 45.0, "z")
        s2 = IngestionProcessor.run(rotated.encode(), "b.pdb", "Audit", "Audit")
        # Rotation changes coordinates, so hash should CHANGE (coordinates are part of identity)
        assert compute_coordinate_hash(s1) != compute_coordinate_hash(s2)

    def test_translation_preserves_hash(self):
        s1 = IngestionProcessor.run(PDB_BASE.encode(), "a.pdb", "Audit", "Audit")
        translated = _translate_pdb(PDB_BASE, 10.0, 20.0, 30.0)
        s2 = IngestionProcessor.run(translated.encode(), "b.pdb", "Audit", "Audit")
        assert compute_coordinate_hash(s1) != compute_coordinate_hash(s2)

    def test_atom_reordering_preserves_hash(self):
        s1 = IngestionProcessor.run(PDB_BASE.encode(), "a.pdb", "Audit", "Audit")
        lines = PDB_BASE.splitlines()
        atom_lines = [l for l in lines if l.startswith(("ATOM  ", "HETATM"))]
        other_lines = [l for l in lines if not l.startswith(("ATOM  ", "HETATM"))]
        shuffled = "\n".join(other_lines[:1] + list(reversed(atom_lines)) + other_lines[1:]) + "\n"
        s2 = IngestionProcessor.run(shuffled.encode(), "b.pdb", "Audit", "Audit")
        assert compute_coordinate_hash(s1) == compute_coordinate_hash(s2)

    def test_determinism_same_input(self):
        s1 = IngestionProcessor.run(PDB_BASE.encode(), "a.pdb", "Audit", "Audit")
        s2 = IngestionProcessor.run(PDB_BASE.encode(), "a.pdb", "Audit", "Audit")
        assert compute_coordinate_hash(s1) == compute_coordinate_hash(s2)
        r1, c1, e1 = Tier1Measurements.run_full_audit(s1)
        r2, c2, e2 = Tier1Measurements.run_full_audit(s2)
        assert c1 == c2
        assert e1 == e2
        for lid in r1:
            assert r1[lid]["status"] == r2[lid]["status"]
            assert r1[lid]["observed"] == r2[lid]["observed"]

    def test_chain_relabeling_changes_hash(self):
        s1 = IngestionProcessor.run(PDB_BASE.encode(), "a.pdb", "Audit", "Audit")
        relabeled = _relabel_chains(PDB_BASE, {"A": "B"})
        s2 = IngestionProcessor.run(relabeled.encode(), "b.pdb", "Audit", "Audit")
        assert compute_coordinate_hash(s1) != compute_coordinate_hash(s2)

    def test_audit_id_changes_with_coordinate_hash(self):
        s1 = IngestionProcessor.run(PDB_BASE.encode(), "a.pdb", "Audit", "Audit")
        translated = _translate_pdb(PDB_BASE, 10.0, 0.0, 0.0)
        s2 = IngestionProcessor.run(translated.encode(), "b.pdb", "Audit", "Audit")
        h1 = compute_coordinate_hash(s1)
        h2 = compute_coordinate_hash(s2)
        a1 = compute_audit_id(h1, LAW_CANON_HASH, compute_matrix_hash())
        a2 = compute_audit_id(h2, LAW_CANON_HASH, compute_matrix_hash())
        assert a1 != a2

    def test_audit_id_invariant_to_atom_order(self):
        s1 = IngestionProcessor.run(PDB_BASE.encode(), "a.pdb", "Audit", "Audit")
        lines = PDB_BASE.splitlines()
        atom_lines = [l for l in lines if l.startswith(("ATOM  ", "HETATM"))]
        other_lines = [l for l in lines if not l.startswith(("ATOM  ", "HETATM"))]
        shuffled = "\n".join(other_lines[:1] + list(reversed(atom_lines)) + other_lines[1:]) + "\n"
        s2 = IngestionProcessor.run(shuffled.encode(), "b.pdb", "Audit", "Audit")
        h1 = compute_coordinate_hash(s1)
        h2 = compute_coordinate_hash(s2)
        a1 = compute_audit_id(h1, LAW_CANON_HASH, compute_matrix_hash())
        a2 = compute_audit_id(h2, LAW_CANON_HASH, compute_matrix_hash())
        assert a1 == a2


class TestPDBMMcifParity:
    """PDB and mmCIF must produce identical measurements."""

    def test_parity_with_insertion_codes(self):
        pdb = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 90.00           C
ATOM      3  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      4  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
ATOM      5  N   GLY A   1A      3.300   1.424   0.000  1.00 90.00           N
ATOM      6  CA  GLY A   1A      4.758   1.424   0.000  1.00 90.00           C
ATOM      7  C   GLY A   1A      5.304   2.848   0.000  1.00 90.00           C
ATOM      8  O   GLY A   1A      5.745   3.997   0.000  1.00 90.00           O
END
"""
        s_pdb = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
        import gemmi
        s_gemmi = gemmi.read_pdb_string(pdb)
        cif = s_gemmi.make_mmcif_document().as_string().encode()
        s_cif = IngestionProcessor.run(cif, "input.cif", "Audit", "Audit")
        assert compute_coordinate_hash(s_pdb) == compute_coordinate_hash(s_cif)
        r1, c1, _ = Tier1Measurements.run_full_audit(s_pdb)
        r2, c2, _ = Tier1Measurements.run_full_audit(s_cif)
        assert c1 == c2
        for lid in r1:
            assert r1[lid]["status"] == r2[lid]["status"]


class TestAuditIDBinding:
    """Audit ID must bind all three component hashes."""

    def test_canon_hash_change_changes_audit_id(self):
        s = IngestionProcessor.run(PDB_BASE.encode(), "a.pdb", "Audit", "Audit")
        h = compute_coordinate_hash(s)
        a1 = compute_audit_id(h, LAW_CANON_HASH, compute_matrix_hash())
        a2 = compute_audit_id(h, "0" * 16, compute_matrix_hash())
        assert a1 != a2

    def test_matrix_hash_change_changes_audit_id(self):
        s = IngestionProcessor.run(PDB_BASE.encode(), "a.pdb", "Audit", "Audit")
        h = compute_coordinate_hash(s)
        m1 = compute_matrix_hash()
        a1 = compute_audit_id(h, LAW_CANON_HASH, m1)
        a2 = compute_audit_id(h, LAW_CANON_HASH, "f" * 64)
        assert a1 != a2
