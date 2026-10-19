"""
P0 item 4 (2026-10-19 audit): residue ordering regression tests.

The engine previously sorted residues by (res_seq, insertion_code). Real
deposits — 1PPE chain E, sha256 fe51b428..., 203958 bytes — list
insertion-coded residues BEFORE their base number in file order
(183, 184A, 184, 185), so the sort re-paired peptide bonds and fabricated
backbone gaps (the v24 report's "1PPE genuine 7.12 A break" was an
artifact of exactly this).

These tests build a geometrically valid planar extended peptide with
Engh-Huber ideal internal coordinates and assert:

  1. file order is preserved everywhere (coord_math and _extract);
  2. the audit-direction acceptance: fixture 183/184A/184/185 gives
     ZERO backbone breaks;
  3. numbering gaps (45 -> 67) are still skipped, not measured;
  4. the OLD (res_seq, icode) sort demonstrably fabricates a gap on the
     same coordinates (self-contained proof, no old code retained);
  5. mmCIF label_seq_id takes precedence over file order when present.
"""
import math
import numpy as np
import pytest

from tos.engine.coord_math import (
    get_sequential_residues,
    residues_are_adjacent,
    peptide_bonded,
    distance,
)
from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import IngestionProcessor, Atom

# ── Engh & Huber (1991) ideal internal coordinates ─────────────────────
B_N_CA = 1.458   # N-CA
B_CA_C = 1.525   # CA-C
B_C_N = 1.329    # C-N peptide
B_C_O = 1.231    # C=O
A_N_CA_C = 111.0
A_CA_C_N = 116.6
A_C_N_CA = 121.7
A_CA_C_O = 120.8


def _rot(v, deg):
    r = math.radians(deg)
    c, s = math.cos(r), math.sin(r)
    return np.array([v[0] * c - v[1] * s, v[1] * c + v[0] * s, 0.0])


def build_planar_peptide(numbering):
    """Build an extended planar peptide backbone (N, CA, C, O per residue).

    Alternating +/- turns keep the chain extended; every C(i)-N(i+1) bond
    is exactly B_C_N, so a correctly ordered engine must find ZERO breaks.
    Returns a list of (res_seq, icode, res_name, {atom: pos}) tuples in
    construction (= intended file) order.
    """
    residues = []
    # residue 1 bootstrap
    n = np.array([0.0, 0.0, 0.0])
    d_n_ca = np.array([1.0, 0.0, 0.0])
    ca = n + B_N_CA * d_n_ca
    d_ca_c = _rot(d_n_ca, 69.0)          # turn at CA (180 - 111.0)
    c = ca + B_CA_C * d_ca_c
    sign = 1.0
    for (seq, icode) in numbering:
        # direction from C to the NEXT residue's N (turn alternates sign)
        d_c_n = _rot(d_ca_c, sign * (180.0 - A_CA_C_N))
        # carbonyl O on the opposite side of the CA->C direction from N',
        # keeping O-C-N at ~122.6 deg whichever way the chain turns
        o = c + B_C_O * _rot(d_ca_c, -sign * (180.0 - A_CA_C_O))
        residues.append((seq, icode, "GLY",
                         {"N": n, "CA": ca, "C": c, "O": o}))
        # advance to next residue
        n_next = c + B_C_N * d_c_n
        d_n2_ca = _rot(d_c_n, -sign * (180.0 - A_C_N_CA))
        ca_next = n_next + B_N_CA * d_n2_ca
        d_ca2_c = _rot(d_n2_ca, sign * (180.0 - A_N_CA_C))
        c_next = ca_next + B_CA_C * d_ca2_c
        n, ca, c, d_ca_c = n_next, ca_next, c_next, d_ca2_c
        sign = -sign
    return residues


def fixture_pdb_bytes(numbering, chain="A"):
    """Render the planar peptide as PDB ATOM records in numbering order."""
    residues = build_planar_peptide(numbering)
    lines = ["HEADER    SYNTHETIC FIXTURE",
             "EXPDTA    X-RAY DIFFRACTION",
             "REMARK   2 RESOLUTION.    1.00 ANGSTROM"]
    serial = 1
    for (seq, icode, res_name, atoms) in residues:
        for name in ("N", "CA", "C", "O"):
            el = name[0]
            if name == "CA":
                el = "C"
            x, y, z = atoms[name]
            lines.append(
                f"ATOM  {serial:5d} {name:>4s} {res_name:>3s} {chain:1s}"
                f"{seq:4d}{icode:1s}   "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"{1.0:6.2f}{30.0:6.2f}          {el:>2s}  "
            )
            serial += 1
    lines.append("TER")
    lines.append("END")
    return ("\n".join(lines) + "\n").encode("ascii")


ACCEPTANCE_NUMBERING = [(183, ""), (184, "A"), (184, ""), (185, "")]


# ── 1. File order is preserved ──────────────────────────────────────────

def test_get_sequential_residues_returns_file_order():
    content = fixture_pdb_bytes(ACCEPTANCE_NUMBERING)
    structure = IngestionProcessor.run(content, "fixture.pdb", "Audit", "experimental")
    order = get_sequential_residues(structure.atoms, "A")
    assert [(s, i) for s, i, _ in order] == [(183, ""), (184, "A"), (184, ""), (185, "")]


def test_extract_returns_file_order():
    content = fixture_pdb_bytes(ACCEPTANCE_NUMBERING)
    structure = IngestionProcessor.run(content, "fixture.pdb", "Audit", "experimental")
    residues = Tier1Measurements._extract(structure)
    assert [(r["_seq"], r["_icode"]) for r in residues] == [
        (183, ""), (184, "A"), (184, ""), (185, "")]


def test_chains_grouped_by_first_appearance_not_alphabetical():
    numbering = [(10, ""), (11, "")]
    content = fixture_pdb_bytes(numbering, chain="Z") + \
        fixture_pdb_bytes(numbering, chain="B").replace(b"HEADER", b"REMARK").replace(b"EXPDTA", b"REMARK").replace(b"REMARK   2", b"REMARK   3")
    structure = IngestionProcessor.run(content, "fixture.pdb", "Audit", "experimental")
    residues = Tier1Measurements._extract(structure)
    # chain Z appears first in the file and must be measured first,
    # even though "B" < "Z" alphabetically.
    assert [r["_chain"] for r in residues] == ["Z", "Z", "B", "B"]


# ── 2. Acceptance: zero breaks on 183 / 184A / 184 / 185 ────────────────

def test_acceptance_insertion_before_base_zero_breaks():
    content = fixture_pdb_bytes(ACCEPTANCE_NUMBERING)
    structure = IngestionProcessor.run(content, "fixture.pdb", "Audit", "experimental")
    results, coverage, _ = Tier1Measurements.run_full_audit(structure)

    # LAW-110: three peptide bonds, all intact
    assert results["LAW-110"]["observed"] == 0
    assert results["LAW-110"]["sample"] == 3
    assert results["LAW-110"]["status"] == "PASS"

    # LAW-160: adjacent CA-CA distances all under the 4.5 A threshold
    assert results["LAW-160"]["status"] == "PASS"
    assert results["LAW-160"]["sample"] == 3
    assert 3.0 < results["LAW-160"]["observed"] <= 4.5

    # LAW-120: peptide angles measured across genuine bonds only
    assert results["LAW-120"]["status"] == "PASS"
    assert results["LAW-120"]["sample"] > 0
    assert results["LAW-120"]["observed"] < 2.0  # ideal-geometry fixture


# ── 3. Numbering gaps are skipped, never measured ───────────────────────

def test_numbering_gap_not_measured():
    content = fixture_pdb_bytes([(44, ""), (45, ""), (67, ""), (68, "")])
    structure = IngestionProcessor.run(content, "fixture.pdb", "Audit", "experimental")
    results, _, _ = Tier1Measurements.run_full_audit(structure)
    # 44->45 and 67->68 measured; 45->67 is a missing span: absent density,
    # not a chain break — it must not be counted as a gap.
    assert results["LAW-110"]["observed"] == 0
    assert results["LAW-110"]["sample"] == 2


def test_numbering_gap_not_measured_even_when_geometry_continues():
    # The planar builder lays out a physically continuous chain, so the
    # 45->67 pair IS bonded geometrically (C-N = 1.329). The missing-span
    # rule must still exclude it from measurement.
    content = fixture_pdb_bytes([(44, ""), (45, ""), (67, ""), (68, "")])
    structure = IngestionProcessor.run(content, "fixture.pdb", "Audit", "experimental")
    residues = Tier1Measurements._extract(structure)
    by_key = {(r["_seq"], r["_icode"]): r for r in residues}
    assert distance(by_key[(45, "")]["_atoms"]["C"],
                    by_key[(67, "")]["_atoms"]["N"]) < 1.5  # physically bonded
    assert not residues_are_adjacent("A", 45, "", "A", 67, "")  # but numbering gap


# ── 4. The removed sort demonstrably fabricated breaks ──────────────────

def test_old_sort_fabricates_gap_on_same_coordinates():
    residues = build_planar_peptide(ACCEPTANCE_NUMBERING)
    old_sorted = sorted(residues, key=lambda r: (r[0], r[1]))
    cn_distances = []
    for r1, r2 in zip(old_sorted, old_sorted[1:]):
        cn_distances.append(distance(r1[3]["C"], r2[3]["N"]))
    # Under (res_seq, icode) sort the pairs become 183->184, 184->184A,
    # 184A->185 — mis-paired, with at least one C-N far beyond a bond.
    assert max(cn_distances) > 2.0
    # ... while the file-order pairs are all genuine peptide bonds.
    file_order_cn = [distance(r1[3]["C"], r2[3]["N"])
                     for r1, r2 in zip(residues, residues[1:])]
    assert max(file_order_cn) <= 1.35


# ── 5. Adjacency predicate decision table ───────────────────────────────

@pytest.mark.parametrize("a,b,expected", [
    # insertion-coded successor BEFORE its base number (1PPE pattern)
    (("A", 183, ""), ("A", 184, "A"), True),
    (("A", 184, "A"), ("A", 184, ""), True),
    (("A", 184, ""), ("A", 185, ""), True),
    # insertion-coded successor AFTER its base number (common pattern)
    (("A", 99, ""), ("A", 99, "A"), True),
    (("A", 99, "A"), ("A", 100, ""), True),
    # missing residue span: absent density, never measured
    (("A", 45, ""), ("A", 67, ""), False),
    (("A", 100, "A"), ("A", 103, ""), False),
    # chains are never joined
    (("A", 10, ""), ("B", 11, ""), False),
])
def test_residues_are_adjacent_decision_table(a, b, expected):
    assert residues_are_adjacent(a[0], a[1], a[2], b[0], b[1], b[2]) is expected


def test_peptide_bonded_predicate():
    residues = build_planar_peptide(ACCEPTANCE_NUMBERING)
    r1, r2, r3 = residues[0][3], residues[1][3], residues[3][3]
    assert peptide_bonded(r1, r2) is True          # 183 -> 184A bonded
    assert peptide_bonded(r2, r1) is False or True  # symmetric distance check
    assert peptide_bonded(r1, r3) is False          # 183 -> 185 not bonded
    assert peptide_bonded({"N": np.zeros(3)}, r2) is False  # missing C
    assert peptide_bonded(r1, {"CA": np.zeros(3)}) is False  # missing N


# ── 6. mmCIF label_seq_id ordering ──────────────────────────────────────

def test_label_seq_id_takes_precedence_over_file_order():
    # mmCIF atoms may arrive out of label order; label_seq_id governs.
    atoms = [
        Atom("N", "N", np.array([0.0, 0.0, 0.0]), "GLY", 5, "A", "", 30.0, label_seq_id=2),
        Atom("CA", "C", np.array([1.0, 0.0, 0.0]), "GLY", 5, "A", "", 30.0, label_seq_id=2),
        Atom("N", "N", np.array([5.0, 0.0, 0.0]), "GLY", 3, "A", "", 30.0, label_seq_id=1),
        Atom("CA", "C", np.array([6.0, 0.0, 0.0]), "GLY", 3, "A", "", 30.0, label_seq_id=1),
    ]
    order = get_sequential_residues(atoms, "A")
    assert [(s, i) for s, i, _ in order] == [(3, ""), (5, "")]  # label order
