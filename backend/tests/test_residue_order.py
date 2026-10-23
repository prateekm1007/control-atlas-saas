"""Residue identity / order / adjacency contract tests (P1).

Covers the insertion-code residue-ordering defect: the engine must preserve
authoritative source order (not alphabetical insertion-code sorting) when
pairing peptide geometry, so deposits such as 1PPE/5T4O with insertion-coded
residues cannot fabricate backbone gaps or false vetoes.
"""
import numpy as np
import pytest

from tos.engine.coord_math import get_sequential_residues
from tos.engine.tier1_measurements import INSERTION_ALPHABET, Tier1Measurements
from tos.ingestion.atom_selection import residue_source_order_key
from tos.ingestion.processor import Atom, ConfidenceProfile, IngestionProcessor, Structure


def _atom(serial, name, res, chain, seq, icode, x, y, z, element):
    return (
        f"ATOM  {serial:5d} {name:<4s} {res:>3s} {chain}{seq:4d}{icode:1s}"
        f"   {x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{90.00:6.2f}          {element:>2s}\n"
    )


def _pdb(atoms):
    return "".join(atoms) + "END\n"


def _seq(serial, name, res, seq, icode, x, y, z):
    return _atom(serial, name, res, "A", seq, icode, x, y, z, "C")


def test_insertion_coded_residue_keeps_source_order():
    """1PPE-class: an insertion-coded residue deposited before its base residue
    must keep its source position, not be re-sorted alphabetically."""
    # Source file order: 1A, then 1, then 2. Peptide bonds: 1A->1, 1->2.
    pdb = _pdb([
        _seq(1, "N", "ALA", 1, "A", 0.0, 0.0, 0.0),
        _seq(2, "CA", "ALA", 1, "A", 1.4, 0.0, 0.0),
        _seq(3, "C", "ALA", 1, "A", 2.0, 1.4, 0.0),
        _seq(4, "O", "ALA", 1, "A", 2.5, 2.5, 0.0),
        _seq(5, "N", "ALA", 1, "", 3.0, 0.0, 0.0),
        _seq(6, "CA", "ALA", 1, "", 4.4, 0.0, 0.0),
        _seq(7, "C", "ALA", 1, "", 5.0, 1.4, 0.0),
        _seq(8, "O", "ALA", 1, "", 5.5, 2.5, 0.0),
        _seq(9, "N", "GLY", 2, "", 6.0, 0.0, 0.0),
        _seq(10, "CA", "GLY", 2, "", 7.4, 0.0, 0.0),
        _seq(11, "C", "GLY", 2, "", 8.0, 1.4, 0.0),
        _seq(12, "O", "GLY", 2, "", 8.5, 2.5, 0.0),
    ])
    s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
    # Source order indices increase with file order: 1A < 1 < 2.
    ca_order = {a.res_seq: a.insertion_code for a in s.atoms if a.atom_name == "CA"}
    assert ca_order[1] == "A"
    assert residue_source_order_key(next(a for a in s.atoms if a.atom_name == "CA" and a.insertion_code == "A")) < \
        residue_source_order_key(next(a for a in s.atoms if a.atom_name == "CA" and a.insertion_code == ""))

    residues = Tier1Measurements._extract(s)
    ordered = [(r["_seq"], r["_icode"]) for r in residues]
    # Source order preserved: 1A before 1 before 2.
    assert ordered == [(1, "A"), (1, ""), (2, "")]

    seq_res = get_sequential_residues(s.atoms, "A")
    assert seq_res == [(1, "A", "ALA"), (1, "", "ALA"), (2, "", "GLY")]

    # The exact expected peptide pairs are measured: (1A -> 1) and (1 -> 2).
    assert Tier1Measurements._is_sequential(residues[0], residues[1]) is True
    assert Tier1Measurements._is_sequential(residues[1], residues[2]) is True
    assert residues[1]["_pos"] == residues[0]["_pos"] + 1
    assert residues[2]["_pos"] == residues[1]["_pos"] + 1

    # The peptide geometry must NOT fabricate a gap; both pairs counted.
    results, _, _ = Tier1Measurements.run_full_audit(s)
    assert results["LAW-110"]["observed"] == 0
    assert results["LAW-110"]["status"] == "PASS"
    assert results["LAW-110"]["sample"] == 2


def test_alphaical_insertion_sort_would_fabricate_gap():
    """Demonstrate the ordering-specific verdict change: the OLD (seq, icode)
    sort + legacy number/insertion-code adjacency fabricates a gap; the
    source-order contract measures both correct pairs and reports zero gaps."""
    pdb = _pdb([
        _seq(1, "N", "ALA", 1, "A", 0.0, 0.0, 0.0),
        _seq(2, "CA", "ALA", 1, "A", 1.4, 0.0, 0.0),
        _seq(3, "C", "ALA", 1, "A", 2.0, 1.4, 0.0),
        _seq(4, "O", "ALA", 1, "A", 2.5, 2.5, 0.0),
        _seq(5, "N", "ALA", 1, "", 3.0, 0.0, 0.0),
        _seq(6, "CA", "ALA", 1, "", 4.4, 0.0, 0.0),
        _seq(7, "C", "ALA", 1, "", 5.0, 1.4, 0.0),
        _seq(8, "O", "ALA", 1, "", 5.5, 2.5, 0.0),
        _seq(9, "N", "GLY", 2, "", 6.0, 0.0, 0.0),
        _seq(10, "CA", "GLY", 2, "", 7.4, 0.0, 0.0),
        _seq(11, "C", "GLY", 2, "", 8.0, 1.4, 0.0),
        _seq(12, "O", "GLY", 2, "", 8.5, 2.5, 0.0),
    ])
    s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
    residues = Tier1Measurements._extract(s)

    # Emulate the previous engine: alphabetical (seq, icode, name) ordering plus
    # the legacy number/insertion-code adjacency rule.
    old_order = sorted(residues, key=lambda r: (r["_seq"], r["_icode"], r["_name"]))
    assert [(r["_seq"], r["_icode"]) for r in old_order] == [(1, ""), (1, "A"), (2, "")]

    def legacy_sequential(r1, r2):
        if r1["_chain"] != r2["_chain"]:
            return False
        s1, s2 = int(r1["_seq"]), int(r2["_seq"])
        i1 = r1.get("_icode", "")
        i2 = r2.get("_icode", "")
        rank = {code: idx for idx, code in enumerate(INSERTION_ALPHABET)}
        if s2 == s1:
            rk1, rk2 = rank.get(i1), rank.get(i2)
            return rk1 is not None and rk2 is not None and rk2 == rk1 + 1
        if s2 == s1 + 1:
            return (not i1 and not i2) or (bool(i1) and not i2)
        return False

    def legacy_gap_count(ordered):
        from tos.engine.coord_math import distance
        gaps = 0
        for i in range(len(ordered) - 1):
            r1, r2 = ordered[i], ordered[i + 1]
            if legacy_sequential(r1, r2):
                if "C" in r1["_atoms"] and "N" in r2["_atoms"]:
                    if distance(r1["_atoms"]["C"], r2["_atoms"]["N"]) > 2.0:
                        gaps += 1
        return gaps

    assert legacy_gap_count(old_order) > 0    # old algorithm fabricated a gap
    results, _, _ = Tier1Measurements.run_full_audit(s)
    assert results["LAW-110"]["observed"] == 0
    assert results["LAW-110"]["sample"] == 2


def test_sequence_gap_is_not_fabricated_as_favorable():
    """A genuine missing residue (1 and 3, no 2) must not be turned into a
    favorable zero: the per is not paired, so LAW-110 is honestly
    NOT_MEASURED rather than PASS."""
    pdb = _pdb([
        _seq(1, "N", "ALA", 1, "", 0.0, 0.0, 0.0),
        _seq(2, "CA", "ALA", 1, "", 1.4, 0.0, 0.0),
        _seq(3, "C", "ALA", 1, "", 2.0, 1.4, 0.0),
        _seq(4, "O", "ALA", 1, "", 2.5, 2.5, 0.0),
        _seq(5, "N", "GLY", 3, "", 20.0, 0.0, 0.0),
        _seq(6, "CA", "GLY", 3, "", 21.4, 0.0, 0.0),
        _seq(7, "C", "GLY", 3, "", 22.0, 1.4, 0.0),
        _seq(8, "O", "GLY", 3, "", 22.5, 2.5, 0.0),
    ])
    s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
    residues = Tier1Measurements._extract(s)
    ordered = [(r["_seq"], r["_icode"]) for r in residues]
    assert ordered == [(1, ""), (3, "")]          # source order preserved; missing 2
    results, _, _ = Tier1Measurements.run_full_audit(s)
    assert results["LAW-110"]["status"] == "NOT_MEASURED"


def test_missing_backbone_is_excluded_from_core_but_in_denominator():
    pdb = _pdb([
        _seq(1, "N", "ALA", 1, "", 0.0, 0.0, 0.0),
        _seq(2, "CA", "ALA", 1, "", 1.4, 0.0, 0.0),
        _seq(3, "C", "ALA", 1, "", 2.0, 1.4, 0.0),
        _seq(4, "O", "ALA", 1, "", 2.5, 2.5, 0.0),
        _seq(5, "N", "GLY", 2, "", 3.0, 0.0, 0.0),
        _seq(6, "C", "GLY", 2, "", 5.0, 1.4, 0.0),
        _seq(7, "O", "GLY", 2, "", 5.5, 2.5, 0.0),
    ])
    s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
    results, _, _ = Tier1Measurements.run_full_audit(s)
    law105 = results["LAW-105"]
    # No source sequence is present: coverage fails closed (never an
    # authoritative numeric percentage). The observed-only incomplete-backbone
    # diagnostic is exposed separately.
    assert law105["observed"] is None
    assert law105["status"] == "NOT_MEASURED"
    assert law105["source_backed"] is False
    assert law105["complete_backbone_residues"] == 1
    assert law105["observed_only_diagnostic_pct"] == 50.0


def test_multiple_chains_keep_independent_order():
    pdb = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 90.00           C\n"
        "ATOM      3  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C\n"
        "ATOM      4  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O\n"
        "ATOM      5  N   GLY B   1       10.000   0.000   0.000  1.00 90.00           N\n"
        "ATOM      6  CA  GLY B   1       11.458   0.000   0.000  1.00 90.00           C\n"
        "ATOM      7  C   GLY B   1       12.004   1.424   0.000  1.00 90.00           C\n"
        "ATOM      8  O   GLY B   1       12.445   2.573   0.000  1.00 90.00           O\n"
        "END\n"
    )
    s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
    seq_a = get_sequential_residues(s.atoms, "A")
    seq_b = get_sequential_residues(s.atoms, "B")
    assert seq_a == [(1, "", "ALA")]
    assert seq_b == [(1, "", "GLY")]


def _mk_atom(name, res, chain, seq, icode, x, y, z, label_seq=None, src=None):
    return Atom(
        atom_name=name, element="C", pos=np.array([x, y, z]), res_name=res,
        res_seq=seq, chain_id=chain, insertion_code=icode, b_iso=90.0,
        label_seq_id=label_seq, source_order=src,
    )


def _ala(chain, seq, icode, x, y, z, label_seq=None, src=None):
    return [
        _mk_atom("N", "ALA", chain, seq, icode, x, y, z, label_seq, src),
        _mk_atom("CA", "ALA", chain, seq, icode, x + 1.4, y, z, label_seq, src),
        _mk_atom("C", "ALA", chain, seq, icode, x + 2.0, y + 1.4, z, label_seq, src),
        _mk_atom("O", "ALA", chain, seq, icode, x + 2.5, y + 2.5, z, label_seq, src),
    ]


def _fake_structure(atoms, source_sequences=None):
    from tos.ingestion.processor import IngestionMetadata
    meta = IngestionMetadata()
    if source_sequences:
        meta.source_sequences = source_sequences
    return Structure(
        atoms=atoms,
        audit_id="",
        confidence=ConfidenceProfile([], True, "xray", None),
        metadata=meta,
    )


def test_label_seq_id_basis_governs_ordering_and_adjacency():
    """When valid label_seq_id is present, it is the canonical basis and
    overrides source_order for ordering/adjacency."""
    # Source order says R1 then R2; label_seq_id says R2(1) then R1(2).
    atoms = _ala("A", 1, "", 0.0, 0.0, 0.0, label_seq=2, src=0)
    atoms += _ala("A", 2, "", 3.0, 0.0, 0.0, label_seq=1, src=1)
    s = _fake_structure(atoms)
    residues = Tier1Measurements._extract(s)
    assert [r["_lseq"] for r in residues] == [1, 2]
    assert [r["_basis"] for r in residues] == ["label_seq_id", "label_seq_id"]
    assert [r["_gap_known"] for r in residues] == [True, True]
    # R2 (seq 2) is pos 1, R1 (seq 1) is pos 2.
    assert residues[0]["_pos"] == 1
    assert residues[1]["_pos"] == 2
    assert Tier1Measurements._is_sequential(residues[0], residues[1]) is True
    assert Tier1Measurements._is_sequential(residues[1], residues[0]) is False


def test_absent_sequence_evidence_falls_back_to_source_order():
    atoms = _ala("A", 1, "", 0.0, 0.0, 0.0, label_seq=None, src=0)
    atoms += _ala("A", 2, "", 3.0, 0.0, 0.0, label_seq=None, src=1)
    s = _fake_structure(atoms)
    residues = Tier1Measurements._extract(s)
    assert [r["_basis"] for r in residues] == ["source_order", "source_order"]
    assert [r["_gap_known"] for r in residues] == [False, False]


def test_label_seq_id_gap_yields_explicit_unknown():
    """A genuine missing internal polymer position (label_seq_id 1 then 3)
    must yield an explicit non-measured LAW-110, never a favorable PASS."""
    atoms = _ala("A", 1, "", 0.0, 0.0, 0.0, label_seq=1, src=0)
    atoms += _ala("A", 3, "", 20.0, 0.0, 0.0, label_seq=3, src=1)
    s = _fake_structure(atoms)
    residues = Tier1Measurements._extract(s)
    assert [r["_pos"] for r in residues] == [1, 3]   # pos gap of 2
    results, _, _ = Tier1Measurements.run_full_audit(s)
    assert results["LAW-110"]["status"] == "NOT_MEASURED"
    assert "missing polymer positions" in results["LAW-110"].get("reason", "")


def test_chain_boundary_is_not_sequential():
    atoms = _ala("A", 1, "", 0.0, 0.0, 0.0, src=0)
    atoms += _ala("B", 1, "", 10.0, 0.0, 0.0, src=1)
    s = _fake_structure(atoms)
    residues = Tier1Measurements._extract(s)
    assert Tier1Measurements._is_sequential(residues[0], residues[1]) is False


def test_insertion_coded_pair_is_sequential_in_source_order():
    """An insertion-coded residue placed before its base-numbered residue in
    source order (184A then 184) must still be recognized as adjacent — equal
    residue numbers with distinct insertion codes are sequential in source order."""
    atoms = _ala("A", 184, "A", 0.0, 0.0, 0.0, src=0)
    atoms += _ala("A", 184, "", 3.0, 0.0, 0.0, src=1)
    s = _fake_structure(atoms)
    residues = Tier1Measurements._extract(s)
    assert [(r["_seq"], r["_icode"]) for r in residues] == [(184, "A"), (184, "")]
    assert Tier1Measurements._is_sequential(residues[0], residues[1]) is True


def test_skipped_insertion_code_is_not_sequential():
    """184A -> 184C skips 184B: the pair must NOT be treated as adjacent."""
    atoms = _ala("A", 184, "A", 0.0, 0.0, 0.0, src=0)
    atoms += _ala("A", 184, "C", 3.0, 0.0, 0.0, src=1)
    s = _fake_structure(atoms)
    residues = Tier1Measurements._extract(s)
    assert Tier1Measurements._is_sequential(residues[0], residues[1]) is False
    # Consecutive codes (A,B) are adjacent.
    atoms2 = _ala("A", 184, "A", 0.0, 0.0, 0.0, src=0)
    atoms2 += _ala("A", 184, "B", 3.0, 0.0, 0.0, src=1)
    r2 = Tier1Measurements._extract(_fake_structure(atoms2))
    assert Tier1Measurements._is_sequential(r2[0], r2[1]) is True


def test_label_seq_id_and_source_order_disagree_share_one_order():
    """Both engine consumers must order by the SAME canonical basis when
    label_seq_id and source_order disagree (label_seq_id wins)."""
    atoms = _ala("A", 1, "", 0.0, 0.0, 0.0, label_seq=2, src=0)
    atoms += _ala("A", 2, "", 3.0, 0.0, 0.0, label_seq=1, src=1)
    s = _fake_structure(atoms)
    ext = [(r["_seq"], r["_icode"]) for r in Tier1Measurements._extract(s)]
    seq = [(rs, ic) for rs, ic, _ in get_sequential_residues(s.atoms, "A")]
    assert ext == [(2, ""), (1, "")]
    assert seq == ext


def test_skipped_numbering_pair_cannot_pass_law110():
    """Positions 1,2,4,5: the 2->4 skip is incomplete adjacency accounting and
    must NOT yield an unqualified PASS from the measured 1->2 and 4->5 pairs."""
    atoms = _ala("A", 1, "", 0.0, 0.0, 0.0, src=0)
    atoms += _ala("A", 2, "", 3.0, 0.0, 0.0, src=1)
    atoms += _ala("A", 4, "", 6.0, 0.0, 0.0, src=2)
    atoms += _ala("A", 5, "", 9.0, 0.0, 0.0, src=3)
    s = _fake_structure(atoms)
    results, _, _ = Tier1Measurements.run_full_audit(s)
    assert results["LAW-110"]["status"] == "NOT_MEASURED"


def test_p2_source_backed_coverage_counts_absent_positions():
    """A missing internal source position must count against coverage: two
    observed complete residues at label_seq_id 1 and 3 of a 3-position polymer
    give 2/3 = 66.67%, not a fabricated 100%."""
    atoms = _ala("A", 1, "", 0.0, 0.0, 0.0, label_seq=1, src=0)
    atoms += _ala("A", 3, "", 3.0, 0.0, 0.0, label_seq=3, src=1)
    source_sequences = {"A": {"positions": 3, "basis": "label_seq_id"}}
    s = _fake_structure(atoms, source_sequences=source_sequences)
    results, coverage, is_exp = Tier1Measurements.run_full_audit(s)
    law105 = results["LAW-105"]
    assert law105["source_backed"] is True
    assert law105["coverage_basis"] == "source_polymer_positions"
    assert law105["coverage_denominator"] == 3
    assert law105["coverage_numerator"] == 2
    assert law105["observed"] == pytest.approx(66.67, abs=0.01)
    assert coverage == pytest.approx(66.67, abs=0.01)
    assert law105["status"] == "FAIL"   # below 70 because the absent position counts


def test_p2_source_backed_coverage_all_positions_present():
    atoms = _ala("A", 1, "", 0.0, 0.0, 0.0, label_seq=1, src=0)
    atoms += _ala("A", 2, "", 3.0, 0.0, 0.0, label_seq=2, src=1)
    source_sequences = {"A": {"positions": 2, "basis": "label_seq_id"}}
    s = _fake_structure(atoms, source_sequences=source_sequences)
    results, coverage, _ = Tier1Measurements.run_full_audit(s)
    law105 = results["LAW-105"]
    assert law105["source_backed"] is True
    assert law105["coverage_denominator"] == 2
    assert law105["coverage_numerator"] == 2
    assert law105["observed"] == 100.0
    assert law105["status"] == "PASS"


def test_p2_unavailable_source_is_non_authoritative():
    """Without authoritative source positions coverage is NOT_MEASURED and
    fail-closed; no observed-only numeric coverage becomes authoritative."""
    atoms = _ala("A", 1, "", 0.0, 0.0, 0.0, src=0)
    atoms += _ala("A", 2, "", 3.0, 0.0, 0.0, src=1)
    s = _fake_structure(atoms)   # no source_sequences
    results, coverage, _ = Tier1Measurements.run_full_audit(s)
    law105 = results["LAW-105"]
    assert law105["source_backed"] is False
    assert law105["coverage_basis"] == "unavailable"
    assert law105["observed"] is None
    assert law105["status"] == "NOT_MEASURED"
    assert coverage is None
