"""Residue identity / order / adjacency contract tests (P1).

Covers the insertion-code residue-ordering defect: the engine must preserve
authoritative source order (not alphabetical insertion-code sorting) when
pairing peptide geometry, so deposits such as 1PPE/5T4O with insertion-coded
residues cannot fabricate backbone gaps or false vetoes.
"""
from tos.engine.coord_math import get_sequential_residues
from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.atom_selection import residue_source_order_key
from tos.ingestion.processor import IngestionProcessor


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

    # The peptide geometry must NOT fabricate a gap.
    results, _, _ = Tier1Measurements.run_full_audit(s)
    assert results["LAW-110"]["observed"] == 0
    assert results["LAW-110"]["status"] == "PASS"


def test_alphaical_insertion_sort_would_fabricate_gap():
    """Demonstrate the ordering-specific verdict change: the OLD (seq, icode)
    sort pairs 1<->1A and 1A<->2 using wrong termini and counts a gap; the
    source-order contract does not."""
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
    # Emulate the previous alphabetical (seq, icode) ordering.
    old_order = sorted(residues, key=lambda r: (r["_seq"], r["_icode"], r["_name"]))
    old_residue_ordered = [(r["_seq"], r["_icode"]) for r in old_order]
    assert old_residue_ordered == [(1, ""), (1, "A"), (2, "")]  # re-ordered

    def old_gap_count(ordered):
        gaps = 0
        for i in range(len(ordered) - 1):
            r1, r2 = ordered[i], ordered[i + 1]
            if Tier1Measurements._is_sequential(r1, r2):
                if "C" in r1["_atoms"] and "N" in r2["_atoms"]:
                    from tos.engine.coord_math import distance
                    if distance(r1["_atoms"]["C"], r2["_atoms"]["N"]) > 2.0:
                        gaps += 1
        return gaps

    assert old_gap_count(old_order) > 0      # false gap under old sort
    assert old_gap_count(residues) == 0      # no gap under source-order contract


def test_sequence_gap_is_not_fabricated_as_favorable():
    """A genuine missing residue (1 and 3, no 2) must not be turned into a
    favorable zero: the residues are non-sequential, so LAW-110 is honestly
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
    # Non-sequential residues are not paired as a favorable zero.
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
    assert law105["coverage_denominator"] == 2        # includes incomplete residue
    assert law105["coverage_numerator"] == 1          # only complete backbone
    assert law105["complete_backbone_residues"] == 1
    assert law105["observed"] < 100.0                 # never fabricates 100%


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
