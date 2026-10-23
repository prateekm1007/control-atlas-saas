"""P1 regression tests: residue ordering and peptide adjacency contract.

Source-file residue order is authoritative. Alphabetical insertion-code
sorting (chain, res_seq, insertion_code, res_name) re-pairs peptide geometry
for deposits whose source order places an insertion-coded residue before its
base-numbered residue (e.g. 1PPE: 183, 184A, 184, 185). These tests pin the
corrected ordering contract (P1): order by authoritative sequence position
(mmcCIF label_seq_id when valid, else res_seq) with a STABLE sort, so
same-position residues keep source order; never alphabetical insertion-code
sort.
"""
import os
from pathlib import Path
import numpy as np
import pytest

os.environ.setdefault("TOSCANINI_API_KEYS", "test-key-xyz")
import main as main_mod
import tos.security.auth as _auth

if not _auth.ALLOWED_KEYS:
    _auth.ALLOWED_KEYS = {"test-key-xyz"}


def _atom(atom_name, element, pos, res_name, res_seq, chain_id="A",
          insertion_code="", b_iso=90.0, label_seq_id=None, source_order=-1):
    from tos.ingestion.processor import Atom
    return Atom(
        atom_name=atom_name,
        element=element,
        pos=np.array(pos, dtype=float),
        res_name=res_name,
        res_seq=res_seq,
        chain_id=chain_id,
        insertion_code=insertion_code,
        b_iso=b_iso,
        label_seq_id=label_seq_id,
        source_order=source_order,
    )


def _backbone_atoms(seq_order, chain_id="A"):
    """Build N/CA/C/O atoms for residues in the exact order given by
    ``seq_order`` — a list of (res_seq, insertion_code). The source (file)
    order is the order of this list, mirroring what the parser sees."""
    atoms = []
    n = 0
    for res_seq, icode in seq_order:
        x = float(n)
        atoms.append(_atom("N", "N", [x, 0.0, 0.0], "ALA", res_seq,
                           chain_id=chain_id, insertion_code=icode, source_order=n))
        atoms.append(_atom("CA", "C", [x + 1.2, 0.0, 0.0], "ALA", res_seq,
                           chain_id=chain_id, insertion_code=icode, source_order=n + 1))
        atoms.append(_atom("C", "C", [x + 2.4, 0.0, 0.0], "ALA", res_seq,
                           chain_id=chain_id, insertion_code=icode, source_order=n + 2))
        atoms.append(_atom("O", "O", [x + 2.4, 0.0, 1.2], "ALA", res_seq,
                           chain_id=chain_id, insertion_code=icode, source_order=n + 3))
        n += 4
    return atoms


class TestGetSequentialResiduesOrder:

    def test_preserves_source_order_for_insertion_coded_residues(self):
        """1PPE pattern: source order 183, 184A, 184, 185 must NOT be
        re-sorted to 183, 184, 184A, 185 by alphabetical insertion-code
        sorting."""
        from tos.engine.coord_math import get_sequential_residues
        seq_order = [(183, ""), (184, "A"), (184, ""), (185, "")]
        atoms = _backbone_atoms(seq_order)
        residues = get_sequential_residues(atoms, "A")
        # (res_seq, insertion_code) in source order
        assert [(r[0], r[1]) for r in residues] == [
            (183, ""), (184, "A"), (184, ""), (185, "")
        ], "insertion-coded residue order must follow source file order"

    def test_label_seq_id_wins_when_present(self):
        """When a valid mmcIF label_seq_id is present it must drive order,
        even if auth res_seq numbers are shared (e.g. insertion-coded
        residues)."""
        from tos.engine.coord_math import get_sequential_residues
        # source order lists 184A first, but label_seq_id says 184 comes
        # before 184A
        atoms = [
            _atom("N", "N", [0.0, 0.0, 0.0], "ALA", 184, insertion_code="A",
                  label_seq_id=31, source_order=0),
            _atom("CA", "C", [1.2, 0.0, 0.0], "ALA", 184, insertion_code="A",
                  label_seq_id=31, source_order=1),
            _atom("N", "N", [2.4, 0.0, 0.0], "ALA", 184, insertion_code="",
                  label_seq_id=30, source_order=2),
            _atom("CA", "C", [3.6, 0.0, 0.0], "ALA", 184, insertion_code="",
                  label_seq_id=30, source_order=3),
        ]
        residues = get_sequential_residues(atoms, "A")
        # label_seq_id 30 (base 184) must sort before 31 (184A)
        assert [r[1] for r in residues] == ["", "A"]

    def test_distinct_seq_preserves_order(self):
        from tos.engine.coord_math import get_sequential_residues
        atoms = _backbone_atoms([(1, ""), (2, ""), (3, "")])
        residues = get_sequential_residues(atoms, "A")
        assert [r[0] for r in residues] == [1, 2, 3]


class TestExtractOrder:

    def test_extract_preserves_source_order_for_insertion_codes(self):
        from tos.engine.tier1_measurements import Tier1Measurements
        seq_order = [(183, ""), (184, "A"), (184, ""), (185, "")]
        atoms = _backbone_atoms(seq_order)
        residues = Tier1Measurements._extract(_structure(atoms))
        assert [(r["_seq"], r["_icode"]) for r in residues] == [
            (183, ""), (184, "A"), (184, ""), (185, "")
        ], "_extract must preserve insertion-coded source order"


def _structure(atoms):
    from tos.ingestion.processor import Structure, ConfidenceProfile
    conf = ConfidenceProfile([90.0] * len(atoms), is_experimental=False,
                             method="predicted", resolution=None)
    return Structure(atoms=list(atoms), audit_id="", confidence=conf)


class TestPeptideAdjacencyNoFalseVeto:

    def test_insertion_code_order_does_not_fabricate_peptide_assignment(self):
        """Build a structure where alphabetical sorting would re-pair the
        peptide bonds (base before insertion-coded residue). The corrected
        ordering keeps source order, so the C-N peptide distance gate is
        evaluated across the true bonded pairs, not a fabricated one. We
        assert the sequential residue order the engine uses, which is the
        precondition for honest LAW-120/135/160 measurement."""
        from tos.engine.tier1_measurements import Tier1Measurements
        seq_order = [(183, ""), (184, "A"), (184, ""), (185, "")]
        atoms = _backbone_atoms(seq_order)
        residues = Tier1Measurements._extract(_structure(atoms))
        # Peptide pairs are evaluated on consecutive residues in this order.
        order = [(r["_seq"], r["_icode"]) for r in residues]
        # The bonded pair (184A -> 184) must be consecutive in source order,
        # NOT re-paired to (184 -> 184A) by an alphabetical sort.
        assert (184, "A") in order
        i184a = order.index((184, "A"))
        assert order[i184a + 1] == (184, ""), \
            "insertion-coded residue must be immediately followed by its base residue"


class TestSharedAdjacencyContract:

    def _res(self, chain, seq, icode="", label_seq=None):
        return {
            "_chain": chain, "_seq": seq, "_icode": icode,
            "_label_seq": label_seq,
            "_has_label_seq": label_seq is not None,
            "_seq_pos": label_seq if label_seq is not None else seq,
        }

    def test_label_seq_consecutive_is_sequential(self):
        from tos.engine.tier1_measurements import Tier1Measurements
        r1 = self._res("A", 5, label_seq=10)
        r2 = self._res("A", 6, label_seq=11)
        assert Tier1Measurements._is_sequential(r1, r2) is True

    def test_label_seq_missing_internal_position_not_sequential(self):
        """label_seq gap > 1 means a missing internal polymer position; the
        pair must NOT be treated as sequential (no fabricated adjacency)."""
        from tos.engine.tier1_measurements import Tier1Measurements
        r1 = self._res("A", 5, label_seq=10)
        r3 = self._res("A", 7, label_seq=12)  # position 11 missing
        assert Tier1Measurements._is_sequential(r1, r3) is False

    def test_label_seq_reverse_not_sequential(self):
        from tos.engine.tier1_measurements import Tier1Measurements
        r1 = self._res("A", 5, label_seq=10)
        r2 = self._res("A", 6, label_seq=11)
        # reverse order (11 then 10) is not consecutive forward
        assert Tier1Measurements._is_sequential(r2, r1) is False

    def test_mixed_label_seq_availability_is_ambiguous(self):
        """One residue has label_seq_id, the other does not -> ambiguous -> not
        sequential (fail-closed)."""
        from tos.engine.tier1_measurements import Tier1Measurements
        r1 = self._res("A", 5, label_seq=10)
        r2 = self._res("A", 6)  # no label_seq (author fallback)
        assert Tier1Measurements._is_sequential(r1, r2) is False
        assert Tier1Measurements._is_sequential(r2, r1) is False

    def test_author_numbering_fallback_unchanged(self):
        from tos.engine.tier1_measurements import Tier1Measurements
        r1 = self._res("A", 1)
        r2 = self._res("A", 2)
        assert Tier1Measurements._is_sequential(r1, r2) is True
        # missing span under author numbering (1 -> 3) is not sequential
        r3 = self._res("A", 3)
        assert Tier1Measurements._is_sequential(r1, r3) is False
        # chain boundary is never sequential
        rB = self._res("B", 2)
        assert Tier1Measurements._is_sequential(r1, rB) is False

    def test_chain_boundary_is_not_sequential(self):
        from tos.engine.tier1_measurements import Tier1Measurements
        assert Tier1Measurements._is_sequential(
            self._res("A", 1, label_seq=1),
            self._res("B", 2, label_seq=2)) is False
