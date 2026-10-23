"""P2 regression tests: source-backed LAW-105 coverage denominator.

Coverage must use the authoritative polymer sequence (PDB SEQRES / mmCIF
entity_poly_seq), including polymer positions absent from the coordinate
records. When the source sequence cannot be established, LAW-105 must
fail closed (MEASUREMENT_FAILED / observed=None) instead of silently
manufacturing 100% coverage on observed coordinates alone.
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
          insertion_code="", b_iso=90.0, label_seq_id=None):
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
    )


def _backbone_atoms(count, chain_id="A", b_iso=90.0, start_seq=1):
    atoms = []
    for i in range(start_seq, start_seq + count):
        x = float(i - start_seq)
        atoms.append(_atom("N", "N", [x, 0.0, 0.0], "ALA", i,
                           chain_id=chain_id, b_iso=b_iso))
        atoms.append(_atom("CA", "C", [x + 1.2, 0.0, 0.0], "ALA", i,
                           chain_id=chain_id, b_iso=b_iso))
        atoms.append(_atom("C", "C", [x + 2.4, 0.0, 0.0], "ALA", i,
                           chain_id=chain_id, b_iso=b_iso))
        atoms.append(_atom("O", "O", [x + 2.4, 0.0, 1.2], "ALA", i,
                           chain_id=chain_id, b_iso=b_iso))
    return atoms


def _structure(atoms, source_by_chain=None, source_available=None,
               method="predicted"):
    from tos.ingestion.processor import Structure, ConfidenceProfile, IngestionMetadata
    conf = ConfidenceProfile([90.0] * len(atoms), is_experimental=(method != "predicted"),
                             method=method, resolution=None)
    md = IngestionMetadata(
        source_atom_records=len(atoms),
        selected_atoms=len(atoms),
        total_residues=len({(a.chain_id, a.res_seq) for a in atoms}),
        source_polymer_residues_by_chain=source_by_chain if source_by_chain is not None else {},
        source_sequence_available=bool(source_available),
    )
    return Structure(atoms=list(atoms), audit_id="", confidence=conf, metadata=md)


class TestSourceBackedDenominator:

    def test_denominator_uses_source_sequence_not_observed(self):
        """Observed coordinates cover 5 residues but the source sequence has
        6 positions -> coverage must be 5/6 = 83.33, not 100."""
        from tos.engine.tier1_measurements import Tier1Measurements
        atoms = _backbone_atoms(5)
        s = _structure(atoms, source_by_chain={"A": 6}, source_available=True)
        r, _, _ = Tier1Measurements.run_full_audit(s)
        law105 = r["LAW-105"]
        assert law105["coverage_denominator"] == 6
        assert law105["coverage_numerator"] == 5
        assert abs(law105["observed"] - (5 / 6 * 100)) < 0.01
        assert law105["coverage_basis"] == "source_polymer_sequence"

    def test_full_coverage_when_source_equals_observed(self):
        from tos.engine.tier1_measurements import Tier1Measurements
        atoms = _backbone_atoms(5)
        s = _structure(atoms, source_by_chain={"A": 5}, source_available=True)
        r, _, _ = Tier1Measurements.run_full_audit(s)
        law105 = r["LAW-105"]
        assert law105["coverage_denominator"] == 5
        assert law105["coverage_numerator"] == 5
        assert law105["observed"] == pytest.approx(100.0)
        assert law105["status"] == "PASS"

    def test_missing_source_sequence_fails_closed(self):
        """Without source-sequence metadata, LAW-105 must not claim 100% on
        observed coordinates alone: it fails closed (MEASUREMENT_FAILED,
        observed=None), never PASS. The returned coverage (consumed by the
        adjudicator/API/PDF) must be None — the same fail-closed state — not
        a stale observed percentage."""
        from tos.engine.tier1_measurements import Tier1Measurements
        atoms = _backbone_atoms(5)
        s = _structure(atoms, source_by_chain={}, source_available=False)
        r, coverage_returned, _ = Tier1Measurements.run_full_audit(s)
        law105 = r["LAW-105"]
        assert law105["observed"] is None
        assert law105["status"] == "MEASUREMENT_FAILED"
        assert law105["source_sequence_available"] is False
        assert law105["coverage_basis"] == "source_polymer_sequence_unavailable"
        # returned coverage agrees (fail-closed), never an observed-only %
        assert coverage_returned is None

    def test_explicit_source_and_method_reported(self):
        from tos.engine.tier1_measurements import Tier1Measurements
        atoms = _backbone_atoms(4)
        s = _structure(atoms, source_by_chain={"A": 4, "B": 3}, source_available=True)
        r, _, _ = Tier1Measurements.run_full_audit(s)
        law105 = r["LAW-105"]
        assert law105["coverage_denominator"] == 7
        assert law105["source_polymer_residues_by_chain"] == {"A": 4, "B": 3}
        assert law105["coverage_basis"] == "source_polymer_sequence"
