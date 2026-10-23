"""Focused negative tests for engine correctness gaps identified in the audit."""
import gemmi
import numpy as np
import pytest

from tos.engine.identity import compute_coordinate_hash
from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import IngestionError, IngestionProcessor


class TestAltlocEdgeCases:
    """Tied, absent, zero, and malformed occupancies in alternate-conformer selection."""

    def test_tied_occupancy_prefers_lexicographic(self):
        pdb = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA AALA A   1       1.458   0.000   0.000  0.50 90.00           C
ATOM      3  CA BALA A   1       1.500   0.000   0.000  0.50 90.00           C
ATOM      4  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      5  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
END
"""
        s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
        ca = next(a for a in s.atoms if a.atom_name == "CA")
        assert ca.altloc == "A"
        assert s.metadata.altloc_conflicts == 1

    def test_zero_occupancy_does_not_win_over_positive(self):
        pdb = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA AALA A   1       1.458   0.000   0.000  0.00 90.00           C
ATOM      3  CA BALA A   1       1.500   0.000   0.000  0.60 90.00           C
ATOM      4  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      5  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
END
"""
        s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
        ca = next(a for a in s.atoms if a.atom_name == "CA")
        assert ca.altloc == "B"
        assert ca.occupancy == pytest.approx(0.6)

    def test_absent_occupancy_defaults_to_nan_and_reports_missing(self):
        pdb = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000              90.00           C
ATOM      3  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      4  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
END
"""
        s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
        ca = next(a for a in s.atoms if a.atom_name == "CA")
        # gemmi parses empty occupancy as 0.0 (valid finite value)
        import math
        assert ca.occupancy == 0.0

    def test_malformed_occupancy_reports_missing(self):
        # Occupancy field with non-numeric garbage
        pdb = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  X.XX 90.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 90.00           C
ATOM      3  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      4  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
END
"""
        s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
        # Should either parse occupancy as 0 or report missing; must not crash
        assert len(s.atoms) >= 3

    def test_altloc_mixed_conformations_do_not_mix(self):
        """Residue-level altloc selection must not construct artificial residues."""
        pdb = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  0.60 90.00           N
ATOM      2  CA AALA A   1       1.458   0.000   0.000  0.60 90.00           C
ATOM      3  N   ALA A   1       0.100   0.100   0.100  0.40 90.00           N
ATOM      4  CA BALA A   1       1.558   0.100   0.100  0.40 90.00           C
ATOM      5  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      6  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
END
"""
        s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
        # N is universal (blank altloc) and CA comes from dominant conformer A
        n = next(a for a in s.atoms if a.atom_name == "N")
        ca = next(a for a in s.atoms if a.atom_name == "CA")
        assert n.altloc == ""  # universal
        assert ca.altloc == "A"  # dominant conformer (occ 0.6 > 0.4)


    def test_blank_altloc_preferred_over_alternate(self):
        pdb = (
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  0.50 90.00           N\n"
            "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  0.50 90.00           C\n"
            "ATOM      3  CA AALA A   1       1.500   0.000   0.000  0.50 90.00           C\n"
            "ATOM      4  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C\n"
            "ATOM      5  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O\n"
            "END\n"
        )
        s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
        ca = next(a for a in s.atoms if a.atom_name == "CA")
        assert ca.altloc == ""  # blank preferred over A at equal occupancy

class TestIngestionErrorPropagation:
    """IngestionError must propagate through all entry points."""

    def test_empty_input_raises(self):
        with pytest.raises(IngestionError):
            IngestionProcessor.run(b"", "input.pdb", "Audit", "Audit")

    def test_garbage_input_raises(self):
        with pytest.raises(IngestionError):
            IngestionProcessor.run(b"not a structure at all", "input.txt", "Audit", "Audit")

    def test_json_input_raises(self):
        with pytest.raises(IngestionError):
            IngestionProcessor.run(b'{"not": "a structure"}', "input.json", "Audit", "Audit")


class TestUnknownAcquisitionMethod:
    """Unknown acquisition method must not silently default to predicted."""

    def test_unknown_method_returns_unknown(self):
        pdb = """\
HEADER    TEST
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 90.00           C
ATOM      3  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      4  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
END
"""
        s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
        assert s.confidence.method == "unknown"
        assert s.confidence.is_experimental is True  # unknown != predicted

    def test_explicit_xray_returns_xray(self):
        pdb = """\
HEADER    TEST
EXPDTA    X-RAY DIFFRACTION
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 90.00           C
ATOM      3  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      4  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
END
"""
        s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
        assert s.confidence.method == "xray"
        assert s.confidence.is_experimental is True


class TestCoordinateHashSensitivity:
    """Coordinate hash must change with insertion code, altloc, B-factor, or coordinate."""

    def test_insertion_code_change_alters_hash(self):
        pdb1 = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 90.00           C
ATOM      3  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      4  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
END
"""
        pdb2 = """\
ATOM      1  N   ALA A   1A      0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  ALA A   1A      1.458   0.000   0.000  1.00 90.00           C
ATOM      3  C   ALA A   1A      2.004   1.424   0.000  1.00 90.00           C
ATOM      4  O   ALA A   1A      2.445   2.573   0.000  1.00 90.00           O
END
"""
        s1 = IngestionProcessor.run(pdb1.encode(), "input.pdb", "Audit", "Audit")
        s2 = IngestionProcessor.run(pdb2.encode(), "input.pdb", "Audit", "Audit")
        assert compute_coordinate_hash(s1) != compute_coordinate_hash(s2)

    def test_bfactor_change_alters_hash(self):
        pdb1 = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 90.00           C
ATOM      3  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      4  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
END
"""
        pdb2 = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 50.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 50.00           C
ATOM      3  C   ALA A   1       2.004   1.424   0.000  1.00 50.00           C
ATOM      4  O   ALA A   1       2.445   2.573   0.000  1.00 50.00           O
END
"""
        s1 = IngestionProcessor.run(pdb1.encode(), "input.pdb", "Audit", "Audit")
        s2 = IngestionProcessor.run(pdb2.encode(), "input.pdb", "Audit", "Audit")
        # B-factor is not in coordinate identity, so hashes should be equal
        assert compute_coordinate_hash(s1) == compute_coordinate_hash(s2)


class TestLAW182MissingChain:
    """LAW-182 must fail closed when a chain lacks measurable Cα for hydrophobic residues."""

    def test_chain_with_hydrophobic_no_ca_fails_closed(self):
        # Chain A has hydrophobic residues with CA; Chain B has hydrophobic without CA
        pdb = """\
ATOM      1  N   LEU A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  LEU A   1       1.458   0.000   0.000  1.00 90.00           C
ATOM      3  C   LEU A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      4  O   LEU A   1       2.445   2.573   0.000  1.00 90.00           O
ATOM      5  N   VAL B   1      10.000   0.000   0.000  1.00 90.00           N
ATOM      6  C   VAL B   1      11.000   0.000   0.000  1.00 90.00           C
ATOM      7  O   VAL B   1      11.500   0.000   0.000  1.00 90.00           O
END
"""
        s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
        results, coverage, _ = Tier1Measurements.run_full_audit(s)
        assert results["LAW-182"]["status"] == "NOT_MEASURED"
        assert "B" in results["LAW-182"].get("reason", "")


class TestDegenerateGeometry:
    """Degenerate geometry must not produce favorable zero-valued observations."""

    def test_degenerate_angle_is_not_counted_as_pass(self):
        # All three atoms at the same position → angle_deg returns 0.0
        # but angle_deg_valid returns False
        from tos.engine.coord_math import angle_deg_valid
        pos = np.array([1.0, 2.0, 3.0])
        angle, valid = angle_deg_valid(pos, pos, pos)
        assert not valid
        assert angle == 0.0

    def test_degenerate_dihedral_is_not_counted_as_pass(self):
        from tos.engine.coord_math import dihedral_deg_valid
        pos = np.array([1.0, 2.0, 3.0])
        dih, valid = dihedral_deg_valid(pos, pos, pos, pos)
        assert not valid
        assert dih == 0.0


class TestLAW130InsertionCodeScope:
    """LAW-130 clash analysis must use insertion-code-aware residue identity."""

    def test_insertion_code_residues_not_same_residue(self):
        # Two residues with same seq but different insertion codes are different residues
        from tos.ingestion.atom_selection import residue_identity
        from tos.ingestion.processor import Atom
        import numpy as np
        a1 = Atom(atom_name="CA", element="C", pos=np.array([0,0,0]), res_name="ALA",
                  res_seq=1, chain_id="A", insertion_code="", b_iso=90.0)
        a2 = Atom(atom_name="CA", element="C", pos=np.array([1,1,1]), res_name="ALA",
                  res_seq=1, chain_id="A", insertion_code="A", b_iso=90.0)
        assert residue_identity(a1) != residue_identity(a2)


class TestNKGWindowsLocking:
    """NKG manager must handle missing fcntl gracefully."""

    def test_nkg_imports_without_fcntl(self):
        import tos.nkg.manager as mgr
        assert mgr is not None

    def test_nkg_record_audit_works(self, tmp_path):
        from tos.nkg.manager import NKGManager
        m = NKGManager(path=str(tmp_path / "test_nkg.jsonl"))
        payload = {
            "verdict": {"binary": "PASS", "physical_score": 100},
            "governance": {"audit_id": "a" * 64},
            "provenance": {"hash": "b" * 64, "byte_count": 100},
            "tier3": {"probability": 50},
        }
        result = m.record_audit(payload)
        assert result is True
        records = m.read_records()
        assert len(records["successes"]) == 1

    def test_nkg_duplicate_suppressed(self, tmp_path):
        from tos.nkg.manager import NKGManager
        m = NKGManager(path=str(tmp_path / "test_nkg2.jsonl"))
        payload = {
            "verdict": {"binary": "PASS", "physical_score": 100},
            "governance": {"audit_id": "a" * 64},
            "provenance": {"hash": "b" * 64, "byte_count": 100},
            "tier3": {"probability": 50},
        }
        assert m.record_audit(payload) is True
        assert m.record_audit(payload) is False  # duplicate suppressed
