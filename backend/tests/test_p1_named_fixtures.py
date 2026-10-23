"""P1 named-fixture regressions: 1PPE and 5T4O residue ordering.

The auditor named 1PPE (insertion-coded residue 184A before base 184 in
source order) and 5T4O as the reproduction structures for the peptide
re-pairing defect. Alphabetical insertion-code sorting would re-order
these; the corrected contract must preserve source order. These tests pin
the real deposits through the Gemmi ingestion path.
"""
import os
from pathlib import Path
import pytest

os.environ.setdefault("TOSCANINI_API_KEYS", "test-key-xyz")
from tos.ingestion.processor import IngestionProcessor
from tos.engine.tier1_measurements import Tier1Measurements
from tos.engine.coord_math import get_sequential_residues

FIXTURES = Path(__file__).parent / "fixtures"


class Test1PPEInsertionOrder:

    def test_sequence_preserves_184a_before_184(self):
        pdb = (FIXTURES / "1PPE.pdb").read_bytes()
        structure = IngestionProcessor.run(pdb, "1PPE.pdb", "Audit", "experimental")
        # find the chain containing residue 184/184A
        found = None
        for chain in {a.chain_id for a in structure.atoms}:
            residues = get_sequential_residues(structure.atoms, chain)
            seqs = [(r[0], r[1]) for r in residues]
            if (184, "A") in seqs and (184, "") in seqs:
                found = seqs
                break
        assert found is not None, "1PPE must contain both 184A and 184"
        i184a = found.index((184, "A"))
        assert found[i184a + 1] == (184, ""), \
            "1PPE source order must place 184A before 184 (no alphabetical re-sort)"

    def test_1ppe_audits_without_fabrication(self):
        pdb = (FIXTURES / "1PPE.pdb").read_bytes()
        structure = IngestionProcessor.run(pdb, "1PPE.pdb", "Audit", "experimental")
        results, coverage, is_exp = Tier1Measurements.run_full_audit(structure)
        # coverage must be a finite value and the audit must not error
        assert isinstance(coverage, float) and 0.0 <= coverage <= 100.0
        law105 = results["LAW-105"]
        # if source sequence available, denominator is source-backed
        assert law105["coverage_denominator"] > 0
        # Exact measured-pair bookkeeping: LAW-110 must count real sequential
        # C-N peptide pairs and report them honestly (not skipped/fabricated).
        law110 = results["LAW-110"]
        assert law110.get("status") in ("PASS", "VETO")
        assert law110.get("sample", 0) > 0
        # Pinned on this candidate: 1PPE yields 238 measured sequential pairs
        assert law110["sample"] == 238
        # LAW-120 must measure backbone angle pairs for the same adjacency
        law120 = results["LAW-120"]
        assert law120.get("sample", 0) > 0
        assert law120.get("status") in ("PASS", "VETO", "NOT_MEASURED",
                                        "MEASUREMENT_FAILED")


class Test5T4OIngestion:

    def test_5t4o_ingests_and_measures(self):
        pdb = (FIXTURES / "5T4O.pdb").read_bytes()
        structure = IngestionProcessor.run(pdb, "5T4O.pdb", "Audit", "experimental")
        results, coverage, is_exp = Tier1Measurements.run_full_audit(structure)
        assert isinstance(coverage, float) and 0.0 <= coverage <= 100.0
        # measure a deterministic law exists and is honest
        assert "LAW-120" in results
        s = results["LAW-120"]["status"]
        assert s in ("PASS", "VETO", "NOT_MEASURED", "MEASUREMENT_FAILED",
                     "NOT_APPLICABLE")
        # Exact measured-pair bookkeeping for the large deposit.
        law110 = results["LAW-110"]
        assert law110.get("sample", 0) > 0
        assert law110["sample"] == 4755
        law120 = results["LAW-120"]
        assert law120.get("sample", 0) > 0
        assert law120["sample"] == 14288
