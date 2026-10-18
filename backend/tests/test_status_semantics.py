"""
WP-SCI-1 containment and status-semantics tests (directive §1-§2, §7).

Verifies, at the adjudication layer:
  - a measured LAW-120 VETO alone cannot produce an authoritative VETO
    (verdict INDETERMINATE with containment provenance on the row);
  - an independent deterministic veto still produces VETO even when a
    contained law also failed;
  - LAW-195 NOT_APPLICABLE neither passes nor vetoes;
  - missing mandatory measurements prevent authoritative PASS;
  - the measured verdict path through run_full_audit honors the same
    semantics end to end.
"""
import copy
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.engine.adjudicator import adjudicate_laws, AdjudicationInput
from tos.engine.tier1_measurements import Tier1Measurements
from tos.governance.station_sop import (
    LAW_CANON,
    CONTAINED_LAWS,
    STATUS_NOT_APPLICABLE,
    STATUS_NOT_MEASURED,
)
from tos.ingestion.processor import Atom, Structure, ConfidenceProfile


def _measurement(lid, observed, status, sample=5):
    return {"observed": observed, "deviation": "0.0", "sample": sample,
            "status": status}


def _all_pass_measurements():
    return {lid: _measurement(lid, 0.0 if LAW_CANON[lid]["operator"] == "<="
                              else 100.0, "PASS")
            for lid in LAW_CANON}


def _make_input(measurements, coverage=100.0, is_experimental=False,
                method_type="predicted"):
    return AdjudicationInput(
        raw_measurements=measurements, coverage=coverage,
        is_experimental=is_experimental, method_type=method_type,
        architecture="NONE",
    )


class TestContainment:

    def test_law120_is_contained_by_policy(self):
        assert "LAW-120" in CONTAINED_LAWS

    def test_contained_veto_alone_yields_indeterminate(self):
        m = _all_pass_measurements()
        m["LAW-120"] = _measurement("LAW-120", 42.0, "VETO")
        result = adjudicate_laws(_make_input(m))
        assert result.verdict == "INDETERMINATE"
        assert result.contained_law_failures == ["LAW-120"]
        row = next(r for r in result.law_rows if r["law_id"] == "LAW-120")
        assert row["contained"] is True
        assert row["containment_reason"] == "measurement-definition-repair"
        assert row["status"] == "VETO"  # measured value preserved, honestly

    def test_independent_veto_still_vetoes(self):
        m = _all_pass_measurements()
        m["LAW-120"] = _measurement("LAW-120", 42.0, "VETO")
        m["LAW-110"] = _measurement("LAW-110", 3, "VETO")
        result = adjudicate_laws(_make_input(m))
        assert result.verdict == "VETO"
        assert result.contained_law_failures == ["LAW-120"]

    def test_contained_law_not_counted_as_passed(self):
        """A contained law measured VETO must not inflate det_passed."""
        m = _all_pass_measurements()
        result_clean = adjudicate_laws(_make_input(m))
        m["LAW-120"] = _measurement("LAW-120", 42.0, "VETO")
        result_veto = adjudicate_laws(_make_input(m))
        assert result_veto.det_passed == result_clean.det_passed - 1


class TestStatusSemantics:

    def test_not_applicable_does_not_block_pass(self):
        m = _all_pass_measurements()
        m["LAW-195"] = _measurement("LAW-195", None, "NOT_APPLICABLE", sample=0)
        result = adjudicate_laws(_make_input(m))
        assert result.verdict == "PASS"
        assert result.indeterminate_causes == []
        # deterministic PASS count reflects the N/A law not passing
        assert result.det_passed == len(
            [r for r in result.law_rows
             if r["status"] == "PASS" and r["method"] == "deterministic"])

    def test_not_measured_blocks_authoritative_pass(self):
        m = _all_pass_measurements()
        m["LAW-125"] = _measurement("LAW-125", None, "NOT_MEASURED", sample=0)
        result = adjudicate_laws(_make_input(m))
        assert result.verdict == "INDETERMINATE"
        assert result.indeterminate_causes == ["LAW-125"]

    def test_measurement_failed_blocks_authoritative_pass(self):
        m = _all_pass_measurements()
        m["LAW-130"] = _measurement("LAW-130", None, "MEASUREMENT_FAILED", sample=0)
        result = adjudicate_laws(_make_input(m))
        assert result.verdict == "INDETERMINATE"
        assert "LAW-130" in result.indeterminate_causes

    def test_missing_law_row_defaults_to_not_measured(self):
        """A law absent from raw_measurements is NOT_MEASURED, not a
        fabricated 'measured and failed'."""
        m = _all_pass_measurements()
        del m["LAW-170"]
        result = adjudicate_laws(_make_input(m))
        assert result.verdict == "INDETERMINATE"
        row = next(r for r in result.law_rows if r["law_id"] == "LAW-170")
        assert row["status"] == "NOT_MEASURED"
        assert row["observed"] is None

    def test_row_provenance_fields_present(self):
        """Each response row can be checked against the stored threshold,
        operator, units and policy provenance."""
        result = adjudicate_laws(_make_input(_all_pass_measurements()))
        for row in result.law_rows:
            canon = LAW_CANON[row["law_id"]]
            assert row["threshold"] == canon["threshold"]
            assert row["operator"] == canon["operator"]
            assert row["units"] == canon["unit"]
            assert row["method"] in ("deterministic", "advisory_experimental",
                                     "heuristic")

    def test_modality_rules_appear_in_actual_adjudication(self):
        """PIL-CAL-03 reclassification must appear in the classified rows
        of a real adjudication, not only in loader tests."""
        m = _all_pass_measurements()
        m["LAW-100"] = _measurement("LAW-100", 6.0, "VETO")
        result = adjudicate_laws(_make_input(
            m, is_experimental=True, method_type="x_ray"))
        row = next(r for r in result.law_rows if r["law_id"] == "LAW-100")
        assert row["method"] == "advisory_experimental"
        assert row["status"] == "FAIL (Advisory)"
        assert result.verdict == "PASS"  # advisory cannot veto


class TestEndToEndVerdictPath:

    @staticmethod
    def _ideal_chain_with_distorted_ocn():
        """Two-residue backbone built with plain numpy (independent of the
        engine): N-CA-C and CA-C-N exactly at Engh-Huber ideals, correct
        1.329 A C-N peptide distance — but O(i) placed so the peptide
        O-C-N angle is ~40 deg instead of 123 deg. Only LAW-120 (and no
        other law) is violated by construction.
        """
        import math

        def d(a):
            r = math.radians(a)
            return np.array([math.cos(r), math.sin(r), 0.0])

        c1 = np.zeros(3)
        # distorted carbonyl oxygen: 40 deg from the N2 direction
        n2 = 1.329 * d(123.0)
        o1 = 1.231 * d(123.0 - 40.0)
        ca1 = 1.525 * d(123.0 + 116.6)
        n1 = ca1 + 1.470 * d(59.6 + 111.0)
        ca2 = n2 + 1.470 * d(30.0)
        c2 = ca2 + 1.525 * d(210.0 - 111.0)
        o2 = c2 + 1.231 * d(0.0)
        atoms = []
        for seq, (n, ca, c, o) in {1: (n1, ca1, c1, o1),
                                   2: (n2, ca2, c2, o2)}.items():
            atoms += [
                Atom(atom_name="N", element="N", pos=n, res_name="ALA",
                     res_seq=seq, chain_id="A", insertion_code="", b_iso=90.0),
                Atom(atom_name="CA", element="C", pos=ca, res_name="ALA",
                     res_seq=seq, chain_id="A", insertion_code="", b_iso=90.0),
                Atom(atom_name="C", element="C", pos=c, res_name="ALA",
                     res_seq=seq, chain_id="A", insertion_code="", b_iso=90.0),
                Atom(atom_name="O", element="O", pos=o, res_name="ALA",
                     res_seq=seq, chain_id="A", insertion_code="", b_iso=90.0),
            ]
        return atoms

    def test_distorted_ocn_measured_veto_contained_verdict(self):
        """Full path: run_full_audit -> measured LAW-120 VETO -> overall
        INDETERMINATE (containment), never an authoritative PASS, and no
        fabricated independent veto."""
        structure = Structure(
            atoms=self._ideal_chain_with_distorted_ocn(),
            audit_id="E2E",
            confidence=ConfidenceProfile([], False, "predicted"))
        raw, coverage, _ = Tier1Measurements.run_full_audit(structure)
        result = adjudicate_laws(AdjudicationInput(
            raw_measurements=raw, coverage=coverage,
            is_experimental=False, method_type="predicted",
            architecture="NONE"))
        law120 = next(r for r in result.law_rows if r["law_id"] == "LAW-120")
        assert law120["status"] == "VETO"
        assert law120["contained"] is True
        assert result.verdict == "INDETERMINATE"
        assert result.contained_law_failures == ["LAW-120"]
        # missing data never became a favorable zero
        assert law120["observed"] > 18.0 and law120["sample_size"] > 0
        # no uncontained deterministic law vetoed for this construction
        uncontained = [r["law_id"] for r in result.law_rows
                       if r["status"] in ("VETO", "FAIL")
                       and r["method"] == "deterministic"
                       and r["law_id"] != "LAW-120"]
        assert uncontained == []
