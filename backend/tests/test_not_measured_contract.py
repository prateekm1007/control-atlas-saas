"""
Round-24: the explicit NOT_MEASURED state, end to end.

An explicitly unmeasured measurement must be represented truthfully across the
whole contract: the engine row, the adjudicator (neither a pass nor a
scientific failure) and the API response schema (a nullable observation with
an explicit `measured` flag). These tests prove the whole chain.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.engine.adjudicator import AdjudicationInput, adjudicate_laws
from tos.engine.tier1_measurements import Tier1Measurements
from tos.governance.station_sop import LAW_CANON
from tos.ingestion.processor import Atom, ConfidenceProfile, Structure
from tos.schemas.response_v1 import PhysicsMeasurement

import numpy as np


def _atom(name, element, pos, res_name="ALA", res_seq=1, chain_id="A",
          insertion_code="", b_iso=90.0):
    return Atom(atom_name=name, element=element, pos=np.array(pos, dtype=float),
                res_name=res_name, res_seq=res_seq, chain_id=chain_id,
                insertion_code=insertion_code, b_iso=b_iso)


def _unmeasurable_atoms():
    """One residue: no sequence-adjacent C-alpha pair exists at all."""
    return [
        _atom("N",  "N", [0.0, 0.0, 0.0],  res_name="ALA", res_seq=1),
        _atom("CA", "C", [1.47, 0.0, 0.0], res_name="ALA", res_seq=1),
        _atom("C",  "C", [2.5, 1.0, 0.0],  res_name="ALA", res_seq=1),
    ]


def _audit(atoms):
    struct = Structure(atoms=atoms, audit_id="TEST",
                       confidence=ConfidenceProfile([], False, "predicted"))
    results, coverage, _ = Tier1Measurements.run_full_audit(struct)
    return results, coverage


def _all_passing(overrides=None):
    m = {}
    for lid in LAW_CANON:
        m[lid] = {
            "observed": LAW_CANON[lid]["threshold"] * 0.5,
            "deviation": "0.0",
            "sample": 100,
            "status": "PASS",
        }
    for lid, row in (overrides or {}).items():
        m[lid] = row
    return m


def _input(measurements, coverage=100.0):
    return AdjudicationInput(
        raw_measurements=measurements, coverage=coverage,
        is_experimental=False, method_type="predicted", architecture="NONE",
    )


def test_engine_returns_an_explicit_not_measured_row():
    """The engine row is explicit: NOT_MEASURED, measured False, null, sample 0."""
    results, _ = _audit(_unmeasurable_atoms())
    row = results["LAW-160"]
    assert row["status"] == "NOT_MEASURED"
    assert row["measured"] is False
    assert row["observed"] is None
    assert row["deviation"] is None
    assert row["sample"] == 0


def test_the_response_model_serializes_the_unmeasured_row():
    """The API response model accepts the unmeasured row and keeps the numbers."""
    results, _ = _audit(_unmeasurable_atoms())
    m = PhysicsMeasurement(
        law_id="LAW-160", title=LAW_CANON["LAW-160"]["title"],
        status=results["LAW-160"]["status"], method="deterministic",
        measured=results["LAW-160"]["measured"],
        observed=results["LAW-160"]["observed"],
        threshold=LAW_CANON["LAW-160"]["threshold"],
        operator=LAW_CANON["LAW-160"]["operator"],
        units=LAW_CANON["LAW-160"]["unit"],
        deviation=results["LAW-160"]["deviation"],
        sample_size=results["LAW-160"]["sample"],
        scope=LAW_CANON["LAW-160"]["scope"],
    )
    assert m.status == "NOT_MEASURED"
    assert m.measured is False
    assert m.observed is None
    dumped = m.model_dump()
    assert dumped["observed"] is None
    assert dumped["measured"] is False
    # an ordinary measured row still validates with a number
    ordinary = PhysicsMeasurement(
        law_id="LAW-160", title="Chain Integrity", status="PASS", method="deterministic",
        measured=True, observed=3.8, threshold=4.5, operator="<=", units="A",
        deviation="0.0", sample_size=2, scope="sequential C-alphas",
    )
    assert ordinary.observed == 3.8
    assert ordinary.measured is True


def test_an_unmeasured_deterministic_law_yields_indeterminate_not_pass():
    """Otherwise-passing deterministic laws + one NOT_MEASURED => INDETERMINATE."""
    measurements = _all_passing({
        "LAW-160": {"observed": None, "deviation": None, "measured": False,
                    "status": "NOT_MEASURED", "sample": 0},
    })
    res = adjudicate_laws(_input(measurements))
    assert res.verdict == "INDETERMINATE"
    assert "LAW-160" not in res.failing_deterministic
    assert not any("LAW-160" in f for f in res.failing_deterministic)
    # an unmeasured observation never earns a passing point
    row = [r for r in res.law_rows if r["law_id"] == "LAW-160"][0]
    assert row["status"] == "NOT_MEASURED"
    assert row["measured"] is False
    assert row["observed"] is None
    assert row["sample_size"] == 0
    assert res.det_passed == sum(
        1 for r in res.law_rows if r["status"] == "PASS" and r["method"] == "deterministic")


def test_a_real_deterministic_veto_is_retained_next_to_an_unmeasured_law():
    """A genuine veto is still a veto when another law is NOT_MEASURED."""
    measurements = _all_passing({
        "LAW-160": {"observed": None, "deviation": None, "measured": False,
                    "status": "NOT_MEASURED", "sample": 0},
        "LAW-120": {"observed": LAW_CANON["LAW-120"]["threshold"] * 5.0,
                    "deviation": "999.0", "sample": 100, "status": "VETO"},
    })
    res = adjudicate_laws(_input(measurements))
    assert res.verdict == "VETO"
    assert any("LAW-120" in f for f in res.failing_deterministic)
    assert not any("LAW-160" in f for f in res.failing_deterministic)


def test_the_coverage_gate_still_yields_indeterminate():
    """The coverage gate keeps its INDETERMINATE behaviour."""
    measurements = _all_passing({
        "LAW-160": {"observed": None, "deviation": None, "measured": False,
                    "status": "NOT_MEASURED", "sample": 0},
    })
    res = adjudicate_laws(_input(measurements, coverage=10.0))
    assert res.verdict == "INDETERMINATE"
