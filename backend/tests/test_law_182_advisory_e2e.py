"""P4 regression tests: LAW-182 advisory-only behavior end-to-end.

LAW-182 (hydrophobic burial) is a heuristic signal. It must be visibly
advisory, must never independently produce a deterministic VETO, must not
appear among deterministic blocking laws, and a missing/failed LAW-182
measurement must not be misrepresented as favorable evidence. Other
deterministic failures must still veto normally.
"""
import os
from pathlib import Path
import numpy as np
import pytest

os.environ.setdefault("TOSCANINI_API_KEYS", "test-key-xyz")
import main as main_mod
import tos.security.auth as _auth
from fastapi.testclient import TestClient

if not _auth.ALLOWED_KEYS:
    _auth.ALLOWED_KEYS = {"test-key-xyz"}

KEY = {"X-API-Key": "test-key-xyz"}
FIXTURE = Path(__file__).parent / "fixtures" / "real_pdb_1crn.pdb"
client = TestClient(main_mod.app, raise_server_exceptions=False)

from tos.engine.adjudicator import AdjudicationInput, adjudicate_laws
from tos.governance.station_sop import LAW_CANON


def _all_passing(expected_threshold=None):
    """A measurement dict where every deterministic law passes."""
    m = {}
    for lid in LAW_CANON:
        m[lid] = {
            "observed": 0.0,
            "deviation": "0.0",
            "sample": 100,
            "status": "PASS",
        }
    m["LAW-105"] = {
        "observed": 100.0, "status": "PASS", "sample": 100,
        "coverage_basis": "source_polymer_sequence",
        "coverage_numerator": 100, "coverage_denominator": 100,
        "source_sequence_available": True,
    }
    return m


def _adjudicate(measurements, coverage=100.0, is_exp=False, method="predicted"):
    return adjudicate_laws(AdjudicationInput(
        raw_measurements=measurements,
        coverage=coverage,
        is_experimental=is_exp,
        method_type=method,
        architecture="NONE",
    ))


class TestLAW182AdvisoryOnly:

    def test_law182_adverse_alone_does_not_veto(self):
        m = _all_passing()
        m["LAW-182"] = {
            "observed": 0.1,  # below the 0.3 threshold
            "status": "VETO", "sample": 50,
        }
        res = _adjudicate(m)
        # An adverse heuristic alone must not veto
        assert res.verdict != "VETO"
        row = next(r for r in res.law_rows if r["law_id"] == "LAW-182")
        assert row["method"] == "heuristic"
        assert row["status"] == "FAIL (Advisory)"
        # Must not appear among deterministic blocking laws
        assert "LAW-182" not in res.failing_deterministic

    def test_deterministic_failure_still_vetoes(self):
        m = _all_passing()
        m["LAW-120"] = {"observed": 50.0, "status": "VETO", "sample": 50}
        m["LAW-182"] = {"observed": 0.1, "status": "VETO", "sample": 50}
        res = _adjudicate(m)
        assert res.verdict == "VETO"
        assert any("LAW-120" in f for f in res.failing_deterministic)
        # LAW-182 still not a deterministic blocker
        assert not any("LAW-182" in f for f in res.failing_deterministic)

    def test_missing_law182_is_not_favorable(self):
        m = _all_passing()
        m["LAW-182"] = {
            "observed": None, "status": "NOT_MEASURED", "sample": 0,
            "reason": "hydrophobic residues lack CA coordinates",
        }
        res = _adjudicate(m)
        row = next(r for r in res.law_rows if r["law_id"] == "LAW-182")
        # missing heuristic evidence must not be represented as a favorable PASS
        assert row["status"] in ("NOT_MEASURED", "FAIL (Advisory)", "NOT_APPLICABLE")
        assert row["observed"] in (None, 0.0)
        # and must not veto
        assert res.verdict != "VETO"


class TestLAW182ApiContract:

    def test_api_response_exposes_method_and_status(self):
        r = client.post(
            "/ingest",
            data={"mode": "predicted", "candidate_id": "P4_182"},
            files={"file": ("1crn.pdb", FIXTURE.read_bytes(), "chemical/x-pdb")},
            headers=KEY,
        )
        assert r.status_code == 200, r.text
        body = r.json()
        law182 = next(
            (l for l in body.get("tier1", {}).get("laws", [])
             if l.get("law_id") == "LAW-182"),
            None,
        )
        assert law182 is not None, "LAW-182 must be present in the API response"
        assert law182.get("method") == "heuristic"
        # heuristic advisory must never serialize as a deterministic VETO
        assert law182.get("status") != "VETO" or law182.get("method") != "heuristic"
