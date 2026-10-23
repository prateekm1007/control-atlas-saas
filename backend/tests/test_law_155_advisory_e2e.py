"""P5 regression tests: LAW-155 advisory-only, non-blocking behavior.

LAW-155 (voxel occupancy) is defined and advisory (D3 / Part-16), not
retired. It must never independently produce a deterministic VETO, must
not appear among deterministic blocking laws, and a missing/failed
LAW-155 must not be misrepresented as favorable evidence. Other
deterministic failures must still veto normally.
"""
import os
from pathlib import Path
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
from tos.governance.station_sop import LAW_CANON, RETIRED_LAWS


def _all_passing():
    m = {}
    for lid in LAW_CANON:
        m[lid] = {
            "observed": 0.0, "deviation": "0.0", "sample": 100,
            "status": "PASS",
        }
    m["LAW-105"] = {
        "observed": 100.0, "status": "PASS", "sample": 100,
        "coverage_basis": "source_polymer_sequence",
        "coverage_numerator": 100, "coverage_denominator": 100,
        "source_sequence_available": True,
    }
    return m


def _adjudicate(m):
    return adjudicate_laws(AdjudicationInput(
        raw_measurements=m, coverage=100.0, is_experimental=False,
        method_type="predicted", architecture="NONE",
    ))


class TestLAW155AdvisoryOnly:

    def test_not_retired(self):
        assert "LAW-155" not in RETIRED_LAWS

    def test_law155_adverse_alone_does_not_veto(self):
        m = _all_passing()
        m["LAW-155"] = {"observed": 0.1, "status": "VETO", "sample": 10}
        res = _adjudicate(m)
        assert res.verdict != "VETO"
        row = next(r for r in res.law_rows if r["law_id"] == "LAW-155")
        assert row["method"] == "heuristic"
        assert row["status"] == "FAIL (Advisory)"
        assert "LAW-155" not in res.failing_deterministic

    def test_deterministic_failure_still_vetoes(self):
        m = _all_passing()
        m["LAW-120"] = {"observed": 50.0, "status": "VETO", "sample": 50}
        m["LAW-155"] = {"observed": 0.1, "status": "VETO", "sample": 10}
        res = _adjudicate(m)
        assert res.verdict == "VETO"
        assert any("LAW-120" in f for f in res.failing_deterministic)
        assert not any("LAW-155" in f for f in res.failing_deterministic)

    def test_missing_law155_is_not_favorable(self):
        m = _all_passing()
        m["LAW-155"] = {
            "observed": None, "status": "NOT_MEASURED", "sample": 0,
            "reason": "no atomic coordinates available",
        }
        res = _adjudicate(m)
        row = next(r for r in res.law_rows if r["law_id"] == "LAW-155")
        assert row["status"] in ("NOT_MEASURED", "FAIL (Advisory)", "NOT_APPLICABLE")
        assert res.verdict != "VETO"


class TestLAW155ApiContract:

    def test_api_exposes_heuristic_method(self):
        r = client.post(
            "/ingest",
            data={"mode": "predicted", "candidate_id": "P5_155"},
            files={"file": ("1crn.pdb", FIXTURE.read_bytes(), "chemical/x-pdb")},
            headers=KEY,
        )
        assert r.status_code == 200, r.text
        law155 = next(
            (l for l in r.json().get("tier1", {}).get("laws", [])
             if l.get("law_id") == "LAW-155"),
            None,
        )
        assert law155 is not None
        assert law155.get("method") == "heuristic"
        assert law155.get("status") != "VETO"
