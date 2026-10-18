"""
WP-SCI-1 end-to-end adjudication tests (directive §7).

Exercises the CORRECTED engine through the actual supported API pipeline
(FastAPI TestClient -> /ingest -> _run_physics_sync -> Tier1 ->
adjudicate_laws -> response payload), proving that:

  - the production verdict path uses the corrected measurement (not just
    the unit-tested functions);
  - every response row can be checked against the stored threshold,
    operator, units and policy provenance (canon + semantic + matrix
    hashes);
  - LAW-120 containment appears in the real response;
  - statuses obey the vocabulary and empty samples never PASS.
"""
import os
import sys
from pathlib import Path

os.environ.setdefault("TOSCANINI_API_KEYS", "test-key-xyz")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from fastapi.testclient import TestClient

import main as main_mod
from tos.governance.station_sop import LAW_CANON, compute_law_semantic_hash

REPO = Path(__file__).resolve().parents[2]
CRN = REPO / "backend" / "tests" / "fixtures" / "real_pdb_1crn.pdb"

client = TestClient(main_mod.app, raise_server_exceptions=False)
KEY = {"X-API-Key": "test-key-xyz"}

VALID_STATUSES = {"PASS", "VETO", "FAIL", "NOT_MEASURED", "NOT_APPLICABLE",
                  "MEASUREMENT_FAILED", "FAIL (Advisory)"}


@pytest.fixture(scope="module")
def ingest_response():
    resp = client.post(
        "/ingest",
        data={"mode": "experimental", "candidate_id": "1CRN"},
        files={"file": ("real_pdb_1crn.pdb", CRN.read_bytes(),
                        "chemical/x-pdb")},
        headers=KEY,
    )
    assert resp.status_code == 200, resp.text
    return resp.json()


class TestIngestVerdictPath:

    def test_verdict_and_fingerprint_present(self, ingest_response):
        body = ingest_response
        assert body["verdict"]["binary"] in ("PASS", "VETO", "INDETERMINATE")
        fp = body["governance"]["governance_fingerprint"]
        assert fp["canon_hash"] == "6a9cd4b4349b81de"
        assert fp["law_semantic_hash"] == compute_law_semantic_hash()
        assert len(fp["matrix_hash"]) == 64

    def test_rows_carry_checkable_provenance(self, ingest_response):
        rows = ingest_response["tier1"]["laws"]
        assert len(rows) == len(LAW_CANON)
        for row in rows:
            canon = LAW_CANON[row["law_id"]]
            assert row["threshold"] == canon["threshold"]
            assert row["operator"] == canon["operator"]
            assert row["units"] == canon["unit"]
            assert row["status"] in VALID_STATUSES
            # measured rows carry a positive sample; unmeasured rows do not
            if row["status"] in ("PASS", "VETO", "FAIL"):
                assert row["sample_size"] > 0, row["law_id"]
            if row["status"] in ("NOT_MEASURED", "NOT_APPLICABLE"):
                assert row["sample_size"] == 0
                assert row["observed"] is None

    def test_law120_contained_in_real_response(self, ingest_response):
        rows = {r["law_id"]: r for r in ingest_response["tier1"]["laws"]}
        assert rows["LAW-120"]["contained"] is True
        assert rows["LAW-120"]["containment_reason"] == \
            "measurement-definition-repair"
        # corrected metric on crambin: measured, well under the pinned
        # (unvalidated) 18.0° threshold
        assert rows["LAW-120"]["sample_size"] > 0
        assert rows["LAW-120"]["observed"] < 18.0

    def test_disulfide_law_measured_for_crambin(self, ingest_response):
        """Crambin has 3 disulfides: LAW-195 must be genuinely MEASURED,
        never NOT_APPLICABLE, for this structure."""
        rows = {r["law_id"]: r for r in ingest_response["tier1"]["laws"]}
        assert rows["LAW-195"]["status"] in ("PASS", "VETO")
        assert rows["LAW-195"]["sample_size"] >= 3

    def test_no_vacuous_passes_in_response(self, ingest_response):
        """WP-SCI-1: no row may claim a measured PASS with zero samples —
        the exact defect class found in the v23 recorded run."""
        for row in ingest_response["tier1"]["laws"]:
            if row["status"] == "PASS":
                assert row["sample_size"] > 0, row["law_id"]
