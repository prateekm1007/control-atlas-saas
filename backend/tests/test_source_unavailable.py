"""
WP-4 — fabricated structures must never be scored.

The dispatcher raises SourceUnavailable when every source fails; /ingest
returns 503 naming the sources tried, with no deterministic_score in the
body and no NKG writes. The valid path is proven with the REAL 1CRN
structure from backend/tests/fixtures/ — no synthetic scoring anywhere.
"""
import json
import os
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

# auth.py reads env at import time — make this module self-sufficient when
# run standalone (the full suite sets this in test_api_key_guard first).
os.environ.setdefault("TOSCANINI_API_KEYS", "test-key-xyz")

import main as main_mod
import tos.security.auth as _auth
from tos.generation.dispatcher import GenerationDispatcher, SourceUnavailable

if not _auth.ALLOWED_KEYS:  # auth imported before env was set — patch defensively
    _auth.ALLOWED_KEYS = {"test-key-xyz"}

FIXTURE = Path(__file__).parent / "fixtures" / "real_pdb_1crn.pdb"
client = TestClient(main_mod.app, raise_server_exceptions=False)
KEY = {"X-API-Key": "test-key-xyz"}


def _nkg_record_count() -> int:
    """Count NKG records on disk in the current TOSCANINI_DATA_DIR."""
    import os
    nkg_dir = Path(os.environ.get("TOSCANINI_DATA_DIR", "/app/data")) / "nkg"
    if not nkg_dir.exists():
        return 0
    total = 0
    for f in nkg_dir.rglob("*"):
        if f.is_file() and f.suffix in (".json", ".jsonl"):
            total += sum(1 for _ in open(f, encoding="utf-8", errors="ignore"))
    return total


class TestDispatcher:

    def test_rcsb_exhaustion_raises_source_unavailable(self, monkeypatch):
        """Both RCSB URLs failing -> SourceUnavailable listing the sources."""
        import tos.generation.dispatcher as disp
        def fail(url, timeout):
            raise ConnectionError(f"simulated outage for {url}")
        monkeypatch.setattr(disp.requests, "get", fail)
        with pytest.raises(SourceUnavailable) as ei:
            GenerationDispatcher._acquire_rcsb("XXXX")
        assert any("files.rcsb.org" in s for s in ei.value.sources)
        assert ei.value.candidate_id == "XXXX"

    def test_afdb_exhaustion_raises_source_unavailable(self, monkeypatch):
        import tos.generation.dispatcher as disp
        def fail(url, timeout):
            raise ConnectionError(f"simulated outage for {url}")
        monkeypatch.setattr(disp.requests, "get", fail)
        with pytest.raises(SourceUnavailable) as ei:
            GenerationDispatcher._acquire_afdb("P99999")
        assert any("alphafold.ebi.ac.uk" in s for s in ei.value.sources)

    def test_mock_is_gone(self):
        """WP-4: the fabricated two-atom _mock must not exist."""
        assert not hasattr(GenerationDispatcher, "_mock")


class TestIngestRoute:

    def test_unknown_candidate_returns_503_without_score(self, monkeypatch, tmp_path):
        """Unknown candidate -> 503 naming sources, no deterministic_score,
        NKG record count unchanged."""
        import os
        os.environ["TOSCANINI_DATA_DIR"] = str(tmp_path / "wp4")
        def always_unavailable(candidate_id, fallback):
            raise SourceUnavailable(candidate_id, [
                "https://files.rcsb.org/download/ZZZZ.pdb",
                "https://alphafold.ebi.ac.uk/api/prediction/ZZZZ",
            ])
        monkeypatch.setattr(main_mod.GenerationDispatcher, "acquire", always_unavailable)

        before = _nkg_record_count()
        r = client.post("/ingest", data={"mode": "experimental", "candidate_id": "ZZZZ"},
                        headers=KEY)
        assert r.status_code == 503
        body_text = r.text
        assert "deterministic_score" not in body_text
        detail = r.json()["detail"]
        assert detail["error"] == "source_unavailable"
        assert detail["candidate_id"] == "ZZZZ"
        assert any("rcsb.org" in s for s in detail["sources_tried"])
        assert any("alphafold.ebi.ac.uk" in s for s in detail["sources_tried"])
        assert _nkg_record_count() == before, "NKG must not receive records for unavailable sources"

    def test_real_structure_still_scores(self, tmp_path):
        """Control: uploading the real 1CRN fixture still produces a verdict
        (the WP-4 change only blocks fabricated/unavailable sources)."""
        import os
        os.environ["TOSCANINI_DATA_DIR"] = str(tmp_path / "wp4-real")
        r = client.post(
            "/ingest",
            data={"mode": "experimental", "candidate_id": "1CRN_TEST"},
            files={"file": ("real_pdb_1crn.pdb", FIXTURE.read_bytes(), "chemical/x-pdb")},
            headers=KEY,
        )
        assert r.status_code == 200, r.text
        body = r.json()
        assert "verdict" in body
        assert body["verdict"].get("binary") in ("PASS", "VETO", "INDETERMINATE")
