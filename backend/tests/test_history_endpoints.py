"""
Toscanini Phase B.3 — Audit History Endpoint Tests
Verifies /history/{email} and /history/ip/{ip} endpoints.

Run: TOSCANINI_DATA_DIR=/tmp/tos_b3_test pytest backend/tests/test_history_endpoints.py -v
"""
import os
import sys
import json
import pytest

# ── CRITICAL: set env BEFORE any storage module is imported ──────────────────
# audit_store.py and comparisons.py call .mkdir() at module level using this var.
# If it is not set before import, they fall back to /app/data and fail on host.
os.environ["TOSCANINI_DATA_DIR"] = "/tmp/tos_b3_history_test"
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

# Pre-create the base dir so module-level mkdir succeeds on first import
import pathlib
pathlib.Path("/tmp/tos_b3_history_test").mkdir(parents=True, exist_ok=True)


@pytest.fixture(autouse=True)
def clean_storage(tmp_path, monkeypatch):
    """Each test gets isolated fresh storage directories.

    Must patch module-level Path constants AFTER import — monkeypatch.setattr
    replaces the attribute on the already-imported module object, which is
    what the functions close over at call time (not import time).
    """
    import tos.storage.comparisons as comp_mod
    import tos.storage.audit_store as audit_mod

    comp_dir  = tmp_path / "comparisons"
    audit_dir = tmp_path / "audits"
    comp_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)

    # Patch the lazy getter functions so all calls in this test use tmp dirs
    monkeypatch.setattr(comp_mod,  "_get_storage_dir", lambda: comp_dir)
    monkeypatch.setattr(audit_mod, "_get_audit_dir",   lambda: audit_dir)
    yield


# ── Helpers ───────────────────────────────────────────────────────────────────

def _write_comparison(email: str, original_id: str, refined_id: str,
                      protocol: str = "openmm") -> None:
    from tos.storage.comparisons import store_comparison
    store_comparison(original_id, refined_id, {
        "user_email":  email,
        "protocol":    protocol,
        "status":      "complete",
        "refinement_method": "managed_gpu",
    })


def _write_audit(audit_id: str, verdict: str, score: int, coverage: float) -> None:
    from tos.storage.audit_store import store_audit_result
    store_audit_result(audit_id, {
        "verdict": {
            "binary":              verdict,
            "deterministic_score": score,
            "coverage_pct":        coverage,
        },
        "governance": {"audit_id": audit_id},
    })


# ── Group 1: list_user_comparisons (storage layer) ────────────────────────────

def test_empty_history_returns_empty_list():
    from tos.storage.comparisons import list_user_comparisons
    result = list_user_comparisons("nobody@test.com")
    assert result == []

def test_single_comparison_returned():
    from tos.storage.comparisons import list_user_comparisons
    _write_comparison("user@test.com", "BASE01", "REF01")
    result = list_user_comparisons("user@test.com")
    assert len(result) == 1
    assert result[0]["original_audit_id"] == "BASE01"
    assert result[0]["refined_audit_id"]  == "REF01"

def test_multiple_comparisons_sorted_newest_first():
    from tos.storage.comparisons import list_user_comparisons
    import time
    _write_comparison("order@test.com", "BASE01", "REF01")
    time.sleep(0.05)
    _write_comparison("order@test.com", "BASE02", "REF02")
    result = list_user_comparisons("order@test.com")
    assert len(result) == 2
    assert result[0]["original_audit_id"] == "BASE02"  # newest first

def test_email_isolation():
    """User A history must not contain User B entries."""
    from tos.storage.comparisons import list_user_comparisons
    _write_comparison("alice@test.com", "ABASE", "AREF")
    _write_comparison("bob@test.com",   "BBASE", "BREF")
    alice = list_user_comparisons("alice@test.com")
    bob   = list_user_comparisons("bob@test.com")
    assert len(alice) == 1
    assert len(bob)   == 1
    assert alice[0]["original_audit_id"] == "ABASE"
    assert bob[0]["original_audit_id"]   == "BBASE"


# ── Group 2: verdict enrichment ───────────────────────────────────────────────

def test_audit_store_roundtrip():
    from tos.storage.audit_store import store_audit_result, get_audit_result
    _write_audit("AUD01", "VETO", 42, 75.5)
    result = get_audit_result("AUD01")
    assert result["verdict"]["binary"]              == "VETO"
    assert result["verdict"]["deterministic_score"] == 42
    assert result["verdict"]["coverage_pct"]        == 75.5

def test_missing_audit_returns_none():
    from tos.storage.audit_store import get_audit_result
    assert get_audit_result("NONEXISTENT") is None


# ── Group 3: main.py source-level wiring ─────────────────────────────────────

def test_history_endpoint_registered():
    src = open(os.path.join(os.path.dirname(__file__), "../main.py")).read()
    assert '"/history/{user_email:path}"' in src

def test_history_ip_endpoint_registered():
    src = open(os.path.join(os.path.dirname(__file__), "../main.py")).read()
    assert '"/history/ip/{client_ip}"' in src

def test_no_hardcoded_data_path_in_health():
    src = open(os.path.join(os.path.dirname(__file__), "../main.py")).read()
    assert '"/app/data/comparisons"' not in src

def test_verdict_enrichment_present():
    src = open(os.path.join(os.path.dirname(__file__), "../main.py")).read()
    assert "refined_verdict"  in src
    assert "original_verdict" in src

def test_list_user_comparisons_imported_in_history():
    src = open(os.path.join(os.path.dirname(__file__), "../main.py")).read()
    assert "list_user_comparisons" in src
