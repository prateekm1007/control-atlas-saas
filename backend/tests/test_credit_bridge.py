"""
Toscanini Phase B.3 — Credit Bridge Tests
Verifies that API key tier maps correctly to GPU credit allocation.

Run: TOSCANINI_DATA_DIR=/tmp/tos_b3_test pytest backend/tests/test_credit_bridge.py -v
No Docker required. No live Redis required.
"""
import os
import sys
import pytest

os.environ.setdefault("TOSCANINI_DATA_DIR", "/tmp/tos_b3_test")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def fresh_db(tmp_path_factory):
    """Give each test run a clean throwaway keys.db."""
    import shutil
    from tos.saas import key_store as ks
    tmp = tmp_path_factory.mktemp("keys")
    original = ks.DB_PATH
    ks.DB_PATH = tmp / "keys.db"
    ks.init_db()
    yield ks
    ks.DB_PATH = original


# ── Group 1: GPU_TIER_CREDITS structure ───────────────────────────────────────

def test_gpu_tier_credits_exists():
    from tos.saas.key_store import GPU_TIER_CREDITS
    assert "free"       in GPU_TIER_CREDITS
    assert "pro"        in GPU_TIER_CREDITS
    assert "enterprise" in GPU_TIER_CREDITS

def test_gpu_tier_ordering():
    from tos.saas.key_store import GPU_TIER_CREDITS
    assert GPU_TIER_CREDITS["free"]  < GPU_TIER_CREDITS["pro"]
    assert GPU_TIER_CREDITS["pro"]   < GPU_TIER_CREDITS["enterprise"]

def test_free_tier_matches_anon():
    """free tier GPU runs must equal BETA_CREDITS_ANON (3) for consistency."""
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../gpu_worker"))
    from tos.saas.key_store import GPU_TIER_CREDITS
    from worker.credits import BETA_CREDITS_ANON
    assert GPU_TIER_CREDITS["free"] == BETA_CREDITS_ANON

def test_tier_limits_has_gpu_runs():
    from tos.saas.key_store import TIER_LIMITS
    for tier in ("free", "pro", "enterprise"):
        assert "gpu_runs" in TIER_LIMITS[tier], f"{tier} missing gpu_runs"


# ── Group 2: get_gpu_allocation_for_key ───────────────────────────────────────

def test_unknown_key_returns_free_allocation(fresh_db):
    from tos.saas.key_store import get_gpu_allocation_for_key, GPU_TIER_CREDITS
    result = get_gpu_allocation_for_key("nonexistent_hash_xyz")
    assert result == GPU_TIER_CREDITS["free"]

def test_free_key_allocation(fresh_db):
    ks = fresh_db
    key = ks.create_key("free")
    key_hash = ks.hash_key(key)
    assert ks.get_gpu_allocation_for_key(key_hash) == ks.GPU_TIER_CREDITS["free"]

def test_pro_key_allocation(fresh_db):
    ks = fresh_db
    key = ks.create_key("pro")
    key_hash = ks.hash_key(key)
    assert ks.get_gpu_allocation_for_key(key_hash) == ks.GPU_TIER_CREDITS["pro"]

def test_enterprise_key_allocation(fresh_db):
    ks = fresh_db
    key = ks.create_key("enterprise")
    key_hash = ks.hash_key(key)
    assert ks.get_gpu_allocation_for_key(key_hash) == ks.GPU_TIER_CREDITS["enterprise"]

def test_revoked_key_returns_free_allocation(fresh_db):
    ks = fresh_db
    key = ks.create_key("pro")
    key_hash = ks.hash_key(key)
    ks.revoke_key(key_hash)
    assert ks.get_gpu_allocation_for_key(key_hash) == ks.GPU_TIER_CREDITS["free"]


# ── Group 3: get_tier_for_key ─────────────────────────────────────────────────

def test_tier_for_unknown_key(fresh_db):
    from tos.saas.key_store import get_tier_for_key
    assert get_tier_for_key("nonexistent") == "free"

def test_tier_for_pro_key(fresh_db):
    ks = fresh_db
    key = ks.create_key("pro")
    key_hash = ks.hash_key(key)
    assert ks.get_tier_for_key(key_hash) == "pro"

def test_tier_for_revoked_key_is_free(fresh_db):
    ks = fresh_db
    key = ks.create_key("enterprise")
    key_hash = ks.hash_key(key)
    ks.revoke_key(key_hash)
    assert ks.get_tier_for_key(key_hash) == "free"


# ── Group 4: main.py bridge wiring (source-level) ────────────────────────────

def test_main_imports_bridge():
    src = open(os.path.join(os.path.dirname(__file__), "../main.py")).read()
    assert "get_tier_for_key" in src
    assert "get_gpu_allocation_for_key" in src

def test_main_has_402_response():
    src = open(os.path.join(os.path.dirname(__file__), "../main.py")).read()
    assert "status_code=402" in src

def test_main_deducts_after_dispatch():
    src = open(os.path.join(os.path.dirname(__file__), "../main.py")).read()
    assert "deduct_credits(_identifier, protocol, job_id)" in src

def test_main_no_dir_bug():
    src = open(os.path.join(os.path.dirname(__file__), "../main.py")).read()
    assert '"audit_result" in dir()' not in src
