"""
backend/tests/test_b6_regression.py

B.6 Operational Hardening — Full Regression Suite
══════════════════════════════════════════════════════════════════════════════

Proves without manual steps:
  1. Rate limiting fires at ceiling (Group 1)
  2. Rate limit is per-key/per-audit (not global) (Group 1)
  3. Token blacklist still works (Group 2 — B.1 preservation)
  4. Per-audit MAX_CALLBACKS_PER_AUDIT cap still works (Group 2)
  5. /usage endpoint source-level wiring (Group 3)
  6. Usage widget graceful degradation (Group 4)
  7. PDF certificate builds valid bytes for all verdict paths (Group 5)
  8. PDF is deterministic (Group 5)
  9. PDF uses comparison_engine delta shape (Group 5)

Run:
    cd backend
    TOSCANINI_DATA_DIR=/tmp/tos_b6_test pytest tests/test_b6_regression.py -v
"""
import os
import sys
import time
import pytest

# Force test data dir BEFORE any tos imports
os.environ.setdefault("TOSCANINI_DATA_DIR", "/tmp/tos_b6_test")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(autouse=True)
def _isolated(tmp_path, monkeypatch):
    """Each test gets a clean isolated data dir."""
    monkeypatch.setenv("TOSCANINI_DATA_DIR", str(tmp_path))
    yield tmp_path


@pytest.fixture
def audit_id():
    return "TESTAUDIT_B6_001"


@pytest.fixture
def fresh_token(audit_id):
    from tos.security.tokens import create_refinement_token
    return create_refinement_token(audit_id, "test@toscanini.ai")


# ═══════════════════════════════════════════════════════════════════════════════
# GROUP 1: Rate Limiting (B.6 new behaviour)
# ═══════════════════════════════════════════════════════════════════════════════

class TestRateLimiting:

    def test_check_rate_limit_allows_new_key(self, audit_id):
        from tos.security.tokens import check_callback_rate_limit
        result = check_callback_rate_limit(audit_id)
        assert result["allowed"] is True
        assert result["attempts_in_window"] == 0
        assert result["limit"] > 0

    def test_increment_then_check(self, audit_id):
        from tos.security.tokens import (
            check_callback_rate_limit, _increment_rate_counter, RATE_LIMIT_MAX_PER_WINDOW
        )
        for _ in range(3):
            _increment_rate_counter(audit_id)
        result = check_callback_rate_limit(audit_id)
        assert result["attempts_in_window"] == 3
        assert result["allowed"] is True  # 3 < ceiling

    def test_rate_limit_fires_at_ceiling(self, monkeypatch, audit_id):
        from tos.security.tokens import (
            check_callback_rate_limit, _increment_rate_counter, RATE_LIMIT_MAX_PER_WINDOW
        )
        # Exhaust the window
        for _ in range(RATE_LIMIT_MAX_PER_WINDOW):
            _increment_rate_counter(audit_id)
        result = check_callback_rate_limit(audit_id)
        assert result["allowed"] is False
        assert result["attempts_in_window"] == RATE_LIMIT_MAX_PER_WINDOW

    def test_different_audit_ids_are_independent(self, monkeypatch):
        from tos.security.tokens import (
            check_callback_rate_limit, _increment_rate_counter, RATE_LIMIT_MAX_PER_WINDOW
        )
        id_a = "AUDIT_A_RATELIMIT"
        id_b = "AUDIT_B_RATELIMIT"
        for _ in range(RATE_LIMIT_MAX_PER_WINDOW):
            _increment_rate_counter(id_a)
        # id_a exhausted, id_b should still be allowed
        assert check_callback_rate_limit(id_a)["allowed"] is False
        assert check_callback_rate_limit(id_b)["allowed"] is True

    def test_consume_token_increments_rate_counter(self, fresh_token, audit_id):
        from tos.security.tokens import (
            validate_refinement_token, check_callback_rate_limit
        )
        validate_refinement_token(fresh_token, consume=True)
        result = check_callback_rate_limit(audit_id)
        assert result["attempts_in_window"] >= 1

    def test_rate_limit_blocks_consume(self, monkeypatch, audit_id, tmp_path):
        """Consume all rate slots, then verify next consume raises ValueError."""
        from tos.security.tokens import (
            create_refinement_token, validate_refinement_token,
            _increment_rate_counter, RATE_LIMIT_MAX_PER_WINDOW
        )
        # Pre-fill rate counter to ceiling
        for _ in range(RATE_LIMIT_MAX_PER_WINDOW):
            _increment_rate_counter(audit_id)

        # Create a fresh token for this audit_id
        tok = create_refinement_token(audit_id)
        with pytest.raises(ValueError, match="rate limit"):
            validate_refinement_token(tok, consume=True)

    def test_get_token_stats_includes_rate_info(self, audit_id):
        from tos.security.tokens import get_token_stats
        stats = get_token_stats(audit_id)
        assert "rate_limit" in stats
        assert "allowed" in stats["rate_limit"]
        assert "attempts_in_window" in stats["rate_limit"]
        assert "limit" in stats["rate_limit"]


# ═══════════════════════════════════════════════════════════════════════════════
# GROUP 2: B.1 Preservation — blacklist + per-audit cap still work
# ═══════════════════════════════════════════════════════════════════════════════

class TestB1Preservation:

    def test_single_use_blacklist_still_works(self, fresh_token, audit_id):
        from tos.security.tokens import validate_refinement_token
        validate_refinement_token(fresh_token, consume=True)
        with pytest.raises(ValueError, match="already been used"):
            validate_refinement_token(fresh_token, consume=True)

    def test_token_valid_without_consume(self, fresh_token):
        from tos.security.tokens import validate_refinement_token
        payload = validate_refinement_token(fresh_token, consume=False)
        assert "audit_id" in payload
        # Can call again without consume — not blacklisted
        payload2 = validate_refinement_token(fresh_token, consume=False)
        assert payload2["audit_id"] == payload["audit_id"]

    def test_create_returns_string(self, audit_id):
        from tos.security.tokens import create_refinement_token
        tok = create_refinement_token(audit_id, "user@test.com")
        assert isinstance(tok, str)
        assert len(tok) > 40

    def test_revoke_blacklists_token(self, fresh_token, audit_id):
        from tos.security.tokens import revoke_token, validate_refinement_token
        revoke_token(fresh_token, audit_id)
        with pytest.raises(ValueError, match="already been used"):
            validate_refinement_token(fresh_token, consume=False)

    def test_per_audit_callback_cap(self, audit_id):
        """MAX_CALLBACKS_PER_AUDIT distinct tokens for same audit → cap fires."""
        from tos.security.tokens import (
            create_refinement_token, validate_refinement_token,
            MAX_CALLBACKS_PER_AUDIT
        )
        for i in range(MAX_CALLBACKS_PER_AUDIT):
            tok = create_refinement_token(f"{audit_id}_sub{i}")
            # Each is a different audit_id — test the per-audit counter directly
        # Now use the same audit_id MAX times
        from tos.security.tokens import _blacklist_token
        for i in range(MAX_CALLBACKS_PER_AUDIT):
            _blacklist_token(f"faketoken{i}", audit_id, reason="used")
        tok_over = create_refinement_token(audit_id)
        with pytest.raises(ValueError, match="Maximum callbacks"):
            validate_refinement_token(tok_over, consume=False)


# ═══════════════════════════════════════════════════════════════════════════════
# GROUP 3: /usage endpoint source-level wiring
# ═══════════════════════════════════════════════════════════════════════════════

class TestUsageEndpointWiring:

    def _src(self):
        main_path = os.path.join(os.path.dirname(__file__), "../main.py")
        return open(main_path).read()

    def test_usage_endpoint_registered(self):
        assert '"/usage"' in self._src()

    def test_usage_returns_tier(self):
        assert '"tier"' in self._src() or "'tier'" in self._src()

    def test_usage_returns_credits_used(self):
        assert "credits_used" in self._src()

    def test_usage_returns_quotas(self):
        assert "quotas" in self._src()

    def test_usage_requires_api_key(self):
        src = self._src()
        assert "X-API-Key" in src
        assert "401" in src

    def test_usage_imports_key_store(self):
        src = self._src()
        assert "TIER_LIMITS" in src or "get_tier_for_key" in src

    def test_main_has_no_hardcoded_data_paths(self):
        src = self._src()
        assert '"/app/data/comparisons"' not in src
        assert '"/app/data/audits"' not in src


# ═══════════════════════════════════════════════════════════════════════════════
# GROUP 4: Usage widget graceful degradation
# ═══════════════════════════════════════════════════════════════════════════════

class TestUsageWidget:
    """Unit tests for usage_widget.py — no live brain required."""

    @pytest.fixture(autouse=True)
    def _widget_path(self):
        """Ensure dashboard dir is importable."""
        dashboard_dir = os.path.join(os.path.dirname(__file__), "../../dashboard")
        sys.path.insert(0, os.path.abspath(dashboard_dir))
        yield

    def test_widget_file_exists(self):
        dashboard_dir = os.path.join(os.path.dirname(__file__), "../../dashboard")
        assert os.path.exists(os.path.join(dashboard_dir, "usage_widget.py"))

    def test_widget_returns_none_on_connection_error(self, monkeypatch):
        import requests as req
        from usage_widget import _fetch_usage
        _fetch_usage.clear()

        def _fail(*a, **kw):
            raise req.exceptions.ConnectionError("brain offline")
        monkeypatch.setattr("usage_widget.requests.get", _fail)
        assert _fetch_usage("tos_any_key") is None

    def test_widget_returns_none_on_401(self, monkeypatch):
        from usage_widget import _fetch_usage
        _fetch_usage.clear()

        class _R:
            status_code = 401
        monkeypatch.setattr("usage_widget.requests.get", lambda *a, **kw: _R())
        assert _fetch_usage("tos_bad_key") is None

    def test_widget_returns_data_on_200(self, monkeypatch):
        from usage_widget import _fetch_usage
        _fetch_usage.clear()

        expected = {"tier": "pro", "credits_used": 5, "credits_total": 50, "quotas": {}}

        class _R:
            status_code = 200
            def json(self): return expected
        monkeypatch.setattr("usage_widget.requests.get", lambda *a, **kw: _R())
        assert _fetch_usage("tos_good_key") == expected

    def test_widget_empty_key_returns_none(self):
        from usage_widget import _fetch_usage
        _fetch_usage.clear()
        assert _fetch_usage("") is None

    def test_widget_uses_backend_url_not_brain_url(self):
        widget_path = os.path.join(
            os.path.dirname(__file__), "../../dashboard/usage_widget.py"
        )
        widget_src = open(widget_path).read()
        assert "BACKEND_URL" in widget_src
        # Verify env var lookup uses BACKEND_URL not BRAIN_URL
        # Parse by splitting on literal string — no regex needed
        parts = widget_src.split("os.environ.get(")
        env_keys = []
        for part in parts[1:]:
            inner = part.split(")")[0]
            key = inner.split(",")[0].strip().strip("'\"")
            env_keys.append(key)
        assert "BACKEND_URL" in env_keys, f"BACKEND_URL missing from env lookups: {env_keys}"
        assert "BRAIN_URL" not in env_keys, f"BRAIN_URL found in env lookups: {env_keys}"

    def test_widget_has_cache_decorator(self):
        widget_path = os.path.join(
            os.path.dirname(__file__), "../../dashboard/usage_widget.py"
        )
        src = open(widget_path).read()
        assert "st.cache_data" in src

    def test_widget_graceful_timeout(self):
        widget_path = os.path.join(
            os.path.dirname(__file__), "../../dashboard/usage_widget.py"
        )
        src = open(widget_path).read()
        assert "timeout" in src


# ═══════════════════════════════════════════════════════════════════════════════
# GROUP 5: PDF Certificate
# ═══════════════════════════════════════════════════════════════════════════════

class TestPDFCertificate:

    @pytest.fixture(autouse=True)
    def _pdf_path(self):
        dashboard_dir = os.path.join(os.path.dirname(__file__), "../../dashboard")
        sys.path.insert(0, os.path.abspath(dashboard_dir))
        yield

    @pytest.fixture
    def improvement_delta(self):
        """Matches compare_audits() output shape exactly."""
        return {
            "baseline_audit_id":      "BASE_001",
            "refined_audit_id":       "REF_001",
            "verdict_change":         "VETO → PASS",
            "verdict_improved":       True,
            "coverage_delta":         12.5,
            "violation_count_before": 8,
            "violation_count_after":  2,
            "violation_count_delta":  -6,
            "improvements":           ["LAW-125", "LAW-130", "LAW-155"],
            "regressions":            [],
            "law_changes":            [],
        }

    @pytest.fixture
    def regression_delta(self):
        return {
            "baseline_audit_id":      "BASE_002",
            "refined_audit_id":       "REF_002",
            "verdict_change":         "PASS → VETO",
            "verdict_improved":       False,
            "coverage_delta":         -3.0,
            "violation_count_before": 2,
            "violation_count_after":  5,
            "violation_count_delta":  3,
            "improvements":           [],
            "regressions":            ["LAW-145", "LAW-182"],
            "law_changes":            [],
        }

    @pytest.fixture
    def neutral_delta(self):
        return {
            "baseline_audit_id":      "BASE_003",
            "refined_audit_id":       "REF_003",
            "verdict_change":         "VETO → VETO",
            "verdict_improved":       False,
            "coverage_delta":         0.0,
            "violation_count_before": 4,
            "violation_count_after":  4,
            "violation_count_delta":  0,
            "improvements":           [],
            "regressions":            [],
            "law_changes":            [],
        }

    @pytest.fixture
    def meta(self):
        return {
            "structure_name":  "6LU7_chain_A.pdb",
            "api_key_masked":  "tos_****abcd",
            "tier":            "pro",
            "engine":          "openmm",
            "issued_at":       1704067200,
        }

    def test_pdf_file_exists(self):
        dashboard_dir = os.path.join(os.path.dirname(__file__), "../../dashboard")
        assert os.path.exists(os.path.join(dashboard_dir, "pdf_certificate.py"))

    def test_improvement_returns_bytes(self, improvement_delta, meta):
        from pdf_certificate import build_certificate
        result = build_certificate(improvement_delta, meta)
        assert isinstance(result, bytes)
        assert len(result) > 100

    def test_regression_returns_bytes(self, regression_delta, meta):
        from pdf_certificate import build_certificate
        result = build_certificate(regression_delta, meta)
        assert isinstance(result, bytes)

    def test_neutral_returns_bytes(self, neutral_delta, meta):
        from pdf_certificate import build_certificate
        result = build_certificate(neutral_delta, meta)
        assert isinstance(result, bytes)

    def test_pdf_magic_bytes_or_utf8(self, improvement_delta, meta):
        """Accept PDF (%PDF) or plaintext fallback (reportlab absent)."""
        from pdf_certificate import build_certificate
        result = build_certificate(improvement_delta, meta)
        assert result[:4] == b"%PDF" or b"TOSCANINI" in result[:200]

    def test_deterministic_output(self, improvement_delta, meta):
        from pdf_certificate import build_certificate, _suppress_timestamps
        a = build_certificate(improvement_delta, meta)
        b = build_certificate(improvement_delta, meta)
        # _suppress_timestamps already applied inside build_certificate.
        # Verify the outputs are identical after any residual timestamp removal.
        assert _suppress_timestamps(a) == _suppress_timestamps(b), (
            "PDF output must be byte-identical for same inputs after timestamp suppression"
        )

    def test_baseline_id_in_output(self, improvement_delta, meta):
        from pdf_certificate import build_certificate
        result = build_certificate(improvement_delta, meta)
        assert b"BASE_001" in result

    def test_no_pdb_coords_in_output(self, improvement_delta, meta):
        """Certificate must contain no ATOM/HETATM lines."""
        from pdf_certificate import build_certificate
        result = build_certificate(improvement_delta, meta)
        text = result.decode("latin-1")
        assert "ATOM  " not in text
        assert "HETATM" not in text

    def test_improvement_cert_contains_certified_word(self, improvement_delta, meta):
        from pdf_certificate import build_certificate
        result = build_certificate(improvement_delta, meta)
        # PDF text streams are compressed; check /Keywords metadata (uncompressed)
        # OR check the plaintext fallback path
        assert b"CERTIFIED" in result, (
            "CERTIFIED must appear in /Keywords metadata field (uncompressed)"
        )

    def test_regression_cert_contains_regression_word(self, regression_delta, meta):
        from pdf_certificate import build_certificate
        result = build_certificate(regression_delta, meta)
        # PDF text streams are compressed; check /Keywords metadata (uncompressed)
        assert b"REGRESSION" in result, (
            "REGRESSION must appear in /Keywords metadata field (uncompressed)"
        )

    def test_uses_comparison_engine_delta_shape(self):
        """Verify pdf_certificate uses the same field names as compare_audits()."""
        pdf_src = open(os.path.join(
            os.path.dirname(__file__), "../../dashboard/pdf_certificate.py"
        )).read()
        for field in ["baseline_audit_id", "refined_audit_id", "violation_count_before",
                      "violation_count_after", "verdict_improved", "improvements", "regressions"]:
            assert field in pdf_src, f"Missing field reference: {field}"
