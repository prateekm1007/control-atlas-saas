"""
Tests for PIL-CAL-03 modality matrix loader.
Deterministic: tests JSON rules, not main.py wiring.
"""
import pytest
from tos.governance.modality_matrix import (
    compute_matrix_hash,
    get_matrix_meta,
    resolve_method,
)


class TestResolveMethod:

    def test_predicted_structure_is_never_reclassified(self):
        for law in ("LAW-100", "LAW-120", "LAW-160", "LAW-170"):
            result = resolve_method(law, False, "predicted")
            assert result == "deterministic", (
                f"{law} should stay deterministic for predicted, got {result}"
            )

    def test_law100_advisory_for_all_experimental(self):
        for method in ("x_ray", "cryo_em", "nmr"):
            result = resolve_method("LAW-100", True, method)
            assert result == "advisory_experimental", (
                f"LAW-100 should be advisory for {method}, got {result}"
            )

    def test_law160_advisory_for_all_experimental(self):
        for method in ("x_ray", "cryo_em", "nmr"):
            result = resolve_method("LAW-160", True, method)
            assert result == "advisory_experimental", (
                f"LAW-160 should be advisory for {method}, got {result}"
            )

    def test_law170_advisory_for_cryo_em(self):
        assert resolve_method("LAW-170", True, "cryo_em") == "advisory_experimental"

    def test_law170_deterministic_for_xray(self):
        assert resolve_method("LAW-170", True, "x_ray") == "deterministic"

    def test_law170_advisory_for_nmr(self):
        assert resolve_method("LAW-170", True, "nmr") == "advisory_experimental"

    def test_law125_advisory_for_nmr_only(self):
        assert resolve_method("LAW-125", True, "nmr") == "advisory_experimental"
        assert resolve_method("LAW-125", True, "cryo_em") == "deterministic"
        assert resolve_method("LAW-125", True, "x_ray") == "deterministic"

    def test_law150_advisory_for_nmr_only(self):
        assert resolve_method("LAW-150", True, "nmr") == "advisory_experimental"
        assert resolve_method("LAW-150", True, "cryo_em") == "deterministic"
        assert resolve_method("LAW-150", True, "x_ray") == "deterministic"

    def test_law120_never_reclassified(self):
        for method in ("x_ray", "cryo_em", "nmr"):
            result = resolve_method("LAW-120", True, method)
            assert result == "deterministic", (
                f"LAW-120 must stay deterministic for {method}, got {result}"
            )

    def test_heuristic_default_preserved_when_no_rule(self):
        result = resolve_method("LAW-999", True, "x_ray", default_method="heuristic")
        assert result == "heuristic"

    def test_unknown_law_returns_default(self):
        result = resolve_method("LAW-FAKE", True, "nmr", default_method="deterministic")
        assert result == "deterministic"


class TestMatrixHash:

    def test_hash_is_stable_across_calls(self):
        h1 = compute_matrix_hash()
        h2 = compute_matrix_hash()
        assert h1 == h2

    def test_hash_is_sha256_length(self):
        h = compute_matrix_hash()
        assert len(h) == 64, f"Expected 64-char SHA-256 hex, got {len(h)}"

    def test_hash_matches_known_value(self):
        expected = "ed1731ace09a9a3bd24d3659ab6f3fd184da2a2335d935a234fd6a9c14880aa4"
        actual = compute_matrix_hash()
        assert actual == expected, (
            f"Matrix hash drifted!\n  Expected: {expected}\n  Actual:   {actual}\n"
            "If intentional, update this test."
        )


class TestMatrixMeta:

    def test_meta_contains_required_fields(self):
        meta = get_matrix_meta()
        for field in ("schema_version", "policy_ref", "policy_version", "locked_by"):
            assert field in meta, f"Missing meta field: {field}"

    def test_policy_ref_is_PIL_CAL_03(self):
        assert get_matrix_meta()["policy_ref"] == "PIL-CAL-03"

    def test_schema_version_is_semver(self):
        parts = get_matrix_meta()["schema_version"].split(".")
        assert len(parts) == 3, "schema_version should be semver (x.y.z)"
        assert all(p.isdigit() for p in parts), "semver parts must be numeric"
