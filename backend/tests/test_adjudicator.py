"""
Tests for the pure Adjudicator engine.
No HTTP, no files, no Structure objects — just measurements in, verdicts out.
"""
import pytest
from tos.engine.adjudicator import (
    adjudicate_laws,
    AdjudicationInput,
    AdjudicationResult,
    StrategicMath,
    _compute_verdict,
)
from tos.governance.station_sop import (
    LAW_CANON,
    DETERMINISTIC_COUNT,
    HEURISTIC_COUNT,
    LAW_METHOD_CLASSIFICATIONS,
)


def _make_input(
    measurements=None,
    coverage=100.0,
    is_experimental=False,
    method_type="predicted",
    architecture="NONE",
):
    """Helper to build AdjudicationInput with sensible defaults."""
    if measurements is None:
        measurements = {}
    return AdjudicationInput(
        raw_measurements=measurements,
        coverage=coverage,
        is_experimental=is_experimental,
        method_type=method_type,
        architecture=architecture,
    )


def _all_passing_measurements():
    """Build a measurements dict where every law passes."""
    m = {}
    for lid in LAW_CANON:
        m[lid] = {
            "observed": LAW_CANON[lid]["threshold"] * 0.5,
            "deviation": "0.0",
            "sample": 100,
            "status": "PASS",
        }
    return m


def _single_veto_measurements(veto_lid):
    """All laws pass except the specified one, which VETOes."""
    m = _all_passing_measurements()
    m[veto_lid] = {
        "observed": LAW_CANON[veto_lid]["threshold"] * 5.0,
        "deviation": "999.0",
        "sample": 100,
        "status": "VETO",
    }
    return m


class TestAdjudicateBasicVerdicts:

    def test_all_passing_returns_pass(self):
        inp = _make_input(measurements=_all_passing_measurements(), coverage=100.0)
        result = adjudicate_laws(inp)
        assert result.verdict == "PASS"

    def test_empty_measurements_returns_veto(self):
        """No measurements → all laws default to FAIL → VETO."""
        inp = _make_input(measurements={}, coverage=100.0)
        result = adjudicate_laws(inp)
        assert result.verdict == "VETO"

    def test_low_coverage_returns_indeterminate(self):
        inp = _make_input(measurements=_all_passing_measurements(), coverage=5.0)
        result = adjudicate_laws(inp)
        assert result.verdict == "INDETERMINATE"

    def test_single_deterministic_veto_returns_veto(self):
        # Find a deterministic law to veto
        det_law = None
        for lid, cls in LAW_METHOD_CLASSIFICATIONS.items():
            if cls == "deterministic" and lid in LAW_CANON:
                det_law = lid
                break
        assert det_law is not None, "No deterministic law found"

        inp = _make_input(
            measurements=_single_veto_measurements(det_law), coverage=100.0
        )
        result = adjudicate_laws(inp)
        assert result.verdict == "VETO"
        assert any(det_law in f for f in result.failing_deterministic)


class TestAdjudicateModalityInteraction:

    def test_advisory_veto_downgraded(self):
        """LAW-100 VETO on experimental x_ray → FAIL (Advisory), not VETO verdict."""
        m = _all_passing_measurements()
        m["LAW-100"] = {
            "observed": 999, "deviation": "999", "sample": 100, "status": "VETO"
        }
        inp = _make_input(
            measurements=m, coverage=100.0,
            is_experimental=True, method_type="x_ray",
        )
        result = adjudicate_laws(inp)

        law100_row = next(r for r in result.law_rows if r["law_id"] == "LAW-100")
        assert law100_row["status"] == "FAIL (Advisory)"
        assert law100_row["method"] == "advisory_experimental"
        # Overall verdict should still be PASS since LAW-100 is advisory
        assert result.verdict == "PASS"

    def test_predicted_veto_not_downgraded(self):
        """LAW-100 VETO on predicted → stays VETO, verdict is VETO."""
        m = _all_passing_measurements()
        m["LAW-100"] = {
            "observed": 999, "deviation": "999", "sample": 100, "status": "VETO"
        }
        inp = _make_input(
            measurements=m, coverage=100.0,
            is_experimental=False, method_type="predicted",
        )
        result = adjudicate_laws(inp)

        law100_row = next(r for r in result.law_rows if r["law_id"] == "LAW-100")
        assert law100_row["status"] == "VETO"
        assert result.verdict == "VETO"


class TestAdjudicateScoring:

    def test_all_passing_det_score_is_100(self):
        inp = _make_input(measurements=_all_passing_measurements(), coverage=100.0)
        result = adjudicate_laws(inp)
        assert result.deterministic_score == 100

    def test_empty_det_score_is_zero(self):
        inp = _make_input(measurements={}, coverage=100.0)
        result = adjudicate_laws(inp)
        assert result.deterministic_score == 0

    def test_low_coverage_zeroes_strategic_math(self):
        inp = _make_input(measurements=_all_passing_measurements(), coverage=5.0)
        result = adjudicate_laws(inp)
        assert result.strategic_math.s6 == 0.0
        assert result.strategic_math.p_score == 0
        assert result.strategic_math.suppression is not None


class TestAdjudicateResultStructure:

    def test_law_rows_count_matches_canon(self):
        inp = _make_input(measurements=_all_passing_measurements(), coverage=100.0)
        result = adjudicate_laws(inp)
        assert len(result.law_rows) == len(LAW_CANON)

    def test_law_row_has_required_keys(self):
        inp = _make_input(measurements=_all_passing_measurements(), coverage=100.0)
        result = adjudicate_laws(inp)
        required = {
            "law_id", "title", "status", "method", "observed",
            "threshold", "operator", "units", "deviation",
            "sample_size", "scope", "principle",
        }
        for row in result.law_rows:
            missing = required - set(row.keys())
            assert not missing, f"{row['law_id']} missing keys: {missing}"

    def test_result_is_dataclass(self):
        inp = _make_input(measurements=_all_passing_measurements(), coverage=100.0)
        result = adjudicate_laws(inp)
        assert isinstance(result, AdjudicationResult)
        assert isinstance(result.strategic_math, StrategicMath)


class TestComputeVerdict:

    def test_pass_with_no_failures(self):
        rows = [{"method": "deterministic", "status": "PASS"}]
        assert _compute_verdict(rows, 100.0) == "PASS"

    def test_veto_with_deterministic_fail(self):
        rows = [{"method": "deterministic", "status": "FAIL"}]
        assert _compute_verdict(rows, 100.0) == "VETO"

    def test_veto_with_deterministic_veto_status(self):
        rows = [{"method": "deterministic", "status": "VETO"}]
        assert _compute_verdict(rows, 100.0) == "VETO"

    def test_indeterminate_on_low_coverage(self):
        rows = [{"method": "deterministic", "status": "PASS"}]
        assert _compute_verdict(rows, 1.0) == "INDETERMINATE"

    def test_advisory_fail_does_not_veto(self):
        rows = [
            {"method": "deterministic", "status": "PASS"},
            {"method": "advisory_experimental", "status": "FAIL"},
        ]
        assert _compute_verdict(rows, 100.0) == "PASS"

    def test_heuristic_fail_does_not_veto(self):
        rows = [
            {"method": "deterministic", "status": "PASS"},
            {"method": "heuristic", "status": "FAIL"},
        ]
        assert _compute_verdict(rows, 100.0) == "PASS"
