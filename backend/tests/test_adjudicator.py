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

    def test_empty_measurements_is_not_authoritative(self):
        """WP-SCI-1: no measurements -> every law NOT_MEASURED -> the verdict
        is an honest INDETERMINATE, never a fabricated VETO or PASS."""
        inp = _make_input(measurements={}, coverage=100.0)
        result = adjudicate_laws(inp)
        assert result.verdict == "INDETERMINATE"
        # every deterministic law is an explicit missing-measurement cause;
        # heuristic laws cannot block a verdict either way.
        assert len(result.indeterminate_causes) == DETERMINISTIC_COUNT
        assert result.contained_law_failures == []

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
        rows = [{"law_id": "LAW-105", "method": "deterministic", "status": "PASS"}]
        assert _compute_verdict(rows, 100.0) == ("PASS", [], [])

    def test_veto_with_deterministic_fail(self):
        rows = [{"law_id": "LAW-110", "method": "deterministic", "status": "FAIL"}]
        assert _compute_verdict(rows, 100.0) == ("VETO", [], [])

    def test_veto_with_deterministic_veto_status(self):
        rows = [{"law_id": "LAW-160", "method": "deterministic", "status": "VETO"}]
        assert _compute_verdict(rows, 100.0) == ("VETO", [], [])

    def test_indeterminate_on_low_coverage(self):
        rows = [{"law_id": "LAW-105", "method": "deterministic", "status": "PASS"}]
        assert _compute_verdict(rows, 1.0) == ("INDETERMINATE", [], [])

    def test_advisory_fail_does_not_veto(self):
        rows = [
            {"law_id": "LAW-105", "method": "deterministic", "status": "PASS"},
            {"law_id": "LAW-100", "method": "advisory_experimental", "status": "FAIL"},
        ]
        assert _compute_verdict(rows, 100.0) == ("PASS", [], [])

    def test_heuristic_fail_does_not_veto(self):
        rows = [
            {"law_id": "LAW-105", "method": "deterministic", "status": "PASS"},
            {"law_id": "LAW-182", "method": "heuristic", "status": "FAIL"},
        ]
        assert _compute_verdict(rows, 100.0) == ("PASS", [], [])

    @pytest.mark.parametrize("law_id", ["LAW-155", "LAW-182", "LAW-200"])
    def test_heuristic_veto_is_reported_as_advisory(self, law_id):
        measurements = _all_passing_measurements()
        measurements[law_id] = {
            "observed": 0.0, "deviation": "999.0", "sample": 10,
            "status": "VETO",
        }
        result = adjudicate_laws(_make_input(
            measurements=measurements, coverage=100.0))
        row = next(r for r in result.law_rows if r["law_id"] == law_id)
        assert row["method"] == "heuristic"
        assert row["status"] == "FAIL (Advisory)"
        assert result.verdict == "PASS"

    # ── WP-SCI-1: missing-measurement and containment semantics ──────

    def test_not_measured_blocks_authoritative_pass(self):
        """A deterministic NOT_MEASURED row must prevent a PASS."""
        rows = [
            {"law_id": "LAW-105", "method": "deterministic", "status": "PASS"},
            {"law_id": "LAW-120", "method": "deterministic", "status": "NOT_MEASURED"},
        ]
        verdict, contained, causes = _compute_verdict(rows, 100.0)
        assert verdict == "INDETERMINATE"
        assert causes == ["LAW-120"]
        assert contained == []

    def test_measurement_failed_blocks_authoritative_pass(self):
        rows = [
            {"law_id": "LAW-105", "method": "deterministic", "status": "PASS"},
            {"law_id": "LAW-130", "method": "deterministic", "status": "MEASUREMENT_FAILED"},
        ]
        verdict, _, causes = _compute_verdict(rows, 100.0)
        assert verdict == "INDETERMINATE"
        assert causes == ["LAW-130"]

    def test_not_applicable_neither_passes_nor_vetoes(self):
        """LAW-195 with no disulfide pair is recorded, not counted as a
        pass, and does not block the verdict by itself."""
        rows = [
            {"law_id": "LAW-105", "method": "deterministic", "status": "PASS"},
            {"law_id": "LAW-195", "method": "deterministic", "status": "NOT_APPLICABLE"},
        ]
        verdict, contained, causes = _compute_verdict(rows, 100.0)
        assert verdict == "PASS"
        assert contained == []
        assert causes == []
