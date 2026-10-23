"""P8 regression tests: calibration remains diagnostic-only.

No v25 execution, held-out evaluation, or threshold/cohort change is
authorized under OD-16. The evaluation report must declare itself
diagnostic and not_approved with supports_threshold_validation=False, and
the recovered protocol must not be marked as supporting threshold validation.
These tests pin those gates so calibration output cannot be mistaken for an
authoritative threshold validation.
"""
import os
from pathlib import Path

os.environ.setdefault("TOSCANINI_API_KEYS", "test-key-xyz")
REPO = Path(__file__).resolve().parents[2]


class TestCalibrationDiagnosticOnly:

    def test_evaluate_report_is_diagnostic_not_approved(self):
        from calibration.evaluate import evaluate
        cases = [
            {"structure_id": "1CRN", "label": "acceptable", "verdict": "PASS"},
            {"structure_id": "1ABC", "label": "defective", "verdict": "VETO",
             "blocking_laws": ["LAW-120"]},
        ]
        out = evaluate(cases)
        assert out["report_type"] == "diagnostic"
        assert out["calibration_status"] == "not_approved"
        assert out["supports_threshold_validation"] is False

    def test_protocol_documents_v23_as_diagnostic_only(self):
        proto = (REPO / "backend" / "calibration" / "protocol"
                 / "evaluation_protocol_v1.md").read_text(encoding="utf-8")
        # The recovered protocol must not claim a validated threshold: v23 is
        # diagnostic-only and the 18-degree threshold traces to the
        # uncorrected metric.
        assert "diagnostic" in proto.lower()
        assert "uncorrected" in proto.lower()
