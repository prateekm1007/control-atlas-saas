"""
Tests for API v1.0 Pydantic response schema.
Validates that the schema matches the actual response structure exactly.
"""
import json
import pytest
from pathlib import Path
from pydantic import ValidationError
from tos.schemas.response_v1 import (
    ToscaniniResponse,
    GovernanceFingerprint,
    AdjudicationVerdict,
    PhysicsMeasurement,
    Governance,
    Tier1,
)


# Load the sample response captured from the live API
SAMPLE_PATH = Path("/tmp/toscanini_response_sample.json")


class TestSchemaMatchesLiveResponse:
    """Verify the Pydantic models accept the actual API output."""

    @pytest.fixture
    def sample_data(self):
        if not SAMPLE_PATH.exists():
            pytest.skip("No sample response at /tmp/toscanini_response_sample.json")
        with SAMPLE_PATH.open() as f:
            data = json.load(f)
        # Add back fields that were stripped for readability
        data.setdefault("witness_reports", None)
        data.setdefault("pdf_b64", "")
        data.setdefault("pdb_b64", "")
        return data

    def test_full_response_validates(self, sample_data):
        """The complete live response must parse without error."""
        resp = ToscaniniResponse(**sample_data)
        assert resp.verdict.binary == "PASS"

    def test_governance_fingerprint_validates(self, sample_data):
        fp_data = sample_data["governance"]["governance_fingerprint"]
        fp = GovernanceFingerprint(**fp_data)
        assert fp.canon_hash == "6a9cd4b4349b81de"
        assert len(fp.matrix_hash) == 64
        assert fp.policy_ref == "PIL-CAL-03"

    def test_all_15_laws_validate(self, sample_data):
        laws = sample_data["tier1"]["laws"]
        assert len(laws) == 15
        for law_data in laws:
            m = PhysicsMeasurement(**law_data)
            assert m.law_id.startswith("LAW-")

    def test_verdict_enum_constraint(self, sample_data):
        v = AdjudicationVerdict(**sample_data["verdict"])
        assert v.binary in ("PASS", "VETO", "INDETERMINATE")

    def test_scores_are_bounded(self, sample_data):
        v = AdjudicationVerdict(**sample_data["verdict"])
        assert 0 <= v.deterministic_score <= 100
        assert 0 <= v.advisory_score <= 100
        assert 0 <= v.coverage_pct <= 100


class TestSchemaRejectsBadData:
    """Verify the schema catches type errors."""

    def test_rejects_invalid_verdict(self):
        with pytest.raises(ValidationError):
            AdjudicationVerdict(
                binary="MAYBE",  # Not in Literal
                deterministic_score=50, advisory_score=50,
                physical_score=50, confidence_score=90.0,
                det_passed=5, det_total=10, heur_passed=3, heur_total=3,
                coverage_pct=100.0, suppression_reason=None,
            )

    def test_rejects_negative_score(self):
        with pytest.raises(ValidationError):
            AdjudicationVerdict(
                binary="PASS",
                deterministic_score=-1,  # ge=0 violated
                advisory_score=50,
                physical_score=50, confidence_score=90.0,
                det_passed=5, det_total=10, heur_passed=3, heur_total=3,
                coverage_pct=100.0, suppression_reason=None,
            )

    def test_rejects_missing_canon_hash(self):
        with pytest.raises(ValidationError):
            GovernanceFingerprint(
                # canon_hash missing
                matrix_hash="abc123",
                matrix_schema_version="1.0.0",
                policy_ref="PIL-CAL-03",
            )

    def test_rejects_missing_law_id(self):
        with pytest.raises(ValidationError):
            PhysicsMeasurement(
                # law_id missing
                title="Test", status="PASS", method="deterministic",
                observed=1.0, threshold=5.0, operator="<=",
                units="%", deviation="0.0", sample_size=100,
                scope="core residues",
            )
