"""
Toscanini API v1.0 Response Schema
-----------------------------------
Pydantic v2 models enforcing the exact contract of the /ingest response.
"""
from __future__ import annotations
from typing import Literal, Optional
from pydantic import BaseModel, Field, model_validator

MeasurementStatus = Literal[
    "PASS",
    "FAIL",
    "VETO",
    "FAIL (Advisory)",
    "NOT_MEASURED",
    "MEASUREMENT_FAILED",
    "NOT_APPLICABLE",
    "RETIRED",
]


class GovernanceFingerprint(BaseModel):
    canon_hash: str
    matrix_hash: str
    matrix_schema_version: str
    policy_ref: str


class Governance(BaseModel):
    audit_id: str
    station_version: str
    timestamp_utc: str
    governance_fingerprint: GovernanceFingerprint


class AdjudicationVerdict(BaseModel):
    binary: Literal["PASS", "VETO", "INDETERMINATE"]
    reason: Optional[str] = None
    coverage_basis: Optional[str] = None
    deterministic_score: int = Field(..., ge=0, le=100)
    advisory_score: int = Field(..., ge=0, le=100)
    physical_score: int = Field(..., ge=0, le=100)
    confidence_score: Optional[float] = Field(None, ge=0, le=100)
    det_passed: int = Field(..., ge=0)
    det_score_normalized: Optional[int] = Field(None, ge=0, le=100)
    det_passed_normalized: Optional[int] = Field(None, ge=0)
    det_total: int = Field(..., ge=0)
    heur_passed: int = Field(..., ge=0)
    heur_total: int = Field(..., ge=0)
    coverage_pct: float = Field(..., ge=0, le=100)
    suppression_reason: Optional[str] = None
    primary_cause: Optional[str] = None
    recommended_action: Optional[str] = None
    veto_law: Optional[str] = None


class PhysicsMeasurement(BaseModel):
    law_id: str
    title: str
    status: MeasurementStatus
    method: str
    observed: Optional[float] = None
    threshold: float
    operator: str
    units: str
    deviation: Optional[str] = None
    sample_size: int = Field(..., ge=0)
    scope: str
    principle: str = "N/A"

    @model_validator(mode="after")
    def measurement_evidence_is_consistent(self):
        measured = {"PASS", "FAIL", "VETO", "FAIL (Advisory)"}
        missing = {"NOT_MEASURED", "MEASUREMENT_FAILED"}
        if self.status in measured and (self.observed is None or self.sample_size == 0):
            raise ValueError(
                "measured status requires a numeric observation and positive sample size"
            )
        if self.status in missing and self.observed is not None:
            raise ValueError("missing-evidence status must not carry a numeric observation")
        if self.status == "NOT_APPLICABLE" and (self.observed is not None or self.sample_size != 0):
            raise ValueError("NOT_APPLICABLE must have no observation and zero sample size")
        return self


class Tier1(BaseModel):
    laws: list[PhysicsMeasurement]


class Tier3(BaseModel):
    probability: int = Field(..., ge=0, le=100)


class Provenance(BaseModel):
    source: str
    hash: str
    byte_count: int = Field(..., ge=0)


class Characterization(BaseModel):
    total_atoms: int = Field(..., ge=0)
    total_residues: int = Field(..., ge=0)
    source_type: str
    resolution: Optional[float] = None
    method: str = "Standard"


class ConfidenceMeta(BaseModel):
    schema_version: Literal["confidence-meta/1.0.0"] = "confidence-meta/1.0.0"
    source_type: str
    provenance_method: str
    mean: Optional[float] = None
    data_available: bool
    observation_count: int = Field(0, ge=0)

    @model_validator(mode="after")
    def availability_matches_mean(self):
        if self.data_available and self.mean is None:
            raise ValueError("available confidence metadata requires a mean")
        if not self.data_available and self.mean is not None:
            raise ValueError("unavailable confidence metadata must not carry a mean")
        return self


class StrategicMath(BaseModel):
    s6: float
    w_arch: float
    m_s8: float
    architecture: str


class ToscaniniResponse(BaseModel):
    verdict: AdjudicationVerdict
    governance: Governance
    provenance: Provenance
    tier1: Tier1
    tier3: Tier3
    characterization: Characterization
    confidence_meta: ConfidenceMeta
    strategic_math: StrategicMath
    witness_reports: Optional[dict] = None
    ai_model_used: Optional[str] = None
    pdf_b64: Optional[str] = ""
    pdb_b64: Optional[str] = ""
