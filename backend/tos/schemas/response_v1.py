"""
Toscanini API v1.0 Response Schema
-----------------------------------
Pydantic v2 models enforcing the exact contract of the /ingest response.
"""
from __future__ import annotations
from typing import Literal, Optional
from pydantic import BaseModel, Field, model_validator


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
    deterministic_score: int = Field(..., ge=0, le=100)
    advisory_score: int = Field(..., ge=0, le=100)
    physical_score: int = Field(..., ge=0, le=100)
    confidence_score: float
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
    status: str
    method: str
    # Round-24 (the auditor-approved NOT_MEASURED state): a measurement that
    # could not be taken is represented truthfully - `measured` is False, the
    # observation is null (still REQUIRED: a measurement always carries the
    # field) and the deviation is null. An ordinary measurement keeps its
    # number.
    measured: bool = True
    observed: Optional[float]
    threshold: float
    operator: str
    units: str
    deviation: Optional[str]
    sample_size: int = Field(..., ge=0)
    scope: str
    principle: str = "N/A"

    @model_validator(mode="after")
    def _measurement_is_self_consistent(self) -> "PhysicsMeasurement":
        """A measurement row must agree with itself.

        Round-24 (the auditor-approved NOT_MEASURED contract): a law that
        could not be evaluated is reported as NOT MEASURED with a null
        observation, a null deviation, a zero sample and measured=False. A row
        that claims a number while saying it was not measured, or that claims
        to have been measured while carrying no number, is contradictory and
        is rejected here: the contract cannot be violated by accident in a
        response.
        """
        if self.status == "NOT_MEASURED":
            if self.measured:
                raise ValueError(
                    "NOT_MEASURED requires measured=false, got measured=true")
            if self.observed is not None:
                raise ValueError(
                    "NOT_MEASURED requires a null observation, got %r" % (self.observed,))
            if self.deviation is not None:
                raise ValueError(
                    "NOT_MEASURED requires a null deviation, got %r" % (self.deviation,))
            if self.sample_size != 0:
                raise ValueError(
                    "NOT_MEASURED requires sample_size=0, got %r" % (self.sample_size,))
        else:
            if not self.measured:
                raise ValueError(
                    "an unmeasured row must carry status=NOT_MEASURED, got %r" % (self.status,))
            if self.observed is None:
                raise ValueError(
                    "a measured row requires a numeric observation, got None")
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
    source_type: str
    provenance_method: str
    mean: float
    data_available: bool


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
