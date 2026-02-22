"""
Toscanini API v1.0 Response Schema
-----------------------------------
Pydantic v2 models enforcing the exact contract of the /ingest response.
"""
from __future__ import annotations
from typing import Literal, Optional
from pydantic import BaseModel, Field


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
    det_total: int = Field(..., ge=0)
    heur_passed: int = Field(..., ge=0)
    heur_total: int = Field(..., ge=0)
    coverage_pct: float = Field(..., ge=0, le=100)
    suppression_reason: Optional[str] = None


class PhysicsMeasurement(BaseModel):
    law_id: str
    title: str
    status: str
    method: str
    observed: float
    threshold: float
    operator: str
    units: str
    deviation: str
    sample_size: int = Field(..., ge=0)
    scope: str
    principle: str = "N/A"


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
