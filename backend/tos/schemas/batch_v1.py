"""
Batch API v1.0 Response Schema
-------------------------------
Pydantic v2 models for the /v1/batch endpoint.
"""
from __future__ import annotations
from typing import Optional
from pydantic import BaseModel, Field
from .response_v1 import ToscaniniResponse


class BatchStructureResult(BaseModel):
    """Result for a single structure within a batch."""
    filename: str
    candidate_id: str
    success: bool
    response: Optional[ToscaniniResponse] = None
    error: Optional[str] = None


class FailingLawCount(BaseModel):
    law_id: str
    count: int


class BatchSummary(BaseModel):
    """Aggregate statistics across the batch."""
    total: int = Field(..., ge=0)
    passed: int = Field(..., ge=0)
    vetoed: int = Field(..., ge=0)
    indeterminate: int = Field(..., ge=0)
    errors: int = Field(..., ge=0)
    mean_deterministic_score: float = Field(..., ge=0, le=100)
    common_failing_laws: list[FailingLawCount] = Field(default_factory=list)


class BatchResponse(BaseModel):
    """Complete /v1/batch response."""
    results: list[BatchStructureResult]
    summary: BatchSummary
