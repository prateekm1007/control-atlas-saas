"""
Adjudicator: Pure Physics Verdict Engine
-----------------------------------------
Stages 3 (Adjudication) and 4 (Scoring) extracted from main.py.

This module has ZERO knowledge of:
- FastAPI, HTTP, file I/O
- PDF generation (narrative layer removed, H-1)
- PDB bytes or Structure objects

It takes measurements in, returns verdicts out.

All thresholds pulled from station_sop.py (Contract Sovereignty).
All modality logic delegated to modality_matrix.py.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

from tos.governance.station_sop import (
    LAW_CANON,
    DETERMINISTIC_COUNT,
    HEURISTIC_COUNT,
    LAW_METHOD_CLASSIFICATIONS,
    RETIRED_LAWS,
    ARCHITECTURE_WEIGHTS,
)
from tos.governance.modality_matrix import resolve_method


MEASURED_STATUSES = frozenset({"PASS", "FAIL", "VETO"})
MISSING_STATUSES = frozenset({"NOT_MEASURED", "MEASUREMENT_FAILED"})
VALID_STATUSES = MEASURED_STATUSES | MISSING_STATUSES | {"NOT_APPLICABLE"}


# ═══════════════════════════════════════════════════════════════════════════════
# DATA CONTRACTS
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class AdjudicationInput:
    """
    Everything the adjudicator needs. No Structure object, no HTTP context.

    Attributes:
        raw_measurements: lid → {observed, deviation, sample, status} from Tier1
        coverage:         Atom coverage percentage from run_full_audit
        is_experimental:  True if structure came from experimental source
        method_type:      Acquisition method: "x_ray", "cryo_em", "nmr", "predicted"
        architecture:     t3_category passed by caller (e.g. "NONE", "MONOMER")
    """
    raw_measurements: dict
    coverage: float
    is_experimental: bool
    method_type: str
    architecture: str


@dataclass
class StrategicMath:
    """PIL-MTH-12 scoring breakdown."""
    s6: float
    w_arch: float
    m_s8: float
    p_score: int
    suppression: Optional[str]


@dataclass
class AdjudicationResult:
    """
    Complete adjudication output. Pure data, no side effects.

    Every field here was previously computed inline in _run_physics_sync.
    The field names and value semantics are preserved exactly.
    """
    law_rows: list = field(default_factory=list)
    failing_deterministic: list = field(default_factory=list)
    verdict: str = "INDETERMINATE"
    verdict_reason: Optional[str] = None
    deterministic_score: int = 0
    advisory_score: int = 0
    det_passed: int = 0
    heur_passed: int = 0
    det_score_normalized: int = 0
    det_passed_normalized: int = 0
    strategic_math: StrategicMath = field(
        default_factory=lambda: StrategicMath(0.0, 1.0, 0.0, 0, None)
    )


# ═══════════════════════════════════════════════════════════════════════════════
# VERDICT LOGIC (moved verbatim from main.py _compute_verdict)
# ═══════════════════════════════════════════════════════════════════════════════

def _compute_verdict(res_t1: list[dict], coverage: float) -> str:
    """
    Determine binary verdict from classified law rows and coverage.

    Logic preserved exactly from main.py:
    - Coverage below LAW-105 threshold → INDETERMINATE
    - Any deterministic law with FAIL or VETO → VETO
    - Otherwise → PASS
    """
    cov_threshold = float(LAW_CANON["LAW-105"]["threshold"])
    if coverage < cov_threshold:
        return "INDETERMINATE"

    for r in res_t1:
        if r["method"] == "deterministic" and r["status"] in ("FAIL", "VETO"):
            return "VETO"

    for r in res_t1:
        if r["method"] == "deterministic" and r["status"] in ("NOT_MEASURED", "MEASUREMENT_FAILED"):
            return "INDETERMINATE"

    return "PASS"


def _normalize_measurement(m: dict) -> tuple[str, object, int, Optional[str]]:
    """Fail closed on unknown or internally inconsistent measurement rows."""
    raw_status = m.get("status", "NOT_MEASURED")
    status = str(raw_status).strip().upper() if isinstance(raw_status, str) else "INVALID"
    observed = m.get("observed")
    try:
        sample = int(m.get("sample", 0))
    except (TypeError, ValueError):
        return "MEASUREMENT_FAILED", None, 0, "sample size is not an integer"
    if sample < 0:
        return "MEASUREMENT_FAILED", None, 0, "sample size is negative"

    if status not in VALID_STATUSES:
        return "MEASUREMENT_FAILED", None, sample, f"unknown measurement status {raw_status!r}"

    if status in MEASURED_STATUSES and (observed is None or sample == 0):
        return (
            "MEASUREMENT_FAILED",
            None,
            sample,
            "measured status claimed without a numeric observation and positive sample",
        )

    if status in MEASURED_STATUSES:
        try:
            finite_observation = math.isfinite(float(observed))
        except (TypeError, ValueError):
            finite_observation = False
        if not finite_observation:
            return (
                "MEASUREMENT_FAILED",
                None,
                sample,
                "measured status carries a non-finite or non-numeric observation",
            )

    if status == "NOT_APPLICABLE":
        if observed is not None or sample != 0:
            return (
                "MEASUREMENT_FAILED",
                None,
                sample,
                "NOT_APPLICABLE must not carry an observation or nonzero sample",
            )
        return status, None, 0, None

    if status in MISSING_STATUSES and observed is not None:
        return (
            "MEASUREMENT_FAILED",
            None,
            sample,
            "missing-evidence status carried a numeric observation",
        )

    return status, observed, sample, None


# ═══════════════════════════════════════════════════════════════════════════════
# CORE ADJUDICATION (Stages 3 + 4)
# ═══════════════════════════════════════════════════════════════════════════════

def adjudicate_laws(inp: AdjudicationInput) -> AdjudicationResult:
    """
    Pure adjudication: measurements in, verdict out.

    This is the exact logic from main.py lines ~68–130, moved verbatim.
    The loop order, rounding, and conditional structure are preserved.

    Args:
        inp: AdjudicationInput with raw measurements and context.

    Returns:
        AdjudicationResult with classified law rows, scores, and verdict.
    """
    res_t1 = []
    failing_det = []

    # ── Stage 3: Law Classification Loop ─────────────────────────────────
    for lid in sorted(LAW_CANON.keys()):
        m = inp.raw_measurements.get(
            lid,
            {"observed": None, "deviation": None, "sample": 0, "status": "NOT_MEASURED"},
        )
        method = LAW_METHOD_CLASSIFICATIONS.get(lid, "unknown")

        # OD-03: LAW-155 is retained for traceability but retired from
        # enforced semantics. No modality rule may reactivate it.
        if lid in RETIRED_LAWS:
            method = "retired"
        else:
            # PIL-CAL-03: Modality-Aware Enforcement Matrix (v22.5.3)
            # Delegated to modality_matrix.json — see modality_matrix.py for loader.
            method = resolve_method(
                law_id=lid,
                is_experimental=inp.is_experimental,
                method_type=inp.method_type,
                default_method=method,
            )

        # PIL-CAL-03: Advisory/heuristic laws cannot issue a VETO.
        raw_status = m.get("status", "NOT_MEASURED")
        normalized, observed, sample, normalization_error = _normalize_measurement(m)
        if method == "retired":
            display_status = "RETIRED"
        elif method in ("advisory_experimental", "heuristic") and normalized == "VETO":
            display_status = "FAIL (Advisory)"
        else:
            display_status = normalized

        row = {
            "law_id": lid, "title": LAW_CANON[lid]["title"], "status": display_status,
            "raw_status": raw_status, "method": method,
            "enforced": method != "retired",
            "observed": observed, "threshold": LAW_CANON[lid]["threshold"],
            "operator": LAW_CANON[lid]["operator"], "units": LAW_CANON[lid]["unit"],
            "deviation": m.get("deviation"), "sample_size": sample,
            "scope": LAW_CANON[lid]["scope"], "principle": LAW_CANON[lid].get("principle", "N/A")
        }
        # Phase A.5: Propagate residue-level diagnostics if present
        for key in ("granularity", "failing_residues", "failing_residue_pairs",
                    "reason", "plddt_advisory"):
            if key in m:
                row[key] = m[key]
        if normalization_error is not None:
            prior_reason = row.get("reason")
            row["reason"] = (
                f"{prior_reason}; {normalization_error}" if prior_reason else normalization_error
            )
        res_t1.append(row)
        if display_status in ("FAIL", "VETO") and row["method"] == "deterministic":
            failing_det.append(f"{lid}: {LAW_CANON[lid]['title']}")

    # ── Stage 4: Scoring ─────────────────────────────────────────────────
    det_passed = sum(1 for r in res_t1 if r["status"] == "PASS" and r["method"] == "deterministic")
    det_score = int((det_passed / max(DETERMINISTIC_COUNT, 1)) * 100)
    # WP-01: normalized score includes advisory_experimental laws for cross-source comparability
    det_passed_normalized = sum(1 for r in res_t1 if r["status"] == "PASS"
                                and r["method"] in ("deterministic", "advisory_experimental"))
    det_score_normalized = int((det_passed_normalized / max(DETERMINISTIC_COUNT, 1)) * 100)
    heur_passed = sum(1 for r in res_t1 if r["status"] == "PASS" and r["method"] == "heuristic")
    adv_score = int((heur_passed / max(HEURISTIC_COUNT, 1)) * 100)

    verdict = _compute_verdict(res_t1, inp.coverage)

    # PIL-MTH-12: Strategic Success Math
    raw_s6 = round(det_score / 100.0, 4)
    w_arch = ARCHITECTURE_WEIGHTS.get(inp.architecture, 1.0)
    m_s8 = 0.0  # NKG Penalty

    # Coverage Gate Logic
    cov_threshold = float(LAW_CANON["LAW-105"]["threshold"])
    if inp.coverage < cov_threshold:
        s6_final = 0.0
        p_score = 0
        suppression = "Coverage gate enforces zeroing of prioritization index."
    else:
        s6_final = raw_s6
        p_score = int(s6_final * w_arch * (1.0 - m_s8) * 100)
        suppression = None

    if inp.coverage < cov_threshold:
        verdict_reason = (
            f"LAW-105 residue coverage {round(inp.coverage, 2)}% is below "
            f"the {cov_threshold}% minimum"
        )
    elif verdict == "VETO":
        verdict_reason = "Deterministic violation: " + "; ".join(failing_det)
    elif verdict == "INDETERMINATE":
        unavailable = [
            r["law_id"] for r in res_t1
            if r["method"] == "deterministic"
            and r["status"] in ("NOT_MEASURED", "MEASUREMENT_FAILED")
        ]
        verdict_reason = "Mandatory deterministic measurement unavailable: " + ", ".join(unavailable)
    else:
        verdict_reason = None

    return AdjudicationResult(
        law_rows=res_t1,
        failing_deterministic=failing_det,
        verdict=verdict,
        verdict_reason=verdict_reason,
        deterministic_score=det_score,
        advisory_score=adv_score,
        det_passed=det_passed,
        heur_passed=heur_passed,
        det_score_normalized=det_score_normalized,
        det_passed_normalized=det_passed_normalized,
        strategic_math=StrategicMath(
            s6=s6_final, w_arch=w_arch, m_s8=m_s8,
            p_score=p_score, suppression=suppression,
        ),
    )
