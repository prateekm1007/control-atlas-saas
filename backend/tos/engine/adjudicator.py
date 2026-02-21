"""
Adjudicator: Pure Physics Verdict Engine
-----------------------------------------
Stages 3 (Adjudication) and 4 (Scoring) extracted from main.py.

This module has ZERO knowledge of:
- FastAPI, HTTP, file I/O
- PDF generation, Gemini narrative
- PDB bytes or Structure objects

It takes measurements in, returns verdicts out.

All thresholds pulled from station_sop.py (Contract Sovereignty).
All modality logic delegated to modality_matrix.py.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from tos.governance.station_sop import (
    LAW_CANON,
    DETERMINISTIC_COUNT,
    HEURISTIC_COUNT,
    LAW_METHOD_CLASSIFICATIONS,
    ARCHITECTURE_WEIGHTS,
)
from tos.governance.modality_matrix import resolve_method


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
    deterministic_score: int = 0
    advisory_score: int = 0
    det_passed: int = 0
    heur_passed: int = 0
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

    return "PASS"


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
            lid, {"observed": 0, "deviation": "0.0", "sample": 0, "status": "FAIL"}
        )
        method = LAW_METHOD_CLASSIFICATIONS.get(lid, "unknown")

        # PIL-CAL-03: Modality-Aware Enforcement Matrix (v22.5.3)
        # Delegated to modality_matrix.json — see modality_matrix.py for loader.
        method = resolve_method(
            law_id=lid,
            is_experimental=inp.is_experimental,
            method_type=inp.method_type,
            default_method=method,
        )

        # PIL-CAL-03: Advisory laws cannot VETO — normalize status
        raw_status = m.get("status", "FAIL")
        if method == "advisory_experimental" and raw_status == "VETO":
            raw_status = "FAIL (Advisory)"

        row = {
            "law_id": lid, "title": LAW_CANON[lid]["title"], "status": raw_status,
            "method": method,
            "observed": m.get("observed", 0), "threshold": LAW_CANON[lid]["threshold"],
            "operator": LAW_CANON[lid]["operator"], "units": LAW_CANON[lid]["unit"],
            "deviation": m.get("deviation", "0.0"), "sample_size": m.get("sample", 0),
            "scope": LAW_CANON[lid]["scope"], "principle": LAW_CANON[lid].get("principle", "N/A")
        }
        res_t1.append(row)
        if m.get("status", "FAIL") in ("FAIL", "VETO") and row["method"] == "deterministic":
            failing_det.append(f"{lid}: {LAW_CANON[lid]['title']}")

    # ── Stage 4: Scoring ─────────────────────────────────────────────────
    det_passed = sum(1 for r in res_t1 if r["status"] == "PASS" and r["method"] == "deterministic")
    det_score = int((det_passed / max(DETERMINISTIC_COUNT, 1)) * 100)
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

    return AdjudicationResult(
        law_rows=res_t1,
        failing_deterministic=failing_det,
        verdict=verdict,
        deterministic_score=det_score,
        advisory_score=adv_score,
        det_passed=det_passed,
        heur_passed=heur_passed,
        strategic_math=StrategicMath(
            s6=s6_final, w_arch=w_arch, m_s8=m_s8,
            p_score=p_score, suppression=suppression,
        ),
    )
