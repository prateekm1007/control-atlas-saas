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
    BLOCKING_STATUSES,
    INDETERMINATE_CAUSE_STATUSES,
    STATUS_NOT_MEASURED,
    CONTAINED_LAWS,
    LAW_CONTAINMENT,
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

    WP-SCI-1 additions:
        contained_law_failures: law IDs whose measured failure could not
            independently veto because the law is contained (see
            station_sop.LAW_CONTAINMENT); the verdict is INDETERMINATE
            with explicit provenance when this is the only blocking
            evidence.
        indeterminate_causes: law IDs whose status was NOT_MEASURED or
            MEASUREMENT_FAILED on a deterministic law — missing mandatory
            measurements prevent an authoritative PASS.
    """
    law_rows: list = field(default_factory=list)
    failing_deterministic: list = field(default_factory=list)
    verdict: str = "INDETERMINATE"
    deterministic_score: int = 0
    advisory_score: int = 0
    det_passed: int = 0
    heur_passed: int = 0
    det_score_normalized: int = 0
    det_passed_normalized: int = 0
    contained_law_failures: list = field(default_factory=list)
    indeterminate_causes: list = field(default_factory=list)
    strategic_math: StrategicMath = field(
        default_factory=lambda: StrategicMath(0.0, 1.0, 0.0, 0, None)
    )


# ═══════════════════════════════════════════════════════════════════════════════
# VERDICT LOGIC (moved verbatim from main.py _compute_verdict)
# ═══════════════════════════════════════════════════════════════════════════════

def _compute_verdict(res_t1: list[dict], coverage: float) -> tuple[str, list[str], list[str]]:
    """
    Determine the verdict from classified law rows and coverage.

    WP-SCI-1 semantics (directive §1–§2, Constitution §4):
    - Coverage below LAW-105 threshold → INDETERMINATE
    - Any deterministic law measured FAIL or VETO → VETO, EXCEPT a
      contained law cannot independently issue the veto: if contained
      failures exist but no uncontained deterministic law failed, the
      verdict is INDETERMINATE (honest non-success, no physical claim).
    - Any deterministic law NOT_MEASURED or MEASUREMENT_FAILED blocks an
      authoritative PASS → INDETERMINATE (missing evidence is neither
      favorable evidence nor a fabricated violation).
    - NOT_APPLICABLE laws neither pass nor veto.
    - Otherwise → PASS

    Returns:
        (verdict, contained_law_failures, indeterminate_causes)
    """
    cov_threshold = float(LAW_CANON["LAW-105"]["threshold"])
    if coverage < cov_threshold:
        return "INDETERMINATE", [], []

    contained_failures: list[str] = []
    indeterminate_causes: list[str] = []
    uncontained_veto = False

    for r in res_t1:
        if r["method"] != "deterministic":
            continue
        status = r["status"]
        if status in ("FAIL", "VETO"):
            if r.get("law_id", "") in CONTAINED_LAWS:
                contained_failures.append(r["law_id"])
            else:
                uncontained_veto = True
        elif status in INDETERMINATE_CAUSE_STATUSES:
            indeterminate_causes.append(r["law_id"])

    if uncontained_veto:
        return "VETO", contained_failures, indeterminate_causes
    if contained_failures or indeterminate_causes:
        # Missing mandatory measurements and contained-law failures are
        # explicit non-success states: no authoritative PASS.
        return "INDETERMINATE", contained_failures, indeterminate_causes

    return "PASS", [], []


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
        # WP-SCI-1: a law entirely absent from the measurements is NOT
        # MEASURED, not "measured and failed". (run_full_audit currently
        # emits every law, so this default is a fail-closed guard for
        # future callers rather than a live path.)
        m = inp.raw_measurements.get(
            lid,
            {"observed": None, "deviation": "N/A", "sample": 0,
             "status": STATUS_NOT_MEASURED,
             "reason": "no measurement row emitted by the engine"},
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

        # PIL-CAL-03: Advisory laws cannot VETO — normalize status.
        # P0 item 8: a present-but-malformed row (no status key) is
        # NOT_MEASURED, not "measured and failed" — mirroring the
        # absent-law default above. Fail-closed behavior is preserved:
        # NOT_MEASURED on a deterministic law blocks an authoritative
        # PASS via the indeterminate-causes path.
        raw_status = m.get("status", STATUS_NOT_MEASURED)
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
        # Phase A.5: Propagate residue-level diagnostics if present
        if "granularity" in m:
            row["granularity"] = m["granularity"]
        if "failing_residues" in m:
            row["failing_residues"] = m["failing_residues"]
        if "failing_residue_pairs" in m:
            row["failing_residue_pairs"] = m["failing_residue_pairs"]
        # WP-SCI-1: containment provenance and explicit non-success reason
        if lid in CONTAINED_LAWS:
            row["contained"] = True
            row["containment_reason"] = LAW_CONTAINMENT[lid].get("reason", "contained")
        if "reason" in m:
            row["reason"] = m["reason"]
        res_t1.append(row)
        if raw_status in ("FAIL", "VETO") and row["method"] == "deterministic":
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

    # WP-SCI-1: verdict carries containment and missing-measurement provenance.
    verdict, contained_failures, indeterminate_causes = _compute_verdict(res_t1, inp.coverage)

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
        det_score_normalized=det_score_normalized,
        det_passed_normalized=det_passed_normalized,
        contained_law_failures=contained_failures,
        indeterminate_causes=indeterminate_causes,
        strategic_math=StrategicMath(
            s6=s6_final, w_arch=w_arch, m_s8=m_s8,
            p_score=p_score, suppression=suppression,
        ),
    )
