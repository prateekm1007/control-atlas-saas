"""
dashboard/comparison_engine.py
Toscanini Before/After Audit Comparison Engine

Pure utility — no side effects. Takes two audit results and returns
a structured comparison dict.
"""


def compare_audits(baseline, refined):
    """
    Compare two audit results (before refinement vs after).
    
    Returns dict with:
      - coverage_delta (float)
      - det_score_delta (int)
      - verdict_change (str)
      - law_changes (list of dicts)
      - violation_count_delta (int)
      - regressions (list of law_ids)
      - improvements (list of law_ids)
    """
    b_verdict = baseline.get("verdict", {})
    r_verdict = refined.get("verdict", {})
    
    b_laws = {l["law_id"]: l for l in baseline.get("tier1", {}).get("laws", [])}
    r_laws = {l["law_id"]: l for l in refined.get("tier1", {}).get("laws", [])}
    
    all_law_ids = sorted(set(b_laws.keys()) | set(r_laws.keys()))
    
    law_changes = []
    improvements = []
    regressions = []
    
    for lid in all_law_ids:
        b_law = b_laws.get(lid, {})
        r_law = r_laws.get(lid, {})
        
        b_status = b_law.get("status", "MISSING")
        r_status = r_law.get("status", "MISSING")
        
        # Only track deterministic law changes for the summary
        is_det = b_law.get("method") == "deterministic" or r_law.get("method") == "deterministic"
        
        if b_status != r_status:
            change_type = "UNCHANGED"
            if b_status in ("FAIL", "VETO") and r_status == "PASS":
                change_type = "RESOLVED"
                if is_det:
                    improvements.append(lid)
            elif b_status == "PASS" and r_status in ("FAIL", "VETO"):
                change_type = "REGRESSED"
                if is_det:
                    regressions.append(lid)
            elif b_status == "MISSING" or r_status == "MISSING":
                change_type = "STRUCTURAL_CHANGE"
            else:
                change_type = "CHANGED"
            
            law_changes.append({
                "law_id": lid,
                "title": r_law.get("title") or b_law.get("title", ""),
                "method": r_law.get("method") or b_law.get("method", ""),
                "before_status": b_status,
                "after_status": r_status,
                "before_observed": b_law.get("observed", "N/A"),
                "after_observed": r_law.get("observed", "N/A"),
                "threshold": r_law.get("threshold") or b_law.get("threshold", "N/A"),
                "change_type": change_type,
            })
    
    # Coverage delta
    b_cov = b_verdict.get("coverage_pct", 0.0)
    r_cov = r_verdict.get("coverage_pct", 0.0)
    
    # Deterministic score delta
    b_det = b_verdict.get("deterministic_score", 0)
    r_det = r_verdict.get("deterministic_score", 0)
    
    # Verdict change
    b_bin = b_verdict.get("binary", "UNKNOWN")
    r_bin = r_verdict.get("binary", "UNKNOWN")
    verdict_change = f"{b_bin} → {r_bin}"
    
    # Violation counts (deterministic only)
    b_det_total = b_verdict.get("det_total", 12)
    r_det_total = r_verdict.get("det_total", 12)
    b_det_passed = b_verdict.get("det_passed", 0)
    r_det_passed = r_verdict.get("det_passed", 0)
    b_violations = b_det_total - b_det_passed
    r_violations = r_det_total - r_det_passed
    
    return {
        "coverage_delta": round(r_cov - b_cov, 2),
        "det_score_delta": r_det - b_det,
        "verdict_change": verdict_change,
        "verdict_improved": _verdict_rank(r_bin) > _verdict_rank(b_bin),
        "law_changes": law_changes,
        "violation_count_before": b_violations,
        "violation_count_after": r_violations,
        "violation_count_delta": r_violations - b_violations,
        "improvements": improvements,
        "regressions": regressions,
        "baseline_audit_id": baseline.get("governance", {}).get("audit_id", "UNKNOWN"),
        "refined_audit_id": refined.get("governance", {}).get("audit_id", "UNKNOWN"),
    }


def _verdict_rank(verdict):
    """Internal: rank verdicts for comparison. Higher is better."""
    return {"VETO": 0, "INDETERMINATE": 1, "PASS": 2}.get(verdict, -1)
