"""
WP-SCI-1 Evaluation Tooling (directive §5)
===========================================

Computes an independently measured confusion matrix with uncertainty
intervals for a labeled calibration/evaluation run, plus per-law
contribution and missing-measurement statistics. Labels come from the
evaluation protocol (independently defined), NEVER from Toscanini's own
verdicts.

Input format (JSON):
{
  "policy_binding": {"law_semantic_hash": "...", "code_sha": "..."},
  "cases": [
    {
      "structure_id": "1ABC",
      "label": "acceptable" | "defective",       # protocol-defined truth
      "strata": {"modality": "x_ray", "resolution_band": "medium"},
      "verdict": "PASS" | "VETO" | "INDETERMINATE",   # engine output
      "blocking_laws": ["LAW-110"],                   # measured veto causes
      "indeterminate_causes": ["LAW-120"],            # missing measurements
      "not_applicable": ["LAW-195"],
      "not_measured": ["LAW-120"]
    }, ...
  ]
}

Semantics (fixed by protocol v1 — see
backend/calibration/protocol/evaluation_protocol_v1.md):

- positive class = "defective" (a true positive is a VETO on a defective
  structure — the tool's job is detecting defects)
- VETO              -> predicted defective
- PASS              -> predicted acceptable
- INDETERMINATE     -> neither; counted in the indeterminate rate and
  reported separately. It never inflates sensitivity or specificity.
  (Under current policy a contained-law failure alone yields
  INDETERMINATE, so a structure can be defective and INDETERMINATE.)

Statistics: Wilson score intervals (95%) for every rate. No rate is
reported without its interval and denominator.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

DEFAULT_ALPHA = 0.05


def wilson_interval(successes: int, total: int, alpha: float = DEFAULT_ALPHA) -> dict:
    """Wilson score interval for a binomial proportion.

    Returns {"estimate", "lower", "upper", "successes", "total"} or a
    not-computable marker when total == 0 (an empty denominator is never
    reported as a zero rate).
    """
    if total <= 0:
        return {"estimate": None, "lower": None, "upper": None,
                "successes": successes, "total": total,
                "note": "no cases — rate not computable"}
    z = _z(alpha)
    p = successes / total
    denom = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denom
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denom
    return {"estimate": round(p, 4),
            "lower": round(max(0.0, center - half), 4),
            "upper": round(min(1.0, center + half), 4),
            "successes": successes, "total": total}


def _z(alpha: float) -> float:
    # Two-sided normal quantile via complementary error function inversion
    # (no scipy dependency: bisect on erfc).
    lo, hi = 0.0, 10.0
    target = alpha  # erfc(z/sqrt(2)) == alpha
    for _ in range(200):
        mid = (lo + hi) / 2
        if math.erfc(mid / math.sqrt(2)) > target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def confusion_matrix(cases: list[dict]) -> dict:
    """Confusion matrix with the defective class as positive."""
    cm = {"tp": 0, "fp": 0, "tn": 0, "fn": 0,
          "indeterminate_defective": 0, "indeterminate_acceptable": 0,
          "unknown_label": 0}
    for c in cases:
        label = c.get("label")
        verdict = c.get("verdict")
        if label not in ("acceptable", "defective"):
            cm["unknown_label"] += 1
            continue
        if verdict == "INDETERMINATE":
            cm["indeterminate_defective" if label == "defective"
               else "indeterminate_acceptable"] += 1
        elif verdict == "VETO":
            cm["tp" if label == "defective" else "fp"] += 1
        elif verdict == "PASS":
            cm["fn" if label == "defective" else "tn"] += 1
        else:
            raise ValueError(f"unexpected verdict {verdict!r} for {c.get('structure_id')}")
    return cm


def evaluate(cases: list[dict], alpha: float = DEFAULT_ALPHA) -> dict:
    """Full evaluation report: rates with Wilson intervals, per-law
    contributions, missing-measurement statistics, and strata breakdowns."""
    cm = confusion_matrix(cases)
    det_total = cm["tp"] + cm["fn"]
    acc_total = cm["tn"] + cm["fp"]

    report = {
        "confusion_matrix": cm,
        "sensitivity": wilson_interval(cm["tp"], det_total, alpha),
        "specificity": wilson_interval(cm["tn"], acc_total, alpha),
        "false_veto_rate": wilson_interval(cm["fp"], acc_total, alpha),
        "false_accept_rate": wilson_interval(cm["fn"], det_total, alpha),
        "indeterminate_rate": {
            "defective": wilson_interval(cm["indeterminate_defective"],
                                         det_total, alpha),
            "acceptable": wilson_interval(cm["indeterminate_acceptable"],
                                          acc_total, alpha),
        },
        "per_law_contribution": _per_law(cases, alpha),
        "missing_measurement_rate": _missing_rate(cases, alpha),
        "strata": _strata(cases, alpha),
        "n_cases": len(cases),
        "note": ("INDETERMINATE verdicts count toward neither sensitivity nor "
                 "specificity; they are reported separately."),
    }
    return report


def _per_law(cases: list[dict], alpha: float) -> dict:
    """How often each law is the measured blocking (veto) cause, and how
    often it is an indeterminate cause, among cases where it could act."""
    veto_counts: dict[str, int] = {}
    indet_counts: dict[str, int] = {}
    for c in cases:
        for lid in c.get("blocking_laws", []) or []:
            veto_counts[lid] = veto_counts.get(lid, 0) + 1
        for lid in c.get("indeterminate_causes", []) or []:
            indet_counts[lid] = indet_counts.get(lid, 0) + 1
    n = len(cases)
    return {
        "veto_causes": {k: {"count": v, "share_of_cases": wilson_interval(v, n, alpha)}
                        for k, v in sorted(veto_counts.items())},
        "indeterminate_causes": {k: {"count": v}
                                 for k, v in sorted(indet_counts.items())},
    }


def _missing_rate(cases: list[dict], alpha: float) -> dict:
    """Fraction of cases with any NOT_MEASURED on a deterministic law,
    and per-law NOT_MEASURED / NOT_APPLICABLE counts."""
    any_missing = sum(1 for c in cases if c.get("not_measured"))
    per_law_nm: dict[str, int] = {}
    per_law_na: dict[str, int] = {}
    for c in cases:
        for lid in c.get("not_measured", []) or []:
            per_law_nm[lid] = per_law_nm.get(lid, 0) + 1
        for lid in c.get("not_applicable", []) or []:
            per_law_na[lid] = per_law_na.get(lid, 0) + 1
    return {
        "cases_with_missing_measurement": wilson_interval(any_missing, len(cases), alpha),
        "per_law_not_measured": dict(sorted(per_law_nm.items())),
        "per_law_not_applicable": dict(sorted(per_law_na.items())),
    }


def _strata(cases: list[dict], alpha: float) -> dict:
    """Sensitivity/false-veto by protocol strata (modality, resolution band)."""
    out: dict[str, dict] = {}
    strata_keys = ("modality", "resolution_band")
    for key in strata_keys:
        groups: dict[str, list[dict]] = {}
        for c in cases:
            val = (c.get("strata") or {}).get(key)
            if val is None:
                continue
            groups.setdefault(str(val), []).append(c)
        out[key] = {}
        for val, members in sorted(groups.items()):
            sub = evaluate(members, alpha)
            out[key][val] = {
                "n": len(members),
                "sensitivity": sub["sensitivity"],
                "false_veto_rate": sub["false_veto_rate"],
                "indeterminate_rate": sub["indeterminate_rate"],
            }
    return out


def evaluate_file(path: Path, alpha: float = DEFAULT_ALPHA) -> dict:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    report = evaluate(data["cases"], alpha)
    report["policy_binding"] = data.get("policy_binding", {})
    return report


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("usage: python -m calibration.evaluate <labeled_cases.json>",
              file=sys.stderr)
        raise SystemExit(2)
    print(json.dumps(evaluate_file(Path(sys.argv[1])), indent=2))
