"""
Independent evaluation tooling for labeled calibration runs.

The evaluator consumes engine verdicts plus independently defined labels.
It never derives labels from Toscanini verdicts, and it never counts an
INDETERMINATE verdict as a sensitivity or specificity success.

Input contract (JSON)::

    {
      "policy_binding": {
        "code_sha": "<40 or 64 hex>",
        "law_semantic_hash": "<64 hex>"
      },
      "cases": [
        {
          "structure_id": "1ABC",
          "label": "acceptable" | "defective",
          "strata": {"modality": "x_ray", "resolution_band": "medium"},
          "verdict": "PASS" | "VETO" | "INDETERMINATE",
          "blocking_laws": ["LAW-110"],
          "indeterminate_causes": ["LAW-120"],
          "not_applicable": ["LAW-195"],
          "not_measured": ["LAW-120"]
        }
      ]
    }

Validation is fail-closed: duplicate IDs, unknown labels/verdicts, unknown
law IDs, retired laws used as veto causes, and inconsistent cause lists are
errors, not silently ignored data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Iterable

from tos.governance.station_sop import (
    LAW_CANON,
    LAW_METHOD_CLASSIFICATIONS,
    RETIRED_LAWS,
)

DEFAULT_ALPHA = 0.05

VALID_LABELS = {"acceptable", "defective"}
VALID_VERDICTS = {"PASS", "VETO", "INDETERMINATE"}
VALID_LABEL_SOURCES = {
    "wwpdb-validation",
    "perturbed-by-construction",
    "known-defective-deposit",
}
VALID_COHORTS = {"A", "B1", "B2"}
VALID_MODALITIES = {"x_ray", "cryo_em", "nmr", "predicted"}
VALID_RESOLUTION_BANDS = {"ultra_high", "medium", "low", "na"}
SPLIT_NAMES = ("calibration", "validation", "held_out")
HEX40 = re.compile(r"^[0-9a-f]{40}$")
HEX64 = re.compile(r"^[0-9a-f]{64}$")


class EvaluationInputError(ValueError):
    """Raised when a labeled run or split manifest is not usable evidence."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise EvaluationInputError(message)


def _hex(value, length: int, field: str) -> None:
    _require(isinstance(value, str), f"{field} must be a string")
    pattern = HEX40 if length == 40 else HEX64
    _require(bool(pattern.fullmatch(value.lower())), f"{field} must be {length} hex chars")


def _validate_law_list(case: dict, field: str, allowed: set[str]) -> list[str]:
    value = case.get(field, [])
    _require(isinstance(value, list), f"{case.get('structure_id')}: {field} must be a list")
    _require(len(value) == len(set(value)), f"{case.get('structure_id')}: {field} has duplicates")
    for law_id in value:
        _require(isinstance(law_id, str), f"{case.get('structure_id')}: {field} entries must be strings")
        _require(law_id in allowed, f"{case.get('structure_id')}: unknown law {law_id!r} in {field}")
    return value


def validate_cases(cases: Iterable[dict], require_policy: bool = False,
                   payload: dict | None = None) -> list[dict]:
    """Validate a case list and return it as a list. Never repairs data."""
    _require(isinstance(cases, (list, tuple)), "cases must be a list")
    cases = list(cases)
    _require(cases, "cases must not be empty")

    if require_policy:
        _require(isinstance(payload, dict), "payload must be an object")
        binding = payload.get("policy_binding")
        _require(isinstance(binding, dict), "policy_binding is required")
        _require(isinstance(binding.get("code_sha"), str), "policy_binding.code_sha is required")
        code_sha = binding["code_sha"].lower()
        _require(bool(HEX40.fullmatch(code_sha) or HEX64.fullmatch(code_sha)),
                 "policy_binding.code_sha must be a full 40- or 64-hex commit")
        _hex(binding.get("law_semantic_hash"), 64, "policy_binding.law_semantic_hash")

    seen_ids: set[str] = set()
    for case in cases:
        _require(isinstance(case, dict), "each case must be an object")
        sid = case.get("structure_id")
        _require(isinstance(sid, str) and sid.strip(), "case.structure_id must be non-empty")
        sid = sid.strip()
        _require(sid not in seen_ids, f"duplicate structure_id: {sid}")
        seen_ids.add(sid)

        label = case.get("label")
        _require(label in VALID_LABELS, f"{sid}: invalid label {label!r}")
        verdict = case.get("verdict")
        _require(verdict in VALID_VERDICTS, f"{sid}: invalid verdict {verdict!r}")

        strata = case.get("strata", {})
        _require(isinstance(strata, dict), f"{sid}: strata must be an object")
        if "modality" in strata:
            _require(strata["modality"] in VALID_MODALITIES,
                     f"{sid}: invalid strata.modality")
        if "resolution_band" in strata:
            _require(strata["resolution_band"] in VALID_RESOLUTION_BANDS,
                     f"{sid}: invalid strata.resolution_band")

        active_laws = set(LAW_CANON) - RETIRED_LAWS
        blocking = _validate_law_list(case, "blocking_laws", active_laws)
        causes = _validate_law_list(case, "indeterminate_causes", active_laws)
        not_measured = _validate_law_list(case, "not_measured", set(LAW_CANON) - RETIRED_LAWS)
        not_applicable = _validate_law_list(case, "not_applicable", set(LAW_CANON))

        _require(not set(not_measured) & set(not_applicable),
                 f"{sid}: a law cannot be both NOT_MEASURED and NOT_APPLICABLE")
        _require(not set(blocking) & set(not_measured),
                 f"{sid}: a blocking law cannot also be NOT_MEASURED")

        if verdict == "VETO":
            _require(blocking, f"{sid}: VETO requires at least one blocking law")
        else:
            _require(not blocking, f"{sid}: blocking_laws are only valid for VETO")

        if verdict == "INDETERMINATE":
            _require(set(causes) <= set(not_measured),
                     f"{sid}: indeterminate causes must also be listed as NOT_MEASURED")
        else:
            _require(not causes, f"{sid}: indeterminate_causes are only valid for INDETERMINATE")

        for law_id in blocking + causes:
            _require(LAW_METHOD_CLASSIFICATIONS[law_id] == "deterministic",
                     f"{sid}: {law_id} is not a deterministic law and cannot block")

    return cases


def wilson_interval(successes: int, total: int,
                    alpha: float = DEFAULT_ALPHA) -> dict:
    """Wilson score interval for a binomial proportion."""
    _require(0.0 < alpha < 1.0, "alpha must be between 0 and 1")
    _require(isinstance(successes, int) and successes >= 0, "successes must be a non-negative integer")
    _require(isinstance(total, int) and total >= 0, "total must be a non-negative integer")
    _require(successes <= total, "successes cannot exceed total")
    if total == 0:
        return {
            "estimate": None, "lower": None, "upper": None,
            "successes": successes, "total": total,
            "note": "no eligible cases - rate not computable",
        }
    z = _z(alpha)
    p = successes / total
    denom = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denom
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denom
    return {
        "estimate": round(p, 4),
        "lower": round(max(0.0, center - half), 4),
        "upper": round(min(1.0, center + half), 4),
        "successes": successes,
        "total": total,
    }


def _z(alpha: float) -> float:
    lo, hi = 0.0, 10.0
    target = alpha
    for _ in range(200):
        mid = (lo + hi) / 2
        if math.erfc(mid / math.sqrt(2)) > target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def confusion_matrix(cases: list[dict]) -> dict:
    cm = {
        "tp": 0, "fp": 0, "tn": 0, "fn": 0,
        "indeterminate_defective": 0,
        "indeterminate_acceptable": 0,
    }
    for case in cases:
        label, verdict = case["label"], case["verdict"]
        if verdict == "INDETERMINATE":
            key = "indeterminate_defective" if label == "defective" else "indeterminate_acceptable"
            cm[key] += 1
        elif verdict == "VETO":
            cm["tp" if label == "defective" else "fp"] += 1
        else:
            cm["fn" if label == "defective" else "tn"] += 1
    return cm


def _summary(cases: list[dict], alpha: float) -> dict:
    cm = confusion_matrix(cases)
    defective_determinate = cm["tp"] + cm["fn"]
    acceptable_determinate = cm["tn"] + cm["fp"]
    defective_total = defective_determinate + cm["indeterminate_defective"]
    acceptable_total = acceptable_determinate + cm["indeterminate_acceptable"]
    return {
        "n": len(cases),
        "confusion_matrix": cm,
        "sensitivity": wilson_interval(cm["tp"], defective_determinate, alpha),
        "specificity": wilson_interval(cm["tn"], acceptable_determinate, alpha),
        "false_veto_rate": wilson_interval(cm["fp"], acceptable_determinate, alpha),
        "false_accept_rate": wilson_interval(cm["fn"], defective_determinate, alpha),
        "indeterminate_rate": {
            "defective": wilson_interval(cm["indeterminate_defective"], defective_total, alpha),
            "acceptable": wilson_interval(cm["indeterminate_acceptable"], acceptable_total, alpha),
        },
    }


def _per_law(cases: list[dict], alpha: float) -> dict:
    veto_cases = [c for c in cases if c["verdict"] == "VETO"]
    indeterminate_cases = [c for c in cases if c["verdict"] == "INDETERMINATE"]
    veto_counts: dict[str, int] = {}
    indet_counts: dict[str, int] = {}
    for case in veto_cases:
        for law_id in case.get("blocking_laws", []):
            veto_counts[law_id] = veto_counts.get(law_id, 0) + 1
    for case in indeterminate_cases:
        for law_id in case.get("indeterminate_causes", []):
            indet_counts[law_id] = indet_counts.get(law_id, 0) + 1

    return {
        "veto_causes": {
            law_id: {
                "count": count,
                "share_of_veto_cases": wilson_interval(count, len(veto_cases), alpha),
            }
            for law_id, count in sorted(veto_counts.items())
        },
        "indeterminate_causes": {
            law_id: {
                "count": count,
                "share_of_indeterminate_cases": wilson_interval(
                    count, len(indeterminate_cases), alpha),
            }
            for law_id, count in sorted(indet_counts.items())
        },
        "denominators": {
            "veto_cases": len(veto_cases),
            "indeterminate_cases": len(indeterminate_cases),
            "all_cases": len(cases),
        },
    }


def _missing_rate(cases: list[dict], alpha: float) -> dict:
    n_cases = len(cases)
    any_missing = sum(1 for case in cases if case.get("not_measured"))
    active_laws = set(LAW_CANON) - RETIRED_LAWS

    per_law: dict[str, dict] = {}
    missing_opportunities = 0
    eligible_opportunities = 0
    for law_id in sorted(active_laws):
        not_applicable = sum(1 for case in cases if law_id in case.get("not_applicable", []))
        eligible = n_cases - not_applicable
        missing = sum(1 for case in cases if law_id in case.get("not_measured", []))
        eligible = max(eligible, missing)
        missing_opportunities += missing
        eligible_opportunities += eligible
        per_law[law_id] = {
            "not_measured": missing,
            "eligible": eligible,
            "not_applicable": not_applicable,
            "rate": wilson_interval(missing, eligible, alpha),
        }

    return {
        "cases_with_any_missing": wilson_interval(any_missing, n_cases, alpha),
        "law_opportunities": {
            "missing": missing_opportunities,
            "eligible": eligible_opportunities,
            "rate": wilson_interval(missing_opportunities, eligible_opportunities, alpha),
        },
        "per_law": per_law,
    }


def _strata(cases: list[dict], alpha: float) -> dict:
    out: dict[str, dict] = {}
    for key in ("modality", "resolution_band"):
        groups: dict[str, list[dict]] = {}
        unstratified = 0
        for case in cases:
            value = (case.get("strata") or {}).get(key)
            if value is None:
                unstratified += 1
                continue
            groups.setdefault(str(value), []).append(case)
        out[key] = {
            "unstratified": unstratified,
            "groups": {
                value: _summary(members, alpha)
                for value, members in sorted(groups.items())
            },
        }
    return out


def evaluate(cases: list[dict], alpha: float = DEFAULT_ALPHA) -> dict:
    cases = validate_cases(cases, require_policy=False)
    summary = _summary(cases, alpha)
    return {
        "n_cases": len(cases),
        "alpha": alpha,
        "confusion_matrix": summary["confusion_matrix"],
        "sensitivity": summary["sensitivity"],
        "specificity": summary["specificity"],
        "false_veto_rate": summary["false_veto_rate"],
        "false_accept_rate": summary["false_accept_rate"],
        "indeterminate_rate": summary["indeterminate_rate"],
        "per_law_contribution": _per_law(cases, alpha),
        "missing_measurement_rate": _missing_rate(cases, alpha),
        "strata": _strata(cases, alpha),
        "note": (
            "INDETERMINATE verdicts count toward neither sensitivity nor "
            "specificity; they are reported separately."
        ),
    }


def validate_split_manifest(manifest: dict) -> dict:
    """Validate a frozen split manifest (protocol v1, section 3)."""
    _require(isinstance(manifest, dict), "split manifest must be an object")
    _require(manifest.get("schema_version") == "split-manifest/1.0.0",
             "unsupported split manifest schema_version")
    _require(isinstance(manifest.get("protocol"), str) and manifest["protocol"].strip(),
             "split manifest protocol is required")
    binding = manifest.get("policy_binding")
    _require(isinstance(binding, dict), "split manifest policy_binding is required")
    code_sha = binding.get("code_sha")
    _require(isinstance(code_sha, str)
             and bool(HEX40.fullmatch(code_sha.lower()) or HEX64.fullmatch(code_sha.lower())),
             "split manifest policy_binding.code_sha must be a full 40- or 64-hex commit")
    _hex(binding.get("law_semantic_hash"), 64,
         "split manifest policy_binding.law_semantic_hash")
    splits = manifest.get("splits")
    _require(isinstance(splits, dict), "splits must be an object")
    for name in SPLIT_NAMES:
        _require(isinstance(splits.get(name), list), f"splits.{name} must be a list")
        _require(splits[name], f"splits.{name} must not be empty")

    seen: dict[str, str] = {}
    summary: dict[str, int] = {}
    for split_name in SPLIT_NAMES:
        entries = splits[split_name]
        summary[split_name] = len(entries)
        for entry in entries:
            _require(isinstance(entry, dict), f"{split_name}: entry must be an object")
            sid = entry.get("structure_id")
            _require(isinstance(sid, str) and sid.strip(),
                     f"{split_name}: structure_id is required")
            sid = sid.strip()
            _require(sid not in seen,
                     f"structure identity appears in more than one split: {sid} ({seen.get(sid)} and {split_name})")
            seen[sid] = split_name
            _require(entry.get("label") in VALID_LABELS,
                     f"{sid}: invalid label {entry.get('label')!r}")
            _require(entry.get("label_source") in VALID_LABEL_SOURCES,
                     f"{sid}: invalid label_source {entry.get('label_source')!r}")
            cohort = entry.get("cohort")
            _require(cohort in VALID_COHORTS,
                     f"{sid}: invalid cohort {cohort!r}")
            expected_label_source = {
                "A": "wwpdb-validation",
                "B1": "perturbed-by-construction",
                "B2": "known-defective-deposit",
            }[cohort]
            expected_label = "acceptable" if cohort == "A" else "defective"
            _require(entry.get("label_source") == expected_label_source,
                     f"{sid}: cohort {cohort} requires label_source "
                     f"{expected_label_source!r}")
            _require(entry.get("label") == expected_label,
                     f"{sid}: cohort {cohort} requires label {expected_label!r}")
            strata = entry.get("strata")
            _require(isinstance(strata, dict), f"{sid}: strata is required")
            _require(strata.get("modality") in VALID_MODALITIES,
                     f"{sid}: invalid modalities")
            _require(strata.get("resolution_band") in VALID_RESOLUTION_BANDS,
                     f"{sid}: invalid resolution_band")
            _hex(entry.get("coordinate_sha256"), 64, f"{sid}.coordinate_sha256")
            if entry["cohort"] == "B1":
                _require(isinstance(entry.get("defect_manifest_ref"), str)
                         and entry["defect_manifest_ref"].strip(),
                         f"{sid}: B1 requires defect_manifest_ref")

    checks = manifest.get("overlap_checks")
    _require(isinstance(checks, dict), "overlap_checks is required")
    _require(checks.get("identity_overlap_found") is False,
             "identity_overlap_found must be false before freezing")
    _require(checks.get("all_checks_passed") is True,
             "all overlap checks must pass before freezing")
    pairs = checks.get("near_identity_pairs")
    _require(isinstance(pairs, list), "near_identity_pairs must be a list")
    _require(not pairs, "near-identity pairs must be resolved before freezing")
    threshold = checks.get("near_identity_threshold")
    _require(isinstance(threshold, (int, float)) and not isinstance(threshold, bool)
             and 0.9 <= threshold <= 1,
             "near_identity_threshold must be a proportion >= 0.9")

    return {
        "split_counts": summary,
        "total_entries": sum(summary.values()),
        "overlap_checks_passed": True,
    }


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def evaluate_file(path: Path, alpha: float = DEFAULT_ALPHA) -> dict:
    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    _require(isinstance(data, dict), "evaluation file must contain an object")
    cases = validate_cases(data.get("cases"), require_policy=True, payload=data)
    report = evaluate(cases, alpha)
    report["policy_binding"] = data["policy_binding"]
    report["input_file_sha256"] = _file_sha256(path)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("input", help="labeled cases JSON or split manifest JSON")
    parser.add_argument("--alpha", type=float, default=DEFAULT_ALPHA)
    parser.add_argument("--validate-split", action="store_true",
                        help="validate a frozen split manifest instead of cases")
    args = parser.parse_args(argv)

    path = Path(args.input)
    if args.validate_split:
        payload = json.loads(path.read_text(encoding="utf-8"))
        summary = validate_split_manifest(payload)
        summary["input_file_sha256"] = _file_sha256(path)
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0

    report = evaluate_file(path, args.alpha)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
