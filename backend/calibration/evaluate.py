"""
Independent evaluation tooling for labeled calibration runs.

The evaluator consumes engine verdicts plus independently defined labels.
It never derives labels from Toscanini verdicts, and it never counts an
INDETERMINATE verdict as a sensitivity or specificity success.

Input contract (JSON)::

    {
      "policy_binding": {
        "code_sha": "<40 or 64 hex>",
        "canon_hash": "<16 hex>",
        "law_semantic_hash": "<64 hex>",
        "modality_matrix_sha256": "<64 hex>",
        "law_registry_sha256": "<64 hex>"
      },
      "run_inputs": {
        "run_manifest_sha256": "<64 hex>",
        "coordinate_manifest_sha256": "<64 hex>",
        "calibration_results_sha256": "<64 hex>"
      },
      "cases": [
        {
          "structure_id": "1ABC",
          "label": "acceptable" | "defective",
          "label_source": "wwpdb-validation" | "perturbed-by-construction"
                          | "known-defective-deposit",
          "label_evidence_sha256": "<64 hex>",
          "cohort": "A" | "B1" | "B2",
          "strata": {"modality": "x_ray", "resolution_band": "medium"},
          "verdict": "PASS" | "VETO" | "INDETERMINATE",
          "blocking_laws": ["LAW-110"],
          "indeterminate_causes": ["LAW-120"],
          "not_applicable": ["LAW-195"],
          "not_measured": ["LAW-120"],
          "measurement_failed": ["LAW-130"]
        }
      ]
    }

Validation is fail-closed: duplicate IDs, unknown labels/verdicts, unknown
law IDs, retired laws used as veto causes, missing label provenance,
unverifiable policy identities, and inconsistent cause lists are errors,
not silently ignored data. `NOT_MEASURED` and `MEASUREMENT_FAILED` are
distinct states and both count as missing evidence; neither is a pass.
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
    LAW_CANON_HASH,
    LAW_METHOD_CLASSIFICATIONS,
    RETIRED_LAWS,
)
from tos.governance.modality_matrix import compute_matrix_hash

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
COHORT_LABELS = {
    "A": ("acceptable", "wwpdb-validation"),
    "B1": ("defective", "perturbed-by-construction"),
    "B2": ("defective", "known-defective-deposit"),
}
HEX16 = re.compile(r"^[0-9a-f]{16}$")
HEX40 = re.compile(r"^[0-9a-f]{40}$")
HEX64 = re.compile(r"^[0-9a-f]{64}$")
STRUCTURE_ID = re.compile(r"^[A-Z0-9][A-Z0-9_-]{2,24}$")
RUN_INPUT_FILES = {
    "run_manifest_sha256": "run_manifest.json",
    "coordinate_manifest_sha256": "coordinate_manifest.json",
    "calibration_results_sha256": "calibration_results.json",
}


class EvaluationInputError(ValueError):
    """Raised when a labeled run or split manifest is not usable evidence."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise EvaluationInputError(message)


def _hex(value, length: int, field: str) -> None:
    _require(isinstance(value, str), f"{field} must be a string")
    pattern = {16: HEX16, 40: HEX40, 64: HEX64}[length]
    _require(bool(pattern.fullmatch(value.lower())), f"{field} must be {length} hex chars")


def _normalize_structure_id(value) -> str:
    _require(isinstance(value, str), "case.structure_id must be non-empty")
    sid = value.strip().upper()
    _require(bool(sid), "case.structure_id must be non-empty")
    _require(bool(STRUCTURE_ID.fullmatch(sid)),
             f"{sid}: structure_id contains unsupported characters")
    return sid


def _validate_law_list(case: dict, field: str, allowed: set[str]) -> list[str]:
    sid = case["structure_id"]
    value = case.get(field, [])
    _require(isinstance(value, list), f"{sid}: {field} must be a list")
    _require(len(value) == len(set(value)), f"{sid}: {field} has duplicates")
    for law_id in value:
        _require(isinstance(law_id, str), f"{sid}: {field} entries must be strings")
        _require(law_id in allowed, f"{sid}: unknown law {law_id!r} in {field}")
    return value


def _validate_label_provenance(case: dict) -> None:
    sid = case["structure_id"]
    label = case["label"]
    cohort = case.get("cohort")
    _require(cohort in VALID_COHORTS, f"{sid}: invalid or missing cohort")
    expected_label, expected_source = COHORT_LABELS[cohort]
    _require(label == expected_label,
             f"{sid}: cohort {cohort} requires label {expected_label!r}")
    source = case.get("label_source")
    _require(source in VALID_LABEL_SOURCES,
             f"{sid}: invalid or missing label_source")
    _require(source == expected_source,
             f"{sid}: cohort {cohort} requires label_source {expected_source!r}")
    _hex(case.get("label_evidence_sha256"), 64,
         f"{sid}.label_evidence_sha256")
    _require(case.get("label_evidence_ref") and isinstance(case.get("label_evidence_ref"), str),
             f"{sid}: label evidence reference is required for independent verification")
    if cohort == "B1":
        ref = case.get("defect_manifest_ref")
        _require(isinstance(ref, str) and ref.strip(),
                 f"{sid}: B1 requires defect_manifest_ref")
        _hex(case.get("defect_manifest_sha256"), 64,
             f"{sid}.defect_manifest_sha256")


def _validate_policy_binding(payload: dict) -> dict:
    binding = payload.get("policy_binding")
    _require(isinstance(binding, dict), "policy_binding is required")
    code_sha = binding.get("code_sha")
    _require(isinstance(code_sha, str), "policy_binding.code_sha is required")
    code_sha = code_sha.lower()
    _require(bool(HEX40.fullmatch(code_sha) or HEX64.fullmatch(code_sha)),
             "policy_binding.code_sha must be a full 40- or 64-hex commit")
    binding["code_sha"] = code_sha
    for field, length in (
        ("canon_hash", 16),
        ("law_semantic_hash", 64),
        ("modality_matrix_sha256", 64),
        ("law_registry_sha256", 64),
    ):
        _hex(binding.get(field), length, f"policy_binding.{field}")
        binding[field] = binding[field].lower()

    run_inputs = payload.get("run_inputs")
    _require(isinstance(run_inputs, dict), "run_inputs is required")
    for field in RUN_INPUT_FILES:
        _hex(run_inputs.get(field), 64, f"run_inputs.{field}")
        run_inputs[field] = run_inputs[field].lower()
    return binding


def _validate_policy_identity(binding: dict) -> None:
    _require(binding["canon_hash"] == LAW_CANON_HASH,
             "policy_binding.canon_hash does not match the running engine")
    _require(binding["modality_matrix_sha256"] == compute_matrix_hash(),
             "policy_binding.modality_matrix_sha256 does not match the running engine")
    _require(binding["law_semantic_hash"] == _law_semantic_hash(),
             "policy_binding.law_semantic_hash does not match the running engine")
    registry = Path(__file__).resolve().parents[2] / "governance" / "law_registry.yaml"
    _require(registry.is_file(), "governance/law_registry.yaml is missing")
    _require(binding["law_registry_sha256"] == _file_sha256(registry),
             "policy_binding.law_registry_sha256 does not match the live registry")


def validate_cases(cases: Iterable[dict], require_policy: bool = False,
                   payload: dict | None = None,
                   require_label_provenance: bool = False) -> list[dict]:
    """Validate a case list and return it as a list. Never repairs data."""
    _require(isinstance(cases, (list, tuple)), "cases must be a list")
    cases = list(cases)
    _require(cases, "cases must not be empty")

    if require_policy:
        _require(isinstance(payload, dict), "payload must be an object")
        _validate_policy_binding(payload)
        _validate_policy_identity(payload["policy_binding"])

    seen_ids: set[str] = set()
    active_laws = set(LAW_CANON) - RETIRED_LAWS
    for case in cases:
        _require(isinstance(case, dict), "each case must be an object")
        sid = _normalize_structure_id(case.get("structure_id"))
        case["structure_id"] = sid
        _require(sid not in seen_ids, f"duplicate structure_id: {sid}")
        seen_ids.add(sid)

        label = case.get("label")
        _require(label in VALID_LABELS, f"{sid}: invalid label {label!r}")
        verdict = case.get("verdict")
        _require(verdict in VALID_VERDICTS, f"{sid}: invalid verdict {verdict!r}")
        if require_label_provenance:
            _validate_label_provenance(case)

        strata = case.get("strata", {})
        _require(isinstance(strata, dict), f"{sid}: strata must be an object")
        if "modality" in strata:
            _require(strata["modality"] in VALID_MODALITIES,
                     f"{sid}: invalid strata.modality")
        if "resolution_band" in strata:
            _require(strata["resolution_band"] in VALID_RESOLUTION_BANDS,
                     f"{sid}: invalid strata.resolution_band")

        blocking = _validate_law_list(case, "blocking_laws", active_laws)
        causes = _validate_law_list(case, "indeterminate_causes", active_laws)
        not_measured = _validate_law_list(case, "not_measured", active_laws)
        measurement_failed = _validate_law_list(case, "measurement_failed", active_laws)
        not_applicable = _validate_law_list(case, "not_applicable", active_laws)

        missing = set(not_measured) | set(measurement_failed)
        _require(len(not_measured) + len(measurement_failed) == len(missing),
                 f"{sid}: a law cannot be both NOT_MEASURED and MEASUREMENT_FAILED")
        _require(not set(missing) & set(not_applicable),
                 f"{sid}: a law cannot be both missing and NOT_APPLICABLE")
        _require(not set(blocking) & missing,
                 f"{sid}: a blocking law cannot also be missing")
        _require(not set(blocking) & set(not_applicable),
                 f"{sid}: a blocking law cannot also be NOT_APPLICABLE")
        _require(not set(causes) & set(not_applicable),
                 f"{sid}: an indeterminate cause cannot be NOT_APPLICABLE")

        if verdict == "VETO":
            _require(blocking, f"{sid}: VETO requires at least one blocking law")
        else:
            _require(not blocking, f"{sid}: blocking_laws are only valid for VETO")

        if verdict == "INDETERMINATE":
            _require(causes, f"{sid}: INDETERMINATE requires at least one cause")
            _require(set(causes) <= missing,
                     f"{sid}: indeterminate causes must be listed as missing evidence")
        else:
            _require(not causes, f"{sid}: indeterminate_causes are only valid for INDETERMINATE")

        if verdict == "PASS":
            _require(not missing,
                     f"{sid}: PASS is invalid while deterministic evidence is missing")

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
    failure_counts: dict[str, int] = {}
    for case in veto_cases:
        for law_id in case.get("blocking_laws", []):
            veto_counts[law_id] = veto_counts.get(law_id, 0) + 1
    for case in indeterminate_cases:
        for law_id in case.get("indeterminate_causes", []):
            indet_counts[law_id] = indet_counts.get(law_id, 0) + 1
    for case in cases:
        for law_id in case.get("measurement_failed", []):
            failure_counts[law_id] = failure_counts.get(law_id, 0) + 1

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
        "measurement_failed_causes": {
            law_id: {
                "count": count,
                "share_of_all_cases": wilson_interval(count, len(cases), alpha),
            }
            for law_id, count in sorted(failure_counts.items())
        },
        "denominators": {
            "veto_cases": len(veto_cases),
            "indeterminate_cases": len(indeterminate_cases),
            "all_cases": len(cases),
        },
    }


def _missing_rate(cases: list[dict], alpha: float) -> dict:
    n_cases = len(cases)
    any_missing = sum(
        1 for case in cases
        if case.get("not_measured") or case.get("measurement_failed")
    )
    active_laws = set(LAW_CANON) - RETIRED_LAWS

    per_law: dict[str, dict] = {}
    missing_opportunities = 0
    eligible_opportunities = 0
    for law_id in sorted(active_laws):
        not_applicable = sum(1 for case in cases if law_id in case.get("not_applicable", []))
        eligible = n_cases - not_applicable
        not_measured = sum(1 for case in cases if law_id in case.get("not_measured", []))
        measurement_failed = sum(
            1 for case in cases if law_id in case.get("measurement_failed", []))
        missing = not_measured + measurement_failed
        _require(missing <= eligible,
                 f"{law_id}: missing count exceeds eligible opportunities")
        missing_opportunities += missing
        eligible_opportunities += eligible
        per_law[law_id] = {
            "not_measured": not_measured,
            "measurement_failed": measurement_failed,
            "missing": missing,
            "eligible": eligible,
            "not_applicable": not_applicable,
            "rate": wilson_interval(missing, eligible, alpha),
        }

    return {
        "cases_with_any_missing": wilson_interval(any_missing, n_cases, alpha),
        "status_counts": {
            "not_measured": sum(1 for case in cases if case.get("not_measured")),
            "measurement_failed": sum(
                1 for case in cases if case.get("measurement_failed")),
        },
        "law_opportunities": {
            "missing": missing_opportunities,
            "not_measured": sum(
                count["not_measured"] for count in per_law.values()),
            "measurement_failed": sum(
                count["measurement_failed"] for count in per_law.values()),
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
        "report_type": "diagnostic",
        "calibration_status": "not_approved",
        "supports_threshold_validation": False,
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


SEQUENCE_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
IDENTITY_METHOD = "normalized-global-edit-identity"
IDENTITY_TOOL = "calibration.evaluate"
IDENTITY_TOOL_VERSION = "1.0.0"


def _sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _bounded_edit_distance(a: str, b: str, max_errors: int) -> int | None:
    if abs(len(a) - len(b)) > max_errors:
        return None
    if max_errors == 0:
        return 0 if a == b else None
    n, m = len(a), len(b)
    infinity = max_errors + 1
    prev = [infinity] * (m + 1)
    for j in range(min(m, max_errors) + 1):
        prev[j] = j
    for i in range(1, n + 1):
        current = [infinity] * (m + 1)
        if i <= max_errors:
            current[0] = i
        lo = max(1, i - max_errors)
        hi = min(m, i + max_errors)
        row_min = infinity
        for j in range(lo, hi + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            value = min(prev[j] + 1, current[j - 1] + 1, prev[j - 1] + cost)
            current[j] = value
            if value < row_min:
                row_min = value
        if row_min > max_errors:
            return None
        prev = current
    distance = prev[m]
    return distance if distance <= max_errors else None


def _sequence_identity(a: str, b: str, threshold: float) -> float | None:
    if a == b:
        return 1.0
    max_length = max(len(a), len(b))
    max_errors = int(math.floor((1.0 - threshold) * max_length + 1e-9))
    distance = _bounded_edit_distance(a, b, max_errors)
    if distance is None:
        return None
    return round(1.0 - distance / max_length, 9)


def _pair_results(entries_by_split: dict[str, list[dict]], threshold: float) -> tuple:
    pairs = []
    results = []
    for left_index, left_name in enumerate(SPLIT_NAMES):
        for right_name in SPLIT_NAMES[left_index + 1:]:
            for left in entries_by_split[left_name]:
                for right in entries_by_split[right_name]:
                    pair = sorted((left["structure_id"], right["structure_id"]))
                    pairs.append(pair)
                    identity = _sequence_identity(
                        left["_sequence"], right["_sequence"], threshold)
                    if identity is None:
                        identity = "below_threshold"
                    results.append([pair[0], pair[1], identity])
    pairs.sort()
    results.sort(key=lambda row: (row[0], row[1]))
    near = [
        {"a": a, "b": b, "identity": identity}
        for a, b, identity in results
        if isinstance(identity, float) and identity >= threshold
    ]
    pair_digest = hashlib.sha256(
        json.dumps(pairs, sort_keys=True).encode("utf-8")).hexdigest()
    result_digest = hashlib.sha256(
        json.dumps(results, sort_keys=True).encode("utf-8")).hexdigest()
    return pair_digest, result_digest, near


def validate_split_manifest(manifest: dict) -> dict:
    """Validate a frozen split manifest and independently recompute leakage."""
    _require(isinstance(manifest, dict), "split manifest must be an object")
    _require(manifest.get("schema_version") == "split-manifest/1.1.0",
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

    exclusions = manifest.get("historical_exclusions", [])
    _require(isinstance(exclusions, list) and all(
        isinstance(item, str) and item.strip() for item in exclusions),
        "historical_exclusions must be a list of structure IDs")
    excluded = {item.strip().upper() for item in exclusions}
    exclusion_provenance = manifest.get("historical_exclusion_provenance")
    _require(isinstance(exclusion_provenance, dict),
             "historical_exclusion_provenance is required")
    _require(exclusion_provenance.get("inventory_ref") and
             isinstance(exclusion_provenance.get("inventory_sha256"), str),
             "historical exclusion inventory reference and digest are required")
    _hex(exclusion_provenance["inventory_sha256"], 64,
         "historical_exclusion_provenance.inventory_sha256")
    _require(
        {item.strip().upper() for item in exclusion_provenance.get("excluded_structure_ids", [])}
        == excluded,
        "historical_exclusions do not match the verified exclusion inventory")

    seen: dict[str, str] = {}
    seen_digests: dict[str, str] = {}
    seen_sequences: dict[str, str] = {}
    families: dict[str, set[str]] = {}
    entries_by_split: dict[str, list[dict]] = {}
    summary: dict[str, int] = {}
    for split_name in SPLIT_NAMES:
        entries = splits[split_name]
        summary[split_name] = len(entries)
        entries_by_split[split_name] = []
        families[split_name] = set()
        for entry in entries:
            _require(isinstance(entry, dict), f"{split_name}: entry must be an object")
            sid = _normalize_structure_id(entry.get("structure_id"))
            entry["structure_id"] = sid
            _require(sid not in excluded,
                     f"{sid}: appears in historical exclusions")
            _require(sid not in seen,
                     f"structure identity appears in more than one split: {sid} ({seen.get(sid)} and {split_name})")
            seen[sid] = split_name

            _require(entry.get("label") in VALID_LABELS,
                     f"{sid}: invalid label {entry.get('label')!r}")
            _require(entry.get("label_source") in VALID_LABEL_SOURCES,
                     f"{sid}: invalid label_source {entry.get('label_source')!r}")
            _validate_label_provenance(entry)
            cohort = entry["cohort"]
            if cohort == "B1":
                _hex(entry.get("defect_manifest_sha256"), 64,
                     f"{sid}.defect_manifest_sha256")

            strata = entry.get("strata")
            _require(isinstance(strata, dict), f"{sid}: strata is required")
            _require(strata.get("modality") in VALID_MODALITIES,
                     f"{sid}: invalid modalities")
            _require(strata.get("resolution_band") in VALID_RESOLUTION_BANDS,
                     f"{sid}: invalid resolution_band")
            coordinate_digest = entry.get("coordinate_sha256")
            _hex(coordinate_digest, 64, f"{sid}.coordinate_sha256")
            coordinate_digest = coordinate_digest.lower()
            _require(coordinate_digest not in seen_digests,
                     f"{sid}: repeated coordinate digest also used by "
                     f"{seen_digests.get(coordinate_digest)}")
            seen_digests[coordinate_digest] = sid

            sequence = entry.get("sequence")
            _require(isinstance(sequence, str) and sequence,
                     f"{sid}: sequence is required")
            sequence = sequence.strip().upper()
            _require(set(sequence) <= SEQUENCE_ALPHABET,
                     f"{sid}: sequence contains unsupported residues")
            _require(entry.get("sequence_sha256") == _sequence_sha256(sequence),
                     f"{sid}: sequence_sha256 does not match sequence")
            _require(entry.get("sequence_source_ref") and
                     entry.get("sequence_source_sha256"),
                     f"{sid}: sequence source reference and digest are required")
            _hex(entry["sequence_source_sha256"], 64,
                 f"{sid}.sequence_source_sha256")
            _require(sequence not in seen_sequences,
                     f"{sid}: repeated sequence also used by "
                     f"{seen_sequences.get(sequence)}")
            seen_sequences[sequence] = sid
            family = entry.get("family_id")
            _require(isinstance(family, str) and family.strip(),
                     f"{sid}: family_id is required")
            _require(entry.get("family_source_ref") and
                     entry.get("family_source_sha256"),
                     f"{sid}: family source reference and digest are required")
            _hex(entry["family_source_sha256"], 64,
                 f"{sid}.family_source_sha256")
            families[split_name].add(family.strip())

            entries_by_split[split_name].append({
                **entry,
                "_sequence": sequence,
            })

    for index, left_name in enumerate(SPLIT_NAMES):
        for right_name in SPLIT_NAMES[index + 1:]:
            overlap = families[left_name] & families[right_name]
            _require(not overlap,
                     f"family overlap between {left_name} and {right_name}: "
                     f"{sorted(overlap)}")

    checks = manifest.get("overlap_checks")
    _require(isinstance(checks, dict), "overlap_checks is required")
    sequence_checks = checks.get("sequence_identity")
    _require(isinstance(sequence_checks, dict),
             "overlap_checks.sequence_identity is required")
    _require(sequence_checks.get("method") == IDENTITY_METHOD,
             "unsupported sequence identity method")
    _require(sequence_checks.get("tool") == IDENTITY_TOOL,
             "unsupported sequence identity tool")
    _require(sequence_checks.get("tool_version") == IDENTITY_TOOL_VERSION,
             "unsupported sequence identity tool_version")
    threshold = sequence_checks.get("near_identity_threshold")
    _require(isinstance(threshold, (int, float)) and not isinstance(threshold, bool)
             and 0.9 <= threshold <= 1,
             "near_identity_threshold must be a proportion >= 0.9")
    pair_digest, result_digest, near_pairs = _pair_results(
        entries_by_split, float(threshold))
    expected_pairs = (
        len(entries_by_split["calibration"]) * len(entries_by_split["validation"])
        + len(entries_by_split["calibration"]) * len(entries_by_split["held_out"])
        + len(entries_by_split["validation"]) * len(entries_by_split["held_out"])
    )
    _require(sequence_checks.get("checked_pairs") == expected_pairs,
             "sequence_identity.checked_pairs does not match recomputed pairs")
    _require(sequence_checks.get("checked_pair_digest") == pair_digest,
             "sequence_identity.checked_pair_digest does not match recomputed pairs")
    _require(sequence_checks.get("result_digest") == result_digest,
             "sequence_identity.result_digest does not match recomputed identity results")
    _require(checks.get("near_identity_pairs") == near_pairs,
             "near_identity_pairs does not match recomputed identity results")
    _require(not near_pairs,
             "near-identity pairs must be resolved before freezing")
    _require(checks.get("identity_overlap_found") is False,
             "identity_overlap_found must be false before freezing")
    _require(checks.get("family_overlap_found") is False,
             "family_overlap_found must be false before freezing")
    _require(checks.get("all_checks_passed") is True,
             "all overlap checks must pass before freezing")

    return {
        "split_counts": summary,
        "total_entries": sum(summary.values()),
        "overlap_checks_passed": True,
        "checked_pairs": expected_pairs,
        "near_identity_pair_count": len(near_pairs),
        "checked_pair_digest": pair_digest,
        "identity_result_digest": result_digest,
    }


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _law_semantic_hash() -> str:
    payload = {
        "laws": {
            law_id: {
                "title": data["title"],
                "threshold": data["threshold"],
                "unit": data["unit"],
                "operator": data["operator"],
                "scope": data["scope"],
                "type": data["type"],
            }
            for law_id, data in sorted(LAW_CANON.items())
        },
        "method_classifications": dict(sorted(LAW_METHOD_CLASSIFICATIONS.items())),
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def _verify_run_inputs(run_inputs: dict, run_dir: Path,
                       policy_binding: dict) -> None:
    run_dir = Path(run_dir)
    _require(run_dir.is_dir(), f"run directory does not exist: {run_dir}")
    for field, filename in RUN_INPUT_FILES.items():
        path = run_dir / filename
        _require(path.is_file(), f"run input is missing: {filename}")
        _require(_file_sha256(path) == run_inputs[field],
                 f"run input digest mismatch: {filename}")

    integrity_path = run_dir / "run_integrity_manifest.json"
    if integrity_path.is_file():
        integrity = json.loads(integrity_path.read_text(encoding="utf-8"))
        _require(isinstance(integrity, dict), "run integrity manifest must be an object")
        integrity_policy = integrity.get("policy_binding")
        _require(isinstance(integrity_policy, dict),
                 "run integrity manifest policy_binding is missing")
        for field in ("canon_hash", "law_semantic_hash",
                      "modality_matrix_sha256", "law_registry_sha256"):
            _require(integrity_policy.get(field) == policy_binding[field],
                     f"run integrity manifest policy_binding.{field} mismatch")

    manifest = json.loads(
        (run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    _require(isinstance(manifest, dict), "run manifest must be an object")
    _require(manifest.get("code_sha") == policy_binding["code_sha"],
             "run manifest code_sha does not match policy_binding.code_sha")
    manifest_policy = manifest.get("policy_binding")
    _require(isinstance(manifest_policy, dict),
             "run manifest policy_binding is missing")
    for field in ("canon_hash", "law_semantic_hash",
                  "modality_matrix_sha256", "law_registry_sha256"):
        _require(manifest_policy.get(field) == policy_binding[field],
                 f"run manifest policy_binding.{field} does not match evaluation binding")

    results = json.loads(
        (run_dir / "calibration_results.json").read_text(encoding="utf-8"))
    meta = results.get("meta") if isinstance(results, dict) else None
    _require(isinstance(meta, dict), "calibration results meta is missing")
    _require(meta.get("code_sha") == policy_binding["code_sha"],
             "run results code_sha does not match policy_binding.code_sha")
    _require(meta.get("policy_binding") == manifest_policy,
             "run results policy_binding does not match the run manifest")


def evaluate_file(path: Path, alpha: float = DEFAULT_ALPHA,
                  run_dir: Path | None = None) -> dict:
    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    _require(isinstance(data, dict), "evaluation file must contain an object")
    cases = validate_cases(
        data.get("cases"), require_policy=True, payload=data,
        require_label_provenance=True)
    _require(run_dir is not None,
             "file-based evaluation requires a run directory to verify the bound artifacts")
    _verify_run_inputs(data["run_inputs"], run_dir, data["policy_binding"])
    report = evaluate(cases, alpha)
    report["policy_binding"] = data["policy_binding"]
    report["run_inputs"] = data["run_inputs"]
    report["run_inputs_verified"] = True
    report["input_file_sha256"] = _file_sha256(path)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("input", help="labeled cases JSON or split manifest JSON")
    parser.add_argument("--alpha", type=float, default=DEFAULT_ALPHA)
    parser.add_argument("--run-dir", default=None,
                        help="verify the three bound run-input file digests")
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

    report = evaluate_file(path, args.alpha, Path(args.run_dir) if args.run_dir else None)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
