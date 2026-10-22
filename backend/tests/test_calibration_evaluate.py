"""Calibration evaluator and split-manifest validation tests."""
import json

import pytest

from calibration.evaluate import (
    EvaluationInputError,
    evaluate,
    evaluate_file,
    validate_cases,
    validate_split_manifest,
    wilson_interval,
)


def _case(sid, label, verdict, **overrides):
    case = {
        "structure_id": sid,
        "label": label,
        "strata": {"modality": "x_ray", "resolution_band": "medium"},
        "verdict": verdict,
        "blocking_laws": [],
        "indeterminate_causes": [],
        "not_applicable": [],
        "not_measured": [],
    }
    case.update(overrides)
    return case


def _cases():
    return [
        _case("TP1", "defective", "VETO", blocking_laws=["LAW-120"]),
        _case("FN1", "defective", "PASS"),
        _case("ID1", "defective", "INDETERMINATE",
              indeterminate_causes=["LAW-120"], not_measured=["LAW-120"]),
        _case("TN1", "acceptable", "PASS"),
        _case("FP1", "acceptable", "VETO", blocking_laws=["LAW-130"]),
        _case("ID2", "acceptable", "INDETERMINATE",
              indeterminate_causes=["LAW-130"], not_measured=["LAW-130"]),
    ]


def test_wilson_interval_handles_empty_denominator():
    result = wilson_interval(0, 0)
    assert result["estimate"] is None
    assert result["total"] == 0
    assert "not computable" in result["note"]


def test_wilson_interval_uses_two_sided_95_percent_z():
    result = wilson_interval(70, 100)
    assert result["lower"] == 0.6042
    assert result["upper"] == 0.7811


def test_indeterminate_is_excluded_from_sensitivity_and_specificity():
    report = evaluate(_cases())
    cm = report["confusion_matrix"]
    assert (cm["tp"], cm["fn"], cm["fp"], cm["tn"]) == (1, 1, 1, 1)
    assert cm["indeterminate_defective"] == 1
    assert cm["indeterminate_acceptable"] == 1
    assert report["sensitivity"]["total"] == 2
    assert report["specificity"]["total"] == 2
    assert report["sensitivity"]["estimate"] == 0.5
    assert report["specificity"]["estimate"] == 0.5


def test_per_law_denominators_are_their_own_case_sets():
    report = evaluate(_cases())
    per_law = report["per_law_contribution"]
    assert per_law["denominators"] == {
        "veto_cases": 2,
        "indeterminate_cases": 2,
        "all_cases": 6,
    }
    assert per_law["veto_causes"]["LAW-120"]["share_of_veto_cases"]["total"] == 2
    assert per_law["veto_causes"]["LAW-120"]["count"] == 1
    assert per_law["indeterminate_causes"]["LAW-130"]["share_of_indeterminate_cases"]["total"] == 2


def test_missing_rate_uses_eligible_law_opportunities():
    report = evaluate(_cases())
    missing = report["missing_measurement_rate"]
    assert missing["cases_with_any_missing"]["successes"] == 2
    assert missing["cases_with_any_missing"]["total"] == 6
    assert missing["per_law"]["LAW-120"]["not_measured"] == 1
    assert missing["per_law"]["LAW-120"]["eligible"] == 6
    assert missing["law_opportunities"]["missing"] == 2


def test_strata_do_not_recurse_and_report_n():
    report = evaluate(_cases())
    group = report["strata"]["modality"]["groups"]["x_ray"]
    assert group["n"] == 6
    assert "per_law_contribution" not in group
    assert report["strata"]["resolution_band"]["unstratified"] == 0


def test_validation_rejects_duplicate_ids():
    cases = _cases()
    cases.append(_case("TP1", "acceptable", "PASS"))
    with pytest.raises(EvaluationInputError, match="duplicate structure_id"):
        validate_cases(cases)


def test_validation_rejects_retired_or_heuristic_veto_cause():
    cases = [_case("BAD", "defective", "VETO", blocking_laws=["LAW-155"])]
    with pytest.raises(EvaluationInputError, match="unknown law"):
        validate_cases(cases)
    cases = [_case("BAD2", "defective", "VETO", blocking_laws=["LAW-182"])]
    with pytest.raises(EvaluationInputError, match="not a deterministic law"):
        validate_cases(cases)


def test_validation_rejects_inconsistent_cause_lists():
    cases = [_case("BAD", "defective", "INDETERMINATE",
                   indeterminate_causes=["LAW-120"], not_measured=[])]
    with pytest.raises(EvaluationInputError, match="also be listed as NOT_MEASURED"):
        validate_cases(cases)
    cases = [_case("BAD2", "defective", "PASS", blocking_laws=["LAW-120"])]
    with pytest.raises(EvaluationInputError, match="only valid for VETO"):
        validate_cases(cases)


def test_validation_rejects_missing_or_malformed_policy_binding():
    payload = {
        "policy_binding": {"code_sha": "abc", "law_semantic_hash": "0" * 64},
        "cases": _cases(),
    }
    with pytest.raises(EvaluationInputError, match="code_sha"):
        validate_cases(payload["cases"], require_policy=True, payload=payload)
    payload["policy_binding"] = {
        "code_sha": "a" * 40,
        "law_semantic_hash": "not-a-hash",
    }
    with pytest.raises(EvaluationInputError, match="law_semantic_hash"):
        validate_cases(payload["cases"], require_policy=True, payload=payload)


def _split_entry(sid, cohort="A", label="acceptable"):
    entry = {
        "structure_id": sid,
        "label_source": "wwpdb-validation" if cohort == "A" else "perturbed-by-construction",
        "label": label,
        "cohort": cohort,
        "strata": {"modality": "x_ray", "resolution_band": "medium"},
        "coordinate_sha256": "a" * 64,
    }
    if cohort == "B1":
        entry["defect_manifest_ref"] = f"defects/{sid}.json"
    return entry


def _manifest():
    return {
        "schema_version": "split-manifest/1.0.0",
        "protocol": "backend/calibration/protocol/evaluation_protocol_v1.md",
        "created_utc": "2026-10-21T00:00:00Z",
        "policy_binding": {"code_sha": "a" * 40, "law_semantic_hash": "b" * 64},
        "splits": {
            "calibration": [_split_entry("1AAA")],
            "validation": [_split_entry("1BBB")],
            "held_out": [_split_entry("2CCC", cohort="B1", label="defective")],
        },
        "overlap_checks": {
            "identity_overlap_found": False,
            "near_identity_pairs": [],
            "near_identity_threshold": 0.9,
            "checked_by": "test",
            "all_checks_passed": True,
        },
    }


def test_split_manifest_validates_and_counts():
    summary = validate_split_manifest(_manifest())
    assert summary["split_counts"] == {
        "calibration": 1,
        "validation": 1,
        "held_out": 1,
    }
    assert summary["total_entries"] == 3


def test_split_manifest_requires_policy_binding_and_cohort_labels():
    manifest = _manifest()
    manifest.pop("policy_binding")
    with pytest.raises(EvaluationInputError, match="policy_binding"):
        validate_split_manifest(manifest)

    manifest = _manifest()
    manifest["splits"]["held_out"][0]["label_source"] = "wwpdb-validation"
    with pytest.raises(EvaluationInputError, match="requires label_source"):
        validate_split_manifest(manifest)

    manifest = _manifest()
    manifest["splits"]["held_out"][0]["label"] = "acceptable"
    with pytest.raises(EvaluationInputError, match="requires label"):
        validate_split_manifest(manifest)

    manifest = _manifest()
    manifest["overlap_checks"]["near_identity_threshold"] = True
    with pytest.raises(EvaluationInputError, match="proportion"):
        validate_split_manifest(manifest)


def test_split_manifest_rejects_cross_split_identity_and_failed_checks():
    manifest = _manifest()
    manifest["splits"]["held_out"][0]["structure_id"] = "1AAA"
    with pytest.raises(EvaluationInputError, match="more than one split"):
        validate_split_manifest(manifest)

    manifest = _manifest()
    manifest["overlap_checks"]["all_checks_passed"] = False
    with pytest.raises(EvaluationInputError, match="all overlap checks"):
        validate_split_manifest(manifest)

    manifest = _manifest()
    manifest["splits"]["held_out"][0].pop("defect_manifest_ref")
    with pytest.raises(EvaluationInputError, match="defect_manifest_ref"):
        validate_split_manifest(manifest)


def test_evaluate_file_binds_policy_and_input_digest(tmp_path):
    payload = {
        "policy_binding": {"code_sha": "a" * 40, "law_semantic_hash": "b" * 64},
        "cases": _cases(),
    }
    path = tmp_path / "cases.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    report = evaluate_file(path)
    assert report["policy_binding"] == payload["policy_binding"]
    assert len(report["input_file_sha256"]) == 64
    assert report["n_cases"] == 6
