"""Calibration evaluator and split-manifest validation tests."""
import hashlib
import json
from pathlib import Path

import pytest

from calibration.evaluate import (
    EvaluationInputError,
    _law_semantic_hash,
    _pair_results,
    _sequence_sha256,
    evaluate,
    evaluate_file,
    validate_cases,
    validate_split_manifest,
    wilson_interval,
)
from tos.governance.modality_matrix import compute_matrix_hash
from tos.governance.station_sop import LAW_CANON_HASH


def _case(sid, label, verdict, **overrides):
    cohort = "A" if label == "acceptable" else "B1"
    case = {
        "structure_id": sid,
        "label": label,
        "label_source": "wwpdb-validation" if cohort == "A" else "perturbed-by-construction",
        "label_evidence_ref": f"labels/{sid}.json",
        "label_evidence_sha256": "e" * 64,
        "cohort": cohort,
        "strata": {"modality": "x_ray", "resolution_band": "medium"},
        "verdict": verdict,
        "blocking_laws": [],
        "indeterminate_causes": [],
        "not_applicable": [],
        "not_measured": [],
        "measurement_failed": [],
    }
    if cohort == "B1":
        case["defect_manifest_ref"] = f"defects/{sid}.json"
        case["defect_manifest_sha256"] = "d" * 64
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
    with pytest.raises(EvaluationInputError, match="not a deterministic law"):
        validate_cases(cases)
    cases = [_case("BAD2", "defective", "VETO", blocking_laws=["LAW-182"])]
    with pytest.raises(EvaluationInputError, match="not a deterministic law"):
        validate_cases(cases)


def test_validation_rejects_inconsistent_cause_lists():
    cases = [_case("BAD", "defective", "INDETERMINATE",
                   indeterminate_causes=["LAW-120"], not_measured=[])]
    with pytest.raises(EvaluationInputError, match="listed as missing evidence"):
        validate_cases(cases)
    cases = [_case("BAD2", "defective", "PASS", blocking_laws=["LAW-120"])]
    with pytest.raises(EvaluationInputError, match="only valid for VETO"):
        validate_cases(cases)


def test_measurement_failed_is_missing_not_pass_or_measured():
    cases = [_case(
        "FAIL1", "defective", "INDETERMINATE",
        indeterminate_causes=["LAW-120"], measurement_failed=["LAW-120"])]
    report = evaluate(cases)
    assert report["confusion_matrix"]["indeterminate_defective"] == 1
    missing = report["missing_measurement_rate"]
    assert missing["cases_with_any_missing"]["successes"] == 1
    assert missing["per_law"]["LAW-120"]["measurement_failed"] == 1
    assert missing["per_law"]["LAW-120"]["not_measured"] == 0
    assert missing["law_opportunities"]["measurement_failed"] == 1
    assert report["per_law_contribution"][
        "measurement_failed_causes"]["LAW-120"]["count"] == 1


def test_validation_rejects_contradictory_measurement_states():
    cases = [_case("BAD", "defective", "PASS", measurement_failed=["LAW-120"])]
    with pytest.raises(EvaluationInputError, match="PASS is invalid"):
        validate_cases(cases)
    cases = [_case(
        "BAD2", "defective", "INDETERMINATE",
        indeterminate_causes=["LAW-120"],
        not_measured=["LAW-120"], measurement_failed=["LAW-120"])]
    with pytest.raises(EvaluationInputError, match="both NOT_MEASURED"):
        validate_cases(cases)
    cases = [_case(
        "BAD3", "defective", "INDETERMINATE",
        indeterminate_causes=["LAW-120"])]
    with pytest.raises(EvaluationInputError, match="listed as missing evidence"):
        validate_cases(cases)


def _policy_binding():
    registry = Path(__file__).resolve().parents[2] / "governance" / "law_registry.yaml"
    return {
        "code_sha": "a" * 40,
        "canon_hash": LAW_CANON_HASH,
        "law_semantic_hash": _law_semantic_hash(),
        "modality_matrix_sha256": compute_matrix_hash(),
        "law_registry_sha256": hashlib.sha256(registry.read_bytes()).hexdigest(),
    }


def _run_inputs():
    return {
        "run_manifest_sha256": "1" * 64,
        "coordinate_manifest_sha256": "2" * 64,
        "calibration_results_sha256": "3" * 64,
        "run_integrity_manifest_sha256": "4" * 64,
    }


def test_validation_rejects_missing_or_malformed_policy_binding():
    payload = {"policy_binding": {"code_sha": "abc"}, "cases": _cases()}
    with pytest.raises(EvaluationInputError, match="code_sha"):
        validate_cases(payload["cases"], require_policy=True, payload=payload)
    payload = {
        "policy_binding": {**_policy_binding(), "law_semantic_hash": "not-a-hash"},
        "run_inputs": _run_inputs(),
        "cases": _cases(),
    }
    with pytest.raises(EvaluationInputError, match="law_semantic_hash"):
        validate_cases(payload["cases"], require_policy=True, payload=payload)


def test_validation_rejects_policy_identity_drift_and_missing_run_inputs():
    payload = {
        "policy_binding": {**_policy_binding(), "canon_hash": "0" * 16},
        "run_inputs": _run_inputs(),
        "cases": _cases(),
    }
    with pytest.raises(EvaluationInputError, match="canon_hash"):
        validate_cases(payload["cases"], require_policy=True, payload=payload)
    payload["policy_binding"] = _policy_binding()
    payload["run_inputs"] = {"run_manifest_sha256": "1" * 64}
    with pytest.raises(EvaluationInputError, match="run_inputs"):
        validate_cases(payload["cases"], require_policy=True, payload=payload)


def test_label_provenance_is_required_for_evidence_reports():
    cases = _cases()
    cases[0].pop("label_source")
    with pytest.raises(EvaluationInputError, match="label_source"):
        validate_cases(cases, require_label_provenance=True)


def _split_entry(sid, cohort="A", label="acceptable", sequence=None, family=None):
    sequence = sequence or "ACDEFGHIKLMNPQRSTVWY"
    family_id = family or f"family-{sid}"
    entry = {
        "structure_id": sid,
        "label_source": "wwpdb-validation" if cohort == "A" else "perturbed-by-construction",
        "label_evidence_ref": f"labels/{sid}.json",
        "label_evidence_sha256": "e" * 64,
        "label": label,
        "cohort": cohort,
        "strata": {"modality": "x_ray", "resolution_band": "medium"},
        "coordinate_sha256": hashlib.sha256(sid.encode()).hexdigest(),
        "sequence": sequence,
        "sequence_sha256": _sequence_sha256(sequence),
        "sequence_source_ref": f"source/{sid}.fasta",
        "sequence_source_sha256": "f" * 64,
        "family_id": family_id,
        "family_source_ref": "source/families.json",
        "family_source_sha256": "f" * 64,
    }
    if cohort == "B1":
        entry["defect_manifest_ref"] = f"defects/{sid}.json"
        entry["defect_manifest_sha256"] = "d" * 64
    return entry


def _manifest(sequences=None, near_pairs=None):
    sequences = sequences or {
        "1AAA": "ACDEFGHIKLMNPQRSTVWY",
        "1BBB": "YYYYYYYYYYYYYYYYYYYY",
        "2CCC": "WWWWWWWWWWWWWWWWWWWW",
    }
    manifest = {
        "schema_version": "split-manifest/1.1.0",
        "protocol": "backend/calibration/protocol/evaluation_protocol_v1.md",
        "created_utc": "2026-10-21T00:00:00Z",
        "historical_exclusions": ["7OLD"],
        "historical_exclusion_provenance": {
            "inventory_ref": "governance/historical_exclusions.json",
            "inventory_sha256": "c" * 64,
            "excluded_structure_ids": ["7OLD"],
        },
        "policy_binding": {"code_sha": "a" * 40, "law_semantic_hash": "b" * 64},
        "splits": {
            "calibration": [_split_entry("1AAA", sequence=sequences["1AAA"])],
            "validation": [_split_entry("1BBB", sequence=sequences["1BBB"])],
            "held_out": [_split_entry(
                "2CCC", cohort="B1", label="defective",
                sequence=sequences["2CCC"])],
        },
    }
    entries = {
        name: [{**entry, "_sequence": entry["sequence"].upper()}
               for entry in members]
        for name, members in manifest["splits"].items()
    }
    pair_digest, result_digest, computed_near = _pair_results(entries, 0.9)
    manifest["overlap_checks"] = {
        "sequence_identity": {
            "method": "normalized-global-edit-identity",
            "tool": "calibration.evaluate",
            "tool_version": "1.0.0",
            "near_identity_threshold": 0.9,
            "checked_pairs": 3,
            "checked_pair_digest": pair_digest,
            "result_digest": result_digest,
        },
        "identity_overlap_found": bool(computed_near),
        "family_overlap_found": False,
        "near_identity_pairs": near_pairs if near_pairs is not None else computed_near,
        "all_checks_passed": True,
    }
    return manifest


def test_split_manifest_validates_and_counts():
    summary = validate_split_manifest(_manifest())
    assert summary["split_counts"] == {
        "calibration": 1,
        "validation": 1,
        "held_out": 1,
    }
    assert summary["total_entries"] == 3
    assert summary["checked_pairs"] == 3
    assert summary["near_identity_pair_count"] == 0


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
    manifest["overlap_checks"]["sequence_identity"]["near_identity_threshold"] = True
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


def test_split_manifest_recomputes_leakage_instead_of_trusting_flags():
    manifest = _manifest(sequences={
        "1AAA": "ACDEFGHIKLMNPQRSTVWY",
        "1BBB": "ACDEFGHIKLMNPQRSTVWA",
        "2CCC": "WWWWWWWWWWWWWWWWWWWW",
    })
    manifest["overlap_checks"]["near_identity_pairs"] = []
    manifest["overlap_checks"]["identity_overlap_found"] = False
    with pytest.raises(EvaluationInputError, match="near_identity_pairs"):
        validate_split_manifest(manifest)

    manifest = _manifest()
    manifest["overlap_checks"]["sequence_identity"]["checked_pairs"] = 99
    with pytest.raises(EvaluationInputError, match="checked_pairs"):
        validate_split_manifest(manifest)

    manifest = _manifest()
    manifest["overlap_checks"]["sequence_identity"]["result_digest"] = "0" * 64
    with pytest.raises(EvaluationInputError, match="result_digest"):
        validate_split_manifest(manifest)

    manifest = _manifest(sequences={
        "1AAA": "ACDEFGHIKLMNPQRSTVWY",
        "1BBB": "ACDEFGHIKLMNPQRSTVWA",
        "2CCC": "WWWWWWWWWWWWWWWWWWWW",
    })
    with pytest.raises(EvaluationInputError, match="near-identity pairs"):
        validate_split_manifest(manifest)


def test_split_manifest_rejects_repeated_digests_families_and_exclusions():
    manifest = _manifest()
    manifest["splits"]["held_out"][0]["coordinate_sha256"] = \
        manifest["splits"]["calibration"][0]["coordinate_sha256"]
    with pytest.raises(EvaluationInputError, match="repeated coordinate digest"):
        validate_split_manifest(manifest)

    manifest = _manifest()
    manifest["splits"]["held_out"][0]["family_id"] = \
        manifest["splits"]["calibration"][0]["family_id"]
    with pytest.raises(EvaluationInputError, match="family overlap"):
        validate_split_manifest(manifest)

    manifest = _manifest()
    manifest["splits"]["held_out"][0]["structure_id"] = "7OLD"
    with pytest.raises(EvaluationInputError, match="historical exclusions"):
        validate_split_manifest(manifest)

    manifest = _manifest()
    manifest["splits"]["held_out"][0]["sequence"] = "ACDEFGHIKLMNPQRSTVWY"
    manifest["splits"]["held_out"][0]["sequence_sha256"] = _sequence_sha256(
        "ACDEFGHIKLMNPQRSTVWY")
    with pytest.raises(EvaluationInputError, match="repeated sequence"):
        validate_split_manifest(manifest)


def test_evaluate_file_requires_bound_run_inputs(tmp_path):
    payload = {
        "policy_binding": _policy_binding(),
        "run_inputs": _run_inputs(),
        "cases": _cases(),
    }
    path = tmp_path / "cases.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(EvaluationInputError, match="run directory"):
        evaluate_file(path)


def test_evaluate_file_verifies_bound_run_inputs(tmp_path):
    payload = {
        "policy_binding": _policy_binding(),
        "run_inputs": _run_inputs(),
        "cases": _cases(),
    }
    path = tmp_path / "cases.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    with pytest.raises(EvaluationInputError, match="run input is missing"):
        evaluate_file(path, run_dir=run_dir)

    policy = payload["policy_binding"]
    files = {
        "run_manifest.json": {
            "code_sha": policy["code_sha"],
            "policy_binding": policy,
            "run_type": "authoritative",
        },
        "coordinate_manifest.json": {"1ABC": {"status": "acquired"}},
        "run_integrity_manifest.json": {
            "run_status": "complete",
            "run_type": "authoritative",
            "expected_target_count": 1,
            "acquired_count": 1,
            "measured_count": 1,
            "policy_binding": policy,
            "code_sha": policy["code_sha"],
            "artifact_digests": {},
        },
        "calibration_results.json": {
            "meta": {
                "code_sha": policy["code_sha"],
                "policy_binding": policy,
            },
            "summary": {
                "measured": 1,
                "expected_targets": 1,
            },
        },
    }
    files["run_integrity_manifest.json"]["artifact_digests"] = {
        "run_manifest.json": hashlib.sha256(json.dumps(files["run_manifest.json"], sort_keys=True).encode()).hexdigest(),
        "coordinate_manifest.json": hashlib.sha256(json.dumps(files["coordinate_manifest.json"], sort_keys=True).encode()).hexdigest(),
        "calibration_results.json": hashlib.sha256(json.dumps(files["calibration_results.json"], sort_keys=True).encode()).hexdigest(),
    }
    for field, filename in {
        "run_manifest_sha256": "run_manifest.json",
        "coordinate_manifest_sha256": "coordinate_manifest.json",
        "calibration_results_sha256": "calibration_results.json",
    }.items():
        content = json.dumps(files[filename], sort_keys=True).encode()
        (run_dir / filename).write_bytes(content)
        payload["run_inputs"][field] = hashlib.sha256(content).hexdigest()
    integrity_content = json.dumps(files["run_integrity_manifest.json"], sort_keys=True).encode()
    (run_dir / "run_integrity_manifest.json").write_bytes(integrity_content)
    payload["run_inputs"]["run_integrity_manifest_sha256"] = hashlib.sha256(integrity_content).hexdigest()
    manifest_lines = []
    for fname in ["run_manifest.json", "coordinate_manifest.json", "calibration_results.json", "run_integrity_manifest.json"]:
        fpath = run_dir / fname
        manifest_lines.append(f"{hashlib.sha256(fpath.read_bytes()).hexdigest()}  {fname}")
    (run_dir / "MANIFEST.sha256").write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")
    path.write_text(json.dumps(payload), encoding="utf-8")
    report = evaluate_file(path, run_dir=run_dir)
    assert report["run_inputs_verified"] is True

    files["run_manifest.json"]["code_sha"] = "b" * 40
    content = json.dumps(files["run_manifest.json"], sort_keys=True).encode()
    (run_dir / "run_manifest.json").write_bytes(content)
    payload["run_inputs"]["run_manifest_sha256"] = hashlib.sha256(content).hexdigest()
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(EvaluationInputError, match="integrity manifest digest mismatch"):
        evaluate_file(path, run_dir=run_dir)
