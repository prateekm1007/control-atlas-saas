"""
WP-SCI-1 calibration runner tests (directive §4, §7).

Exercises the versioned runner on local fixtures (no network):
  - complete run: every target acquired+measured, digests pinned, manifest
    written and internally verifiable, exit code 0;
  - missing coordinate file: recorded as explicit acquisition failure,
    run_status incomplete, exit code non-zero (fail-closed completeness);
  - digest verification against a prior manifest: a changed coordinate file
    with a stale pinned digest is a hard acquisition failure.
"""
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from calibration.runner import run_calibration, write_manifest_sha256

REPO = Path(__file__).resolve().parents[2]
FIXTURE_PDB = REPO / "backend" / "tests" / "fixtures" / "real_pdb_1crn.pdb"

DATASET = {
    "smoke": [
        {"pdb_id": "1CRN", "resolution": 1.5, "desc": "crambin fixture"},
        {"pdb_id": "1UBQ", "resolution": 1.8, "desc": "ubiquitin fixture"},
    ]
}


@pytest.fixture()
def coord_dir(tmp_path):
    d = tmp_path / "coords"
    d.mkdir()
    shutil.copy(FIXTURE_PDB, d / "1CRN.pdb")
    # a second, minimal valid PDB so the dataset has two targets
    (d / "1UBQ.pdb").write_text(FIXTURE_PDB.read_text().replace(
        "1CRN", "1UBQ"))
    return d


def test_complete_local_run_binds_digests_and_manifest(tmp_path, coord_dir):
    out = tmp_path / "run"
    doc, code = run_calibration(DATASET, out, source_dir=coord_dir)

    assert code == 0
    assert doc["summary"]["run_status"] == "complete"
    assert doc["summary"]["acquired"] == 2
    assert doc["summary"]["measured"] == 2

    # coordinate digests pinned to the actual file bytes
    coord_manifest = json.loads((out / "coordinate_manifest.json").read_text())
    for pdb_id in ("1CRN", "1UBQ"):
        expected = hashlib.sha256((coord_dir / f"{pdb_id}.pdb").read_bytes()).hexdigest()
        assert coord_manifest[pdb_id]["sha256"] == expected
        assert coord_manifest[pdb_id]["acquisition_status"] == "acquired"

    # run manifest binds code and policy identity
    rm = json.loads((out / "run_manifest.json").read_text())
    assert rm["code_sha"], "code SHA must be recorded"
    assert rm["policy_binding"]["canon_hash"] == "6a9cd4b4349b81de"
    assert len(rm["policy_binding"]["law_semantic_hash"]) == 64
    assert len(rm["policy_binding"]["modality_matrix_sha256"]) == 64
    assert rm["run_type"] == "diagnostic"

    # MANIFEST.sha256 covers the run dir and verifies
    manifest_lines = (out / "MANIFEST.sha256").read_text().strip().splitlines()
    covered = {line.split("  ")[1] for line in manifest_lines}
    assert {"run_manifest.json", "coordinate_manifest.json",
            "calibration_results.json"} <= covered
    # every measurement carries full provenance fields
    t = doc["targets"]["1CRN"]
    for lid, row in t["laws"].items():
        assert {"observed", "threshold", "operator", "units", "status",
                "method", "sample_size"} <= set(row)
        assert row["status"] in ("PASS", "VETO", "FAIL", "NOT_MEASURED",
                                 "NOT_APPLICABLE", "MEASUREMENT_FAILED",
                                 "FAIL (Advisory)")


def test_missing_coordinate_is_explicit_failure_not_success(tmp_path, coord_dir):
    (coord_dir / "1UBQ.pdb").unlink()
    out = tmp_path / "run"
    doc, code = run_calibration(DATASET, out, source_dir=coord_dir)

    assert code == 1  # fail-closed completeness: non-zero exit
    assert doc["summary"]["run_status"] == "incomplete"
    assert doc["summary"]["acquisition_failures"] == 1
    failed = doc["targets"]["1UBQ"]
    assert failed["acquisition"]["acquisition_status"] == "failed"
    assert failed["measurement_status"] == "NOT_MEASURED"
    assert failed["verdict"] is None
    assert failed["laws"] is None
    assert any("1UBQ" in f for f in doc["summary"]["failures"])


def test_stale_digest_is_a_hard_failure(tmp_path, coord_dir):
    prior_manifest = tmp_path / "prior.json"
    prior_manifest.write_text(json.dumps({
        "1CRN": {"sha256": "0" * 64},
        "1UBQ": {"sha256": "1" * 64},
    }))
    out = tmp_path / "run"
    doc, code = run_calibration(DATASET, out, source_dir=coord_dir,
                                prior_digests=json.loads(prior_manifest.read_text()))
    assert code == 1
    assert doc["summary"]["acquisition_failures"] == 2
    for pdb_id in ("1CRN", "1UBQ"):
        assert "digest mismatch" in doc["targets"][pdb_id]["acquisition"]["detail"]


def test_manifest_regeneration_is_idempotent(tmp_path, coord_dir):
    out = tmp_path / "run"
    run_calibration(DATASET, out, source_dir=coord_dir)
    first = (out / "MANIFEST.sha256").read_text()
    write_manifest_sha256(out)
    assert (out / "MANIFEST.sha256").read_text() == first


def test_law120_contained_flag_flows_through_runner(tmp_path, coord_dir):
    out = tmp_path / "run"
    doc, _ = run_calibration(DATASET, out, source_dir=coord_dir)
    row = doc["targets"]["1CRN"]["laws"]["LAW-120"]
    assert row.get("contained") is True
    assert row["containment_reason"] == "measurement-definition-repair"
