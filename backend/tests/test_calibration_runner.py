"""Calibration runner evidence-boundary tests. No network access."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from calibration import runner


FIXTURE = Path(__file__).parent / "fixtures" / "real_pdb_1crn.pdb"


def _pdb_bytes() -> bytes:
    return FIXTURE.read_bytes()


def _response(url, content, status=200, final_url=None):
    return SimpleNamespace(
        status_code=status,
        content=content,
        url=final_url or url,
    )


def _clean_provenance():
    return {
        "code_sha": "a" * 40,
        "code_sha_source": "test",
        "working_tree_clean": True,
        "tree_status": [],
    }


def test_measure_target_uses_live_engine_and_counts_deterministic_laws():
    measured = runner.measure_target(_pdb_bytes(), "experimental")
    assert measured["verdict"] in {"PASS", "VETO", "INDETERMINATE"}
    assert measured["det_total"] == runner.DETERMINISTIC_COUNT
    assert measured["coverage"] == 100.0
    assert measured["atoms"] > 0
    assert measured["is_experimental"] is True
    assert set(measured["laws"]) == set(runner.LAW_CANON)
    assert all(
        row["method"] in {
            "deterministic", "advisory", "advisory_experimental",
            "heuristic", "retired",
        }
        for row in measured["laws"].values()
    )


def test_measure_target_forces_predicted_source_regime():
    measured = runner.measure_target(_pdb_bytes(), "predicted")
    assert measured["is_experimental"] is False
    assert measured["method_type"] == "predicted"
    assert measured["coverage"] < 100.0


def test_measure_target_rejects_json_and_empty_atoms():
    with pytest.raises(ValueError, match="JSON payload"):
        runner.measure_target(b'{"not": "a structure"}', "predicted")
    with pytest.raises(ValueError, match="no atoms parsed"):
        runner.measure_target(
            b"HEADER    TEST\n" + b"REMARK   1 FILLER\n" * 40, "experimental")


def test_run_calibration_local_end_to_end(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "_engine_provenance", _clean_provenance)
    source = tmp_path / "structures"
    source.mkdir()
    data = _pdb_bytes()
    (source / "1CRN.pdb").write_bytes(data)
    (source / "1ABC.pdb").write_bytes(data)

    dataset = {
        "xray_high": [
            {"pdb_id": "1CRN", "resolution": 1.5, "desc": "test"},
            {"pdb_id": "1ABC", "resolution": 1.5, "desc": "test"},
        ]
    }
    out = tmp_path / "run"
    document, code = runner.run_calibration(
        dataset, out, source_dir=source,
        prior_digests={
            "1CRN": {"sha256": hashlib.sha256(data).hexdigest()},
            "1ABC": {"sha256": hashlib.sha256(data).hexdigest()},
        },
    )

    assert code == 0
    assert document["summary"]["run_status"] == "complete"
    assert document["summary"]["measured"] == 2
    assert document["meta"]["working_tree_clean"] is True
    assert document["meta"]["policy_binding"]["canon_hash"] == runner.LAW_CANON_HASH
    assert len(document["meta"]["policy_binding"]["law_semantic_hash"]) == 64
    assert len(document["meta"]["policy_binding"]["modality_matrix_sha256"]) == 64

    coordinate_manifest = json.loads(
        (out / "coordinate_manifest.json").read_text(encoding="utf-8"))
    assert coordinate_manifest["1CRN"]["sha256"] == hashlib.sha256(data).hexdigest()
    assert coordinate_manifest["1ABC"]["sha256"] == hashlib.sha256(data).hexdigest()

    results = json.loads(
        (out / "calibration_results.json").read_text(encoding="utf-8"))
    assert results["targets"]["1CRN"]["measurement_status"] == "measured"
    assert results["targets"]["1CRN"]["laws"]["LAW-155"]["status"] == "RETIRED"

    manifest_lines = (out / "MANIFEST.sha256").read_text(encoding="utf-8").splitlines()
    manifest_names = {line.split("  ", 1)[1] for line in manifest_lines}
    assert manifest_names == {
        "run_manifest.json",
        "coordinate_manifest.json",
        "calibration_results.json",
    }
    for line in manifest_lines:
        digest, name = line.split("  ", 1)
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == digest


def test_run_refuses_dirty_tree_before_creating_output(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "_engine_provenance", lambda: {
        "code_sha": "b" * 40,
        "code_sha_source": "test",
        "working_tree_clean": False,
        "tree_status": [" M backend/tos/engine/tier1_measurements.py"],
    })
    out = tmp_path / "run"
    with pytest.raises(runner.CalibrationRunError, match="dirty engine tree"):
        runner.run_calibration({"xray_high": [{"pdb_id": "1CRN"}]}, out)
    assert not out.exists()


def test_run_refuses_existing_output_and_duplicate_ids(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "_engine_provenance", _clean_provenance)
    out = tmp_path / "run"
    out.mkdir()
    with pytest.raises(runner.CalibrationRunError, match="overwrite"):
        runner.run_calibration({"xray_high": [{"pdb_id": "1CRN"}]}, out)

    out2 = tmp_path / "run2"
    with pytest.raises(runner.CalibrationRunError, match="duplicate pdb_id"):
        runner.run_calibration({
            "a": [{"pdb_id": "1CRN"}],
            "b": [{"pdb_id": "1CRN"}],
        }, out2)


def test_missing_target_and_digest_mismatch_are_incomplete(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "_engine_provenance", _clean_provenance)
    source = tmp_path / "structures"
    source.mkdir()
    data = _pdb_bytes()
    (source / "1CRN.pdb").write_bytes(data)
    (source / "1ABC.pdb").write_bytes(data)

    document, code = runner.run_calibration(
        {"xray": [{"pdb_id": "1CRN"}, {"pdb_id": "9ZZZ"}]},
        tmp_path / "missing", source_dir=source,
    )
    assert code == 1
    assert document["summary"]["acquired"] == 1
    assert document["summary"]["measured"] == 1
    assert document["targets"]["9ZZZ"]["measurement_status"] == "NOT_MEASURED"

    document, code = runner.run_calibration(
        {"xray": [{"pdb_id": "1CRN"}]},
        tmp_path / "bad-digest", source_dir=source,
        prior_digests={"1CRN": {"sha256": "0" * 64}},
    )
    assert code == 1
    assert "digest mismatch" in document["targets"]["1CRN"]["detail"]

    document, code = runner.run_calibration(
        {"xray": [{"pdb_id": "1CRN"}]},
        tmp_path / "bad-pin", source_dir=source,
        prior_digests={"1CRN": {"sha256": "not-a-digest"}},
    )
    assert code == 1
    assert "invalid prior digest" in document["targets"]["1CRN"]["detail"]

    document, code = runner.run_calibration(
        {"xray": [{"pdb_id": "1CRN"}, {"pdb_id": "1ABC"}]},
        tmp_path / "partial-pin", source_dir=source,
        prior_digests={"1CRN": {"sha256": hashlib.sha256(data).hexdigest()}},
    )
    assert code == 1
    assert document["targets"]["1CRN"]["measurement_status"] == "measured"
    assert document["targets"]["1ABC"]["measurement_status"] == "NOT_MEASURED"
    assert "invalid prior digest" in document["targets"]["1ABC"]["detail"]


def _afdb_entry(version, entity="AF-Q5VSL9-F1", **overrides):
    entry = {
        "modelEntityId": entity,
        "latestVersion": version,
        "uniprotAccession": "Q5VSL9",
        "globalMetricValue": 88.1,
        "pdbUrl": f"https://alphafold.ebi.ac.uk/files/{entity}-model_v{version}.pdb",
        "cifUrl": f"https://alphafold.ebi.ac.uk/files/{entity}-model_v{version}.cif",
    }
    entry.update(overrides)
    return entry


def test_afdb_selects_highest_exact_entity_version():
    metadata = json.dumps([_afdb_entry(3), _afdb_entry(5)]).encode()
    fetched = []

    def fake_get(url, timeout=None):
        fetched.append(url)
        if "api/prediction" in url:
            return _response(url, metadata)
        return _response(url, _pdb_bytes())

    record, content = runner.acquire_coordinates(
        "AF-Q5VSL9-F1", None, None, fetcher=fake_get)

    assert record["acquisition_status"] == "acquired"
    assert record["selected_entity_id"] == "AF-Q5VSL9-F1"
    assert record["latest_version"] == 5
    assert record["source"].endswith("model_v5.pdb")
    assert record["format"] == "pdb"
    assert record["sha256"] == hashlib.sha256(_pdb_bytes()).hexdigest()
    assert content == _pdb_bytes()
    assert fetched[0] == "https://alphafold.ebi.ac.uk/api/prediction/Q5VSL9"


def test_afdb_fails_closed_on_entity_migration_and_404():
    metadata = json.dumps([_afdb_entry(4, entity="AF-Q5VSL9-F2")]).encode()
    record, content = runner.acquire_coordinates(
        "AF-Q5VSL9-F1", None, None,
        fetcher=lambda url, timeout=None: _response(url, metadata),
    )
    assert content is None
    assert record["acquisition_status"] == "failed"
    assert record["migrated_identifier"] == "AF-Q5VSL9-F2"
    assert "not served" in record["detail"]

    record, content = runner.acquire_coordinates(
        "AF-P62988-F1", None, None,
        fetcher=lambda url, timeout=None: _response(url, b"not found", status=404),
    )
    assert content is None
    assert "HTTP 404" in record["detail"]


def test_afdb_digest_pin_is_enforced():
    metadata = json.dumps([_afdb_entry(4)]).encode()

    def fake_get(url, timeout=None):
        if "api/prediction" in url:
            return _response(url, metadata)
        return _response(url, _pdb_bytes())

    record, content = runner.acquire_coordinates(
        "AF-Q5VSL9-F1", None,
        {"sha256": "0" * 64},
        fetcher=fake_get,
    )
    assert content is None
    assert record["acquisition_status"] == "failed"
    assert "digest mismatch" in record["detail"]


def test_afdb_fails_closed_on_cif_only_and_non_https_source():
    cif = b"data_TEST\n#\n" + b"_atom_site.group_PDB ATOM\n" * 20
    metadata = json.dumps([{
        "modelEntityId": "AF-Q5VSL9-F1",
        "latestVersion": 2,
        "uniprotAccession": "Q5VSL9",
        "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-Q5VSL9-F1-model_v2.cif",
    }]).encode()

    def fake_get(url, timeout=None):
        if "api/prediction" in url:
            return _response(url, metadata)
        return _response(url, cif)

    record, content = runner.acquire_coordinates(
        "AF-Q5VSL9-F1", None, None, fetcher=fake_get)
    assert content is None
    assert "mmCIF" in record["detail"]

    metadata = json.dumps([{
        "modelEntityId": "AF-Q5VSL9-F1",
        "latestVersion": 2,
        "uniprotAccession": "Q5VSL9",
        "pdbUrl": "http://alphafold.ebi.ac.uk/files/insecure.pdb",
    }]).encode()

    def fake_get_insecure(url, timeout=None):
        if "api/prediction" in url:
            return _response(url, metadata)
        return _response(url, _pdb_bytes())

    record, content = runner.acquire_coordinates(
        "AF-Q5VSL9-F1", None, None, fetcher=fake_get_insecure)
    assert content is None
    assert "non-HTTPS" in record["detail"]


def test_http_redirect_to_disallowed_host_is_rejected():
    def fake_get(url, timeout=None):
        return _response(
            url, _pdb_bytes(), final_url="https://example.invalid/evil.pdb")

    with pytest.raises(runner.AcquisitionError, match="host not allowed"):
        runner._download("https://files.rcsb.org/download/1CRN.pdb", fetcher=fake_get)
