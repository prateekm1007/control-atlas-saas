"""
P0 item 15 (2026-10-19 audit): calibration-runner hardening.

  - a run directory that already exists is never overwritten
  - the engine provenance records the commit of the IMPORTED tos package
    and the working-tree state; the CLI refuses dirty trees
  - run_status is 'complete' only when measured == expected and zero
    failures
  - the input filename follows the CONTENT-detected format
"""
import json
from pathlib import Path

import pytest

from calibration import runner as runner_mod
from calibration.runner import main, run_calibration

DATASET = {"smoke": [{"pdb_id": "1CRN", "resolution": 1.5, "desc": "x"}]}


def _coord_dir(tmp_path):
    d = tmp_path / "coords"
    d.mkdir()
    src = Path(__file__).parent / "fixtures" / "real_pdb_1crn.pdb"
    (d / "1CRN.pdb").write_bytes(src.read_bytes())
    return d


def test_existing_run_directory_is_never_overwritten(tmp_path):
    out = tmp_path / "run"
    out.mkdir()
    (out / "calibration_results.json").write_text("sentinel")
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        run_calibration(DATASET, out, source_dir=_coord_dir(tmp_path))
    assert (out / "calibration_results.json").read_text() == "sentinel"


def test_provenance_records_engine_commit_and_tree_state(tmp_path):
    out = tmp_path / "run"
    doc, _ = run_calibration(DATASET, out, source_dir=_coord_dir(tmp_path))
    rm = json.loads((out / "run_manifest.json").read_text())
    assert rm["code_sha"]  # the engine repo's HEAD
    assert "IMPORTED tos engine package" in rm["code_sha_source"]
    assert isinstance(rm["working_tree_clean"], bool)
    assert isinstance(rm["tree_status"], list)


def test_complete_requires_measured_equals_expected(tmp_path):
    coords = _coord_dir(tmp_path)
    # Dataset expects 2 targets; only 1 coordinate file is supplied.
    dataset = {"smoke": [
        {"pdb_id": "1CRN", "resolution": 1.5, "desc": "x"},
        {"pdb_id": "1UBQ", "resolution": 1.8, "desc": "missing"},
    ]}
    out = tmp_path / "run"
    doc, code = run_calibration(dataset, out, source_dir=coords)
    assert code == 1
    assert doc["summary"]["measured"] == 1
    assert doc["summary"]["expected_targets"] == 2
    assert doc["summary"]["run_status"] == "incomplete"
    # And the all-measured case is complete.
    out2 = tmp_path / "run2"
    doc2, code2 = run_calibration(DATASET, out2, source_dir=coords)
    assert code2 == 0
    assert doc2["summary"]["run_status"] == "complete"


def test_cli_refuses_dirty_engine_tree(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(
        runner_mod, "_engine_provenance",
        lambda: {"code_sha": "abc123", "code_sha_source": "test",
                 "working_tree_clean": False,
                 "tree_status": [" M backend/tos/engine/tier1_measurements.py"]})
    rc = main(["--dataset", "v23", "--out", str(tmp_path / "never")])
    assert rc == 2
    assert "dirty working tree" in capsys.readouterr().err
    assert not (tmp_path / "never").exists()  # nothing written


def test_cli_clean_tree_proceeds(tmp_path, monkeypatch):
    monkeypatch.setattr(
        runner_mod, "_engine_provenance",
        lambda: {"code_sha": "abc123", "code_sha_source": "test",
                 "working_tree_clean": True, "tree_status": []})
    rc = main(["--dataset", "v23", "--out", str(tmp_path / "ok"),
               "--source-dir", str(_coord_dir(tmp_path))])
    # v23 dataset has 50 targets; only 1 local file -> incomplete, exit 1
    assert rc == 1
    assert (tmp_path / "ok" / "run_manifest.json").exists()


def test_measure_target_uses_content_detected_format():
    cif = b"""data_AF-TEST
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.B_iso_or_equiv
ATOM N N GLY A 1 0.0 0.0 0.0 90.0
ATOM C CA GLY A 1 1.5 0.0 0.0 90.0
#
"""
    result = runner_mod.measure_target(cif, "predicted")
    assert result["atoms"] == 2
    assert result["is_experimental"] is False
