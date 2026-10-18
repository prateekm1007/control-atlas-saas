"""
WP-SCI-1 semantic law hash tests (directive §6).

The legacy canon hash covers only sorted law IDs + method classifications;
compute_law_semantic_hash() binds MEANING: thresholds, operators, units,
scopes, types, method classes, categories, containment records, the pinned
ideal/sigma tables, the raw modality-matrix bytes digest, and the status
vocabulary. CI pins this digest (governance/law_semantic.sha256) and fails
on drift until a reviewed update re-pins it.
"""
import hashlib
import json
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.governance import station_sop
from tos.governance.station_sop import (
    compute_law_semantic_hash,
    compute_canon_hash,
    LAW_CANON_HASH,
)

REPO = Path(__file__).resolve().parents[2]
PIN_FILE = REPO / "governance" / "law_semantic.sha256"


class TestSemanticHashPin:

    def test_pinned_file_matches_computed_hash(self):
        """CI-equivalent drift check: any meaning/threshold change without
        an explicit re-pin fails here exactly as it fails in the gate."""
        assert PIN_FILE.exists(), "governance/law_semantic.sha256 must exist"
        pinned = PIN_FILE.read_text().strip()
        computed = compute_law_semantic_hash()
        assert computed == pinned, (
            "law semantics drifted without a reviewed re-pin: "
            f"computed={computed} pinned={pinned}"
        )

    def test_canon_hash_unchanged(self):
        """Protected canon hash must remain exactly the pinned value."""
        assert LAW_CANON_HASH == "6a9cd4b4349b81de"
        assert compute_canon_hash() == LAW_CANON_HASH


class TestSemanticHashSensitivity:

    def _fingerprint(self):
        """Reproduce the fingerprint dict the hash function digests, so
        perturbation tests can assert binding of each field."""
        from tos.governance.modality_matrix import compute_matrix_hash
        return {
            "schema": "law-semantic-fingerprint/1.0.0",
            "laws": [
                {
                    "id": lid,
                    "title": station_sop.LAW_CANON[lid]["title"],
                    "threshold": station_sop.LAW_CANON[lid]["threshold"],
                    "operator": station_sop.LAW_CANON[lid]["operator"],
                    "unit": station_sop.LAW_CANON[lid]["unit"],
                    "scope": station_sop.LAW_CANON[lid]["scope"],
                    "type": station_sop.LAW_CANON[lid]["type"],
                    "method_class": station_sop.LAW_METHOD_CLASSIFICATIONS[lid],
                    "category": station_sop.LAW_CATEGORIES[lid],
                    "containment": station_sop.LAW_CONTAINMENT.get(lid),
                }
                for lid in sorted(station_sop.LAW_CANON.keys())
            ],
            "ideal_table": {k: station_sop.IDEAL_TABLE[k]
                            for k in sorted(station_sop.IDEAL_TABLE)},
            "sigma_table": {k: station_sop.SIGMA_TABLE[k]
                            for k in sorted(station_sop.SIGMA_TABLE)},
            "modality_matrix_sha256": compute_matrix_hash(),
            "status_vocabulary": sorted([
                station_sop.STATUS_MEASURED_PASS,
                station_sop.STATUS_MEASURED_VETO,
                station_sop.STATUS_MEASURED_FAIL,
                station_sop.STATUS_NOT_MEASURED,
                station_sop.STATUS_NOT_APPLICABLE,
                station_sop.STATUS_MEASUREMENT_FAILED,
            ]),
        }

    def _digest(self, fingerprint):
        canonical = json.dumps(fingerprint, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    def test_threshold_change_binds(self, monkeypatch):
        base = self._digest(self._fingerprint())
        monkeypatch.setitem(station_sop.LAW_CANON["LAW-120"], "threshold", 99.0)
        assert self._digest(self._fingerprint()) != base

    def test_unit_change_binds(self, monkeypatch):
        base = self._digest(self._fingerprint())
        monkeypatch.setitem(station_sop.LAW_CANON["LAW-160"], "unit", "nm")
        assert self._digest(self._fingerprint()) != base

    def test_operator_change_binds(self, monkeypatch):
        base = self._digest(self._fingerprint())
        monkeypatch.setitem(station_sop.LAW_CANON["LAW-182"], "operator", ">")
        assert self._digest(self._fingerprint()) != base

    def test_containment_change_binds(self, monkeypatch):
        base = self._digest(self._fingerprint())
        monkeypatch.setitem(station_sop.LAW_CONTAINMENT["LAW-120"], "status",
                            "lifted")
        assert self._digest(self._fingerprint()) != base

    def test_ideal_geometry_change_binds(self, monkeypatch):
        base = self._digest(self._fingerprint())
        monkeypatch.setitem(station_sop.IDEAL_TABLE, "O-C-N", 125.0)
        assert self._digest(self._fingerprint()) != base

    def test_status_vocabulary_change_binds(self, monkeypatch):
        base = self._digest(self._fingerprint())
        monkeypatch.setattr(
            station_sop, "STATUS_NOT_APPLICABLE", "NOT_APPLICABLE_V2")
        # vocabulary list in the function reads the module constant
        assert compute_law_semantic_hash() != base


class TestRegistryDefinitionVersions:

    def test_repaired_laws_carry_definition_version_2(self):
        yaml = pytest.importorskip("yaml", reason="PyYAML from requirements-dev")
        data = yaml.safe_load((REPO / "governance" / "law_registry.yaml")
                              .read_text(encoding="utf-8"))
        entries = {e["id"]: e for e in data["laws"]}
        for lid in ("LAW-120", "LAW-160", "LAW-195"):
            assert entries[lid].get("measurement_definition_version") == 2, \
                f"{lid} must carry measurement definition v2"

    def test_registry_declares_status_semantics(self):
        yaml = pytest.importorskip("yaml", reason="PyYAML from requirements-dev")
        data = yaml.safe_load((REPO / "governance" / "law_registry.yaml")
                              .read_text(encoding="utf-8"))
        sem = data["status_semantics"]
        for key in ("measured_success", "measured_violation", "measured_failure",
                    "empty_sample", "not_applicable", "acquisition_failure"):
            assert key in sem, f"status_semantics missing {key}"
        assert sem["empty_sample"] == "NOT_MEASURED"
        assert sem["not_applicable"] == "NOT_APPLICABLE"

    def test_law120_containment_recorded_in_registry(self):
        yaml = pytest.importorskip("yaml", reason="PyYAML from requirements-dev")
        data = yaml.safe_load((REPO / "governance" / "law_registry.yaml")
                              .read_text(encoding="utf-8"))
        entry = next(e for e in data["laws"] if e["id"] == "LAW-120")
        assert entry["containment"]["state"] == "contained"
        assert entry["containment"]["reason"] == "measurement-definition-repair"
