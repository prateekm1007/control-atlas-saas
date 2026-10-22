"""
WP-8 — one law registry, using the pinned code as authority.

Asserts:
  - every LAW id the engine can emit appears in the registry;
  - registry values match the pinned LAW_CANON / LAW_METHOD_CLASSIFICATIONS
    EXACTLY (no drift in either direction);
  - every enforced entry has an implementation and a test file;
  - the dedicated-test gap is exactly the six known laws (new gaps fail CI);
  - retired laws are exactly {190, 210, 220, 230};
  - the registry content hash is pinned (recorded in
    governance/law_registry.sha256 and as a constant here);
  - the glossary describes the pinned measurements (no stale physics).
"""
import hashlib
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
REGISTRY = REPO / "governance" / "law_registry.yaml"
SHA256_FILE = REPO / "governance" / "law_registry.sha256"
MEASUREMENTS = REPO / "backend" / "tos" / "engine" / "tier1_measurements.py"

# Recorded content hash of governance/law_registry.yaml (WP-8 requirement).
REGISTRY_SHA256 = "28e0966705ed97162a2e52228b71b311e7d1bede5d6343825db6a38363420125"

yaml = pytest.importorskip("yaml", reason="PyYAML comes from requirements-dev.txt")

REGISTRY_DATA = yaml.safe_load(REGISTRY.read_text(encoding="utf-8"))
ENTRIES = {e["id"]: e for e in REGISTRY_DATA["laws"]}

from tos.governance.station_sop import LAW_CANON, LAW_METHOD_CLASSIFICATIONS, RETIRED_LAWS

# Laws without a DEDICATED test file — the known gap recorded by the brief.
KNOWN_DEDICATED_GAPS = {"LAW-100", "LAW-105", "LAW-110", "LAW-125", "LAW-130", "LAW-170"}


class TestRegistryCoversMeasurements:

    def test_every_emitted_law_id_is_in_registry(self):
        src = MEASUREMENTS.read_text()
        emitted = set(re.findall(r'results\["(LAW-\d+)"\]', src))
        assert emitted, "no LAW ids found in the measurement module — regex drift?"
        missing = emitted - set(ENTRIES)
        assert not missing, f"measurement results reference laws missing from the registry: {sorted(missing)}"

    def test_registry_enforced_set_matches_pinned_canon(self):
        enforced = {k for k, e in ENTRIES.items() if e["status"] == "enforced"}
        assert enforced == set(LAW_CANON.keys()) - RETIRED_LAWS


class TestRegistryMatchesPinnedValues:

    @pytest.mark.parametrize("law_id", sorted(LAW_CANON.keys()))
    def test_values_copied_exactly(self, law_id):
        canon = LAW_CANON[law_id]
        entry = ENTRIES[law_id]
        assert entry["title"] == canon["title"]
        assert float(entry["threshold"]) == float(canon["threshold"])
        assert entry["units"] == canon["unit"]
        assert entry["operator"] == canon["operator"]
        assert entry["scope"] == canon["scope"]
        assert entry["type"] == canon["type"]
        assert entry["tier"] == LAW_METHOD_CLASSIFICATIONS[law_id]

    def test_retired_laws_exactly_match_decision(self):
        """Brief §4.2 (owner decision): 190/210/220/230 retired, never renumbered."""
        retired = {k for k, e in ENTRIES.items() if e["status"] == "retired"}
        assert retired == {"LAW-155", "LAW-190", "LAW-210", "LAW-220", "LAW-230"}
        for k in retired - {"LAW-155"}:
            assert ENTRIES[k]["implementation"] is None
            assert ENTRIES[k]["test"] is None
        assert ENTRIES["LAW-155"]["implementation"] is not None


class TestEnforcedEntriesHaveImplementationAndTests:

    @pytest.mark.parametrize("law_id", sorted(LAW_CANON.keys()))
    def test_implementation_and_test_file_present(self, law_id):
        entry = ENTRIES[law_id]
        impl = entry["implementation"]
        assert impl, f"{law_id}: enforced law must name its implementation"
        impl_path = REPO / impl.split(" (")[0].split("::")[0]
        assert impl_path.exists(), f"{law_id}: implementation path does not exist: {impl}"

        test_ref = entry["test"]
        assert test_ref and test_ref["file"], f"{law_id}: enforced law must name a test file"
        test_path = REPO / test_ref["file"]
        assert test_path.exists(), f"{law_id}: test file does not exist: {test_ref['file']}"

    def test_dedicated_test_gap_is_exactly_the_known_six(self):
        """If a new law loses its dedicated test (or one of the six gains one),
        this assertion forces a registry + PR update."""
        indirect = {k for k, e in ENTRIES.items()
                    if e["status"] == "enforced" and not e["test"]["dedicated"]}
        assert indirect == KNOWN_DEDICATED_GAPS


class TestContentHashRecorded:

    def test_sha256_file_matches_content(self):
        recorded = SHA256_FILE.read_text().strip()
        actual = hashlib.sha256(REGISTRY.read_bytes()).hexdigest()
        assert recorded == actual

    def test_pinned_constant_matches_content(self):
        actual = hashlib.sha256(REGISTRY.read_bytes()).hexdigest()
        assert REGISTRY_SHA256 == actual


class TestGlossaryMatchesPinnedPhysics:

    def test_glossary_ids_cover_registry(self):
        from tos.glossary.law_glossary import list_all_law_ids
        assert set(list_all_law_ids()) == set(ENTRIES.keys())

    def test_law_155_describes_voxel_occupancy_not_no_overlap(self):
        from tos.glossary.law_glossary import get_law_explanation
        exp = get_law_explanation("LAW-155")
        assert exp["title"] == "Voxel Occupancy"
        assert exp["status"] == "retired"
        assert "voxel" in exp["explanation"].lower()
        assert "retired" in exp["explanation"].lower()
        assert "cannot overlap" not in exp["explanation"].lower()

    def test_law_170_describes_residue_identity_not_pi(self):
        from tos.glossary.law_glossary import get_law_explanation
        exp = get_law_explanation("LAW-170")
        assert exp["title"] == "Residue Identity"

    def test_retired_law_explains_retirement(self):
        from tos.glossary.law_glossary import get_law_explanation
        exp = get_law_explanation("LAW-190")
        assert exp["status"] == "retired"

    def test_stale_glossary_files_rewritten(self):
        """The five conflicting glossary files must carry the pinned titles."""
        expected = {
            "backend/tos/glossary/tier1/law_155_steric.py": "Voxel Occupancy",
            "backend/tos/glossary/tier1/law_120_bond.py": "Bond Angle",
            "backend/tos/glossary/tier1/law_160_backbone.py": "Chain Integrity",
            "backend/tos/glossary/tier1/law_200_cavity.py": "Packing Quality",
            "backend/tos/glossary/tier2/law_170_pi.py": "Residue Identity",
        }
        for rel, title in expected.items():
            src = (REPO / rel).read_text()
            assert f'"{title}"' in src, f"{rel} still carries the stale title (expected {title})"


class TestHandoffMatchesRegistry:

    def test_handoff_carries_pinned_labels(self):
        handoff = (REPO / "SAAS_HANDOFF.md").read_text()
        assert "LAW-155 is **Voxel Occupancy**" in handoff
        assert "LAW-170 is **Residue Identity**" in handoff
        assert "LAW-190/210/220/230 are **retired**" in handoff

    def test_handoff_has_no_conflict_markers(self):
        handoff = (REPO / "SAAS_HANDOFF.md").read_text()
        assert not re.search(r"^<<<<<<<|^>>>>>>>", handoff, flags=re.M)
