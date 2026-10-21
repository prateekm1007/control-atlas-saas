"""
P0 item 7 (2026-10-19 audit): every LAW ID has exactly one definition.

The dead duplicate `Tier1Laws._law_155` in engine/tier1_laws.py (no
callers anywhere) implemented "Steric Clash Prohibition" while the live
engine's LAW-155 is Voxel Occupancy — a silent competing source of
truth. The module is deleted; these tests make the single-definition
invariant structural so a reintroduction fails CI:

  1. every LAW ID in LAW_CANON is emitted by exactly one assignment in
     the measurement engine (no duplicate implementations);
  2. the engine emits exactly the canon's law set (no orphans, no gaps);
  3. no `def _law_*` helper survives anywhere under backend/tos (the
     dead-duplicate pattern is banned);
  4. the governance registry's law set equals the canon's (the registry
     is a projection of the canon, not a second authority).
"""
import re
from pathlib import Path

import pytest

from tos.governance.station_sop import LAW_CANON

BACKEND = Path(__file__).resolve().parents[1]
REPO_ROOT = BACKEND.parent
TOS = BACKEND / "tos"
MEASUREMENTS = TOS / "engine" / "tier1_measurements.py"
REGISTRY = REPO_ROOT / "governance" / "law_registry.yaml"
NESTED_COPY = REPO_ROOT / "control-atlas-saas"


def _emitted_law_ids(source: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for match in re.finditer(r'results\["(LAW-\d+)"\]', source):
        counts[match.group(1)] = counts.get(match.group(1), 0) + 1
    return counts


def _measurement_block_banners(source: str) -> dict[str, int]:
    """One `# LAW-<id>:` banner per law = one measurement block per law."""
    counts: dict[str, int] = {}
    for match in re.finditer(r'# (LAW-\d+):', source):
        counts[match.group(1)] = counts.get(match.group(1), 0) + 1
    return counts


def test_every_law_id_has_exactly_one_measurement_block():
    source = MEASUREMENTS.read_text(encoding="utf-8")
    banners = _measurement_block_banners(source)
    duplicates = {lid: n for lid, n in banners.items() if n != 1}
    assert not duplicates, f"LAW IDs with != 1 measurement block: {duplicates}"
    # every canon law also has at least one assignment site (measured or
    # honest-status branch)
    assigned = set(_emitted_law_ids(source))
    assert set(LAW_CANON) <= assigned


def test_engine_emits_exactly_the_canon_set():
    source = MEASUREMENTS.read_text(encoding="utf-8")
    emitted = set(_emitted_law_ids(source))
    canon = set(LAW_CANON)
    assert emitted == canon, (
        f"missing from engine: {sorted(canon - emitted)}; "
        f"orphan emissions: {sorted(emitted - canon)}")


def test_no_dead_law_helper_definitions_anywhere():
    offenders = []
    for py in TOS.rglob("*.py"):
        for i, line in enumerate(py.read_text(encoding="utf-8",
                                              errors="ignore").splitlines(), 1):
            if re.match(r"\s*def _law_\d+", line):
                offenders.append(f"{py.relative_to(BACKEND)}:{i}: {line.strip()}")
    assert not offenders, (
        "dead `_law_*` implementation found (use the canonical engine):\n"
        + "\n".join(offenders))


def test_deleted_dead_module_stays_deleted():
    assert not (TOS / "engine" / "tier1_laws.py").exists(), (
        "engine/tier1_laws.py was the dead duplicate LAW-155 "
        "(Steric Clash Prohibition) competing with the live Voxel "
        "Occupancy LAW-155; it must not come back")


def test_nested_duplicate_engine_tree_is_absent():
    """The tracked nested control-atlas-saas copy was a second source of
    truth. It must be deleted wholesale, not silently ignored by a scan
    rooted at backend/tos."""
    assert not NESTED_COPY.exists(), (
        "tracked nested control-atlas-saas/ copy is still present")
    nested_engines = [
        p for p in REPO_ROOT.rglob("tier1_measurements.py")
        if ".git" not in p.parts and p != MEASUREMENTS
    ]
    nested_laws = [
        p for p in REPO_ROOT.rglob("tier1_laws.py")
        if ".git" not in p.parts and p != (TOS / "engine" / "tier1_laws.py")
    ]
    assert not nested_engines, f"duplicate engine implementation: {nested_engines}"
    assert not nested_laws, f"dead law implementation: {nested_laws}"


def test_registry_enforced_set_equals_canon():
    yaml = pytest.importorskip("yaml", reason="PyYAML comes from requirements-dev.txt")
    data = yaml.safe_load(REGISTRY.read_text(encoding="utf-8"))
    enforced = {e["id"] for e in data["laws"] if e.get("status") == "enforced"}
    retired = {e["id"] for e in data["laws"] if e.get("status") == "retired"}
    assert enforced == set(LAW_CANON), (
        f"registry-only enforced: {sorted(enforced - set(LAW_CANON))}; "
        f"canon-only: {sorted(set(LAW_CANON) - enforced)}")
    # Retired laws are visibly distinguished (Constitution EC-4.7) and are
    # exactly the four recorded in the WP-8 registry brief.
    assert retired == {"LAW-190", "LAW-210", "LAW-220", "LAW-230"}
