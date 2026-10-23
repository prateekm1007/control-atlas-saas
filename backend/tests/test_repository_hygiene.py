"""Repository hygiene: one engine source and no known-dead routing code."""
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_nested_duplicate_engine_tree_is_absent():
    assert not (REPO_ROOT / "control-atlas-saas").exists(), (
        "stale nested control-atlas-saas/ duplicate engine tree is present"
    )


def test_dead_router_is_absent():
    assert not (REPO_ROOT / "backend" / "tos" / "router").exists(), (
        "dead tos.router package is present"
    )
