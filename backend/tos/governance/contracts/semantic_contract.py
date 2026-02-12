"""
Toscanini Semantic Contract — v21.2.1
Prevents semantic drift between versions.
"""
from typing import Dict, Any

# IMMUTABLE: Only change via formal version bump
SEMANTIC_VERSION = "21.0"
SEMANTIC_CONTRACT = "v21_canonical"


class SemanticViolation(Exception):
    """Raised when payload violates semantic contract."""
    pass


def stamp_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Stamp payload with semantic version."""
    payload["semantic_version"] = SEMANTIC_VERSION
    payload["semantic_contract"] = SEMANTIC_CONTRACT
    return payload


def assert_semantic_compat(payload: Dict[str, Any]) -> None:
    """Validate payload semantic version."""
    if payload.get("semantic_version") != SEMANTIC_VERSION:
        raise SemanticViolation(
            f"Semantic version mismatch: {payload.get('semantic_version')} != {SEMANTIC_VERSION}"
        )


# Frozen legal scope text
LEGAL_SCOPE_TEXT = """THIS DOCUMENT CONFIRMS PHYSICAL PLAUSIBILITY ONLY.
IT DOES NOT PREDICT BIOLOGICAL ACTIVITY, SAFETY, OR EFFICACY.
ALL PROBABILITY VALUES ARE STRATEGIC PRIORITIZATION METRICS."""

PROBABILITY_MEANING = (
    "Structural feasibility under known physical constraints; "
    "not biological success probability."
)


def get_legal_scope() -> str:
    """Return immutable legal scope text."""
    return LEGAL_SCOPE_TEXT


def get_probability_meaning() -> str:
    """Return immutable probability interpretation."""
    return PROBABILITY_MEANING
