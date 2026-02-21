"""
PIL-CAL-03 Modality Matrix Loader
----------------------------------
Single source of truth for method-aware law reclassification.
Replaces the inline conditional block in main.py.

The JSON file is the contract. This module is the reader.
Never hardcode law IDs or method strings here — they live in the JSON.
"""

from __future__ import annotations

import hashlib
import json
import logging
from functools import lru_cache
from pathlib import Path
from typing import Literal

logger = logging.getLogger(__name__)

# Canonical location — resolved relative to this file so imports
# work regardless of the caller's working directory.
_MATRIX_PATH = Path(__file__).parent / "modality_matrix.json"

MethodClassification = Literal["deterministic", "advisory_experimental", "heuristic"]

# The set of method strings that the JSON uses to mean "all experimental"
_ALL_EXPERIMENTAL_SENTINEL = "all_experimental"


@lru_cache(maxsize=1)
def _load_matrix() -> dict:
    """Load and cache the modality matrix JSON. Fails loudly on corruption."""
    if not _MATRIX_PATH.exists():
        raise FileNotFoundError(
            f"Modality matrix not found at {_MATRIX_PATH}. "
            "Run Phase 3 setup or restore from git."
        )
    with _MATRIX_PATH.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    # Basic schema guard
    required_keys = {"_meta", "rules"}
    missing = required_keys - data.keys()
    if missing:
        raise ValueError(
            f"modality_matrix.json is malformed: missing top-level keys {missing}"
        )
    return data


def compute_matrix_hash() -> str:
    """
    Return the SHA-256 hex digest of the raw modality_matrix.json bytes.
    Stamp this into every API response alongside the canon hash.
    """
    raw = _MATRIX_PATH.read_bytes()
    return hashlib.sha256(raw).hexdigest()


def get_matrix_meta() -> dict:
    """Return the _meta block for stamping into API payloads."""
    return _load_matrix()["_meta"]


def resolve_method(
    law_id: str,
    is_experimental: bool,
    method_type: str,
    default_method: MethodClassification = "deterministic",
) -> MethodClassification:
    """
    Determine the enforcement method for a given law and acquisition context.

    Args:
        law_id:          e.g. "LAW-120"
        is_experimental: True if the structure came from an experimental source
                         (x_ray, cryo_em, nmr) rather than a predicted model.
        method_type:     The acquisition method string from
                         structure.confidence.method — e.g. "nmr", "cryo_em",
                         "x_ray", "predicted".
        default_method:  The classification to return when no rule matches.
                         Callers should pass the law's native type from LAW_CANON
                         (deterministic or heuristic).

    Returns:
        The resolved MethodClassification. The caller is responsible for
        downgrading VETO → FAIL (Advisory) when this returns
        "advisory_experimental".
    """
    if not is_experimental:
        return default_method

    rules = _load_matrix()["rules"]

    for rule in rules:
        if rule["law_id"] != law_id:
            continue

        applicable_methods: list[str] = rule["applies_to_methods"]

        matches = (
            _ALL_EXPERIMENTAL_SENTINEL in applicable_methods
            or method_type in applicable_methods
        )

        if matches:
            reclassified = rule["reclassify_as"]
            logger.debug(
                "PIL-CAL-03: %s reclassified to '%s' for method '%s' (rationale: %s)",
                law_id,
                reclassified,
                method_type,
                rule.get("rationale", "N/A"),
            )
            return reclassified  # type: ignore[return-value]

    return default_method
