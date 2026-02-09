"""Bayesian probability engine for cross-tier causal probability."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict


@dataclass(frozen=True)
class ProbabilityInputs:
    plddt: float
    physical_boost: float = 0.25
    base_probability: float = 0.15
    refinement_tax: float = 0.25
    derivation_weight: float = 1.0


def compute_probability(inputs: ProbabilityInputs) -> float:
    """Compute the strategic probability lock described in v21.0.4."""
    p_ml = (inputs.plddt / 100.0) * 0.40
    probability = (inputs.base_probability + inputs.physical_boost + p_ml - inputs.refinement_tax)
    return max(0.0, min(1.0, probability * inputs.derivation_weight))


def compute_probability_breakdown(inputs: ProbabilityInputs) -> Dict[str, float]:
    """Return full probability breakdown for auditability."""
    p_ml = (inputs.plddt / 100.0) * 0.40
    pre_weight = inputs.base_probability + inputs.physical_boost + p_ml - inputs.refinement_tax
    probability = max(0.0, min(1.0, pre_weight * inputs.derivation_weight))
    return {
        "p_base": inputs.base_probability,
        "p_phys": inputs.physical_boost,
        "p_ml": p_ml,
        "tax": inputs.refinement_tax,
        "weight_derivation": inputs.derivation_weight,
        "pre_weight": pre_weight,
        "probability": probability,
    }
