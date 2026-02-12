"""Geometric architecture derivation for Tier-2 architecture lab."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass(frozen=True)
class ArchitectureDecision:
    category: str
    voxel_span: Tuple[float, float, float]
    rationale: str
    symmetry_score: Optional[float] = None
    clash_free_volume: Optional[float] = None


def derive_architecture(
    voxel_span: Tuple[float, float, float],
    symmetry_score: Optional[float] = None,
    clash_free_volume: Optional[float] = None,
) -> ArchitectureDecision:
    """Derive architecture category from voxel span and spatial metadata.

    Categories: Linear, Multivalent, Metal, Helical, Sheet, Compact, Extended.
    """
    x, y, z = voxel_span
    max_dim = max(voxel_span)
    min_dim = min(voxel_span)
    if symmetry_score is not None and symmetry_score >= 0.85 and max_dim <= 50:
        category = "Compact"
        rationale = "High symmetry with bounded span indicates a compact topology."
    elif clash_free_volume is not None and clash_free_volume >= 150_000 and x >= 55 and y >= 55:
        category = "Multivalent"
        rationale = "Large clash-free volume and broad footprint imply multivalent architecture."
    elif max_dim <= 40:
        category = "Compact"
        rationale = "Voxel span is tightly bounded across all axes."
    elif max_dim / max(min_dim, 1.0) >= 3:
        category = "Linear"
        rationale = "One axis dominates the voxel span, indicating a linear architecture."
    elif z >= 60:
        category = "Extended"
        rationale = "Strong extension along the z-axis suggests an extended scaffold."
    elif x >= 55 and y >= 55:
        category = "Multivalent"
        rationale = "Broad footprint across two axes implies multivalent topology."
    elif min_dim <= 15:
        category = "Metal"
        rationale = "Tight core thickness indicates metal-binding compactness."
    elif y >= 45:
        category = "Sheet"
        rationale = "Broad planar span suggests sheet-dominant topology."
    else:
        category = "Helical"
        rationale = "Balanced dimensions imply helical bundle geometry."

    return ArchitectureDecision(category=category, voxel_span=voxel_span, rationale=rationale)
