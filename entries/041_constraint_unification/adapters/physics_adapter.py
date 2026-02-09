"""
Maps Physics layer outputs → Constraints
"""
from ..constraint import Constraint

def volume_constraint(volume, min_volume, max_volume, plddt):
    # Check Min Volume
    if volume < min_volume:
        margin = (volume - min_volume) / min_volume
        status = "FAIL"
    # Check Max Volume
    elif volume > max_volume:
        margin = (max_volume - volume) / max_volume # Negative margin
        status = "FAIL"
    else:
        # Margin relative to optimal center (midpoint)
        mid = (min_volume + max_volume) / 2
        dist = abs(volume - mid)
        max_dist = (max_volume - min_volume) / 2
        margin = 1.0 - (dist / max_dist)
        status = "PASS"

    return Constraint(
        layer="Physics",
        parameter="PocketVolume",
        threshold=min_volume, # Simplified for display
        actual=volume,
        margin=round(margin, 3),
        status=status,
        collapses_space=True,
        provenance={
            "plddt": plddt,
            "min_vol": min_volume,
            "max_vol": max_volume,
            "source": "AlphaFoldDB",
            "logic": "Geometric viability"
        }
    )
