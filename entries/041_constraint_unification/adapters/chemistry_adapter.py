"""
Maps Chemistry metrics → Constraints
"""
from ..constraint import Constraint

def polarity_constraint(logp, pocket_hydro):
    # Heuristic target logP range based on pocket hydrophobicity
    # Polar pocket (hydro < 0.4) -> Needs polar ligand (LogP < 2)
    # Greasy pocket (hydro > 0.6) -> Needs greasy ligand (LogP > 3)
    
    if pocket_hydro > 0.6:
        target_logp = 3.5
        tolerance = 2.5
    elif pocket_hydro < 0.4:
        target_logp = 1.0
        tolerance = 2.0
    else:
        target_logp = 2.5
        tolerance = 3.0
        
    dist = abs(logp - target_logp)
    margin = 1.0 - (dist / tolerance)
    
    # Hard failures for extreme mismatches
    if margin < -0.5:
        status = "FAIL"
        collapses = True
    elif margin < 0:
        status = "WARNING"
        collapses = False
    else:
        status = "PASS"
        collapses = False

    return Constraint(
        layer="Chemistry",
        parameter="PolarityMatch",
        threshold=target_logp,
        actual=logp,
        margin=round(margin, 3),
        status=status,
        collapses_space=collapses,
        provenance={
            "pocket_hydrophobicity": pocket_hydro,
            "target_logp": target_logp,
            "logic": "Energetic compatibility"
        }
    )
