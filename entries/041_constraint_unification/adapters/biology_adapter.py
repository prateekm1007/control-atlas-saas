"""
Maps Biology context → Constraints
"""
from ..constraint import Constraint

def essentiality_constraint(score):
    threshold = 0.5
    # Margin: How far above 0.5?
    # 1.0 -> margin 1.0
    # 0.5 -> margin 0.0
    # 0.0 -> margin -1.0
    
    if score >= threshold:
        margin = (score - threshold) / (1.0 - threshold)
        status = "PASS"
    else:
        margin = (score - threshold) / threshold
        status = "FAIL"

    return Constraint(
        layer="Biology",
        parameter="Essentiality",
        threshold=threshold,
        actual=score,
        margin=round(margin, 3),
        status=status,
        collapses_space=True,
        provenance={
            "source": "DepMap",
            "logic": "Causal relevance"
        }
    )
