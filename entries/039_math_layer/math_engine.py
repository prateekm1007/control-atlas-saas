#!/usr/bin/env python3
"""
Entry 039 — Math Layer: Robustness Calculation
Integrates margins from Physics, Chemistry, and Biology layers.
"""

import sys
from pathlib import Path
from dataclasses import dataclass

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "entries/036_unified_architecture"))
from core_interfaces import Layer

@dataclass
class MathResult:
    robustness: float
    advice: str

class MathLayer:
    """
    Calculates robustness of the survival path.
    Robustness = Product of margins from all layers.
    """
    
    def calculate_robustness(self, physics_metrics, chemistry_metrics, biology_metrics) -> MathResult:
        score = 1.0
        advice = []
        
        # 1. Physics Robustness (Volume match)
        # Ideal: Ligand MW ~ Pocket Volume (density approx)
        vol = physics_metrics.get("volume", 500.0)
        mw = chemistry_metrics.get("mw", 300.0)
        
        if vol > 0:
            fill_ratio = mw / vol 
            # Ideal fill ratio is ~0.8-1.2. Penalize deviations.
            if 0.5 <= fill_ratio <= 1.5:
                score *= 0.9  # Good fit
                advice.append("Good volume match")
            else:
                score *= 0.5  # Poor fit (too small or too big)
                advice.append(f"Volume mismatch (Fill: {fill_ratio:.2f})")
        
        # 2. Chemistry Robustness (Polarity match)
        # Ideal: Pocket Hydrophobicity matches Ligand LogP
        pocket_hydro = physics_metrics.get("hydrophobic_pct", 0.5)
        logp = chemistry_metrics.get("logp", 2.0)
        
        # Heuristic alignment
        if pocket_hydro > 0.6 and logp > 3.0:
            score *= 0.9 # Matched greasy
            advice.append("Hydrophobic match")
        elif pocket_hydro < 0.4 and logp < 2.0:
            score *= 0.9 # Matched polar
            advice.append("Polar match")
        else:
            score *= 0.7 # Neutral/Mixed
            
        # 3. Biology Robustness (Essentiality)
        # High essentiality = High robustness
        ess_score = biology_metrics.get("essentiality_score", 0.0)
        if ess_score > 0.8:
            score *= 1.0 # Essential
            advice.append("High biological relevance")
        elif ess_score > 0.5:
            score *= 0.8
        else:
            score *= 0.5 # Weak relevance
            advice.append("Low biological relevance")
            
        return MathResult(round(score, 3), "; ".join(advice))

if __name__ == "__main__":
    # Test
    eng = MathLayer()
    res = eng.calculate_robustness(
        {"volume": 400.0, "hydrophobic_pct": 0.7},
        {"mw": 350.0, "logp": 3.5},
        {"essentiality_score": 0.9}
    )
    print(f"Robustness: {res.robustness}")
    print(f"Advice: {res.advice}")
