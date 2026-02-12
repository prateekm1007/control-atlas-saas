#!/usr/bin/env python3
"""
Entry 038 — Biology Layer: Relevance Filtering
"""

import sys
from pathlib import Path
from dataclasses import dataclass

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "entries/036_unified_architecture"))
from core_interfaces import Layer, Constraint

from data_loader import get_target_data

@dataclass
class BiologyResult:
    status: str # PASS, FAIL, WARNING
    reasons: list
    metrics: dict

class BiologyLayer(Layer):
    """
    Evaluates:
    1. Essentiality (DepMap proxy)
    2. Disease Relevance (Driver status)
    """
    
    def evaluate(self, context: dict) -> BiologyResult:
        gene_name = context.get("gene_name")
        disease_context = context.get("disease", "cancer")
        
        data = get_target_data(gene_name)
        if not data:
            return BiologyResult("WARNING", ["Target data not found in KB"], {})
            
        reasons = []
        metrics = data
        
        # 1. Driver Status Check
        if not data.get("cancer_driver", False):
            reasons.append(f"{gene_name} is not a known driver in {disease_context}")
            
        # 2. Essentiality Check
        # Score > 0.5 implies dependency in cell lines
        ess_score = data.get("essentiality_score", 0.0)
        if ess_score < 0.3:
            reasons.append(f"Target non-essential (Score: {ess_score})")
            
        # Decision Logic
        status = "FAIL" if reasons else "PASS"
        return BiologyResult(status, reasons, metrics)

if __name__ == "__main__":
    # Test Case: GAPDH (Housekeeping, non-driver)
    engine = BiologyLayer()
    result = engine.evaluate({"gene_name": "GAPDH", "disease": "lung_cancer"})
    print(f"Target: GAPDH | Status: {result.status} | Reasons: {result.reasons}")
    
    # Test Case: KRAS (Driver)
    result = engine.evaluate({"gene_name": "KRAS", "disease": "lung_cancer"})
    print(f"Target: KRAS | Status: {result.status} | Reasons: {result.reasons}")
