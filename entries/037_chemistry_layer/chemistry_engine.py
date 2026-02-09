#!/usr/bin/env python3
"""
Entry 037 — Chemistry Layer: Tractability & Compatibility
"""

import sys
from pathlib import Path
from dataclasses import dataclass

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "entries/036_unified_architecture"))
from core_interfaces import Layer, Constraint

from rdkit import Chem
from rdkit.Chem import QED, Descriptors
from conformer_gen import generate_conformers

@dataclass
class ChemistryResult:
    status: str # PASS, FAIL, WARNING
    reasons: list
    metrics: dict

class ChemistryLayer(Layer):
    """
    Evaluates:
    1. Intrinsic Tractability (QED as Warning)
    2. Pocket Compatibility (Polarity Match)
    3. Conformer Feasibility (Existence check)
    """
    
    def evaluate(self, context: dict) -> ChemistryResult:
        smiles = context.get("smiles")
        pocket_metrics = context.get("pocket_metrics", {})
        
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return ChemistryResult("FAIL", ["Invalid SMILES"], {})
            
        reasons = []
        warnings = []
        metrics = {}
        
        # 1. Intrinsic Tractability (QED) -> WARNING ONLY
        # QED encodes historical bias. We warn but do not reject.
        qed = QED.qed(mol)
        metrics["qed"] = round(qed, 2)
        if qed < 0.2:
            warnings.append(f"Low QED ({qed:.2f}) - potentially intractable")
            
        # 2. Molecular Weight vs Pocket Volume -> HEURISTIC CONSTRAINT
        mw = Descriptors.MolWt(mol)
        metrics["mw"] = round(mw, 1)
        pocket_vol = pocket_metrics.get("volume", 500.0)
        
        # Crude density approximation. 
        # Future: Replace with fragment occupancy or water displacement.
        if mw > pocket_vol * 1.5:
            reasons.append(f"Molecule likely too large for pocket ({mw} vs {pocket_vol})")
            
        # 3. Polarity Mismatch (LogP vs Hydrophobicity) -> HARD CONSTRAINT
        # High LogP in low-hydrophobicity pocket = Energetic penalty
        logp = Descriptors.MolLogP(mol)
        metrics["logp"] = round(logp, 2)
        pocket_hydro = pocket_metrics.get("hydrophobic_pct", 0.5)
        
        if pocket_hydro < 0.3 and logp > 4.5:
            reasons.append(f"Polarity Mismatch: Greasy ligand ({logp}) in polar pocket ({pocket_hydro})")
            
        # 4. Dynamics Check (Conformers) -> EXISTENTIAL CONSTRAINT
        # Checks if molecule can exist in 3D (sanity), not if it fits pocket (yet).
        confs = generate_conformers(mol, num_confs=5)
        if not confs:
            reasons.append("Sterically impossible (no valid conformers generated)")
            
        # Decision Logic
        if reasons:
            status = "FAIL"
        elif warnings:
            status = "WARNING"
        else:
            status = "PASS"
            
        return ChemistryResult(status, reasons + warnings, metrics)

if __name__ == "__main__":
    # Test Case: Sotorasib core vs Generic Pocket
    test_smiles = "CC(C)(C)NC1=NC=NC2=C1N=CN2" 
    test_pocket = {"volume": 400.0, "hydrophobic_pct": 0.6}
    
    engine = ChemistryLayer()
    result = engine.evaluate({"smiles": test_smiles, "pocket_metrics": test_pocket})
    print(f"Status: {result.status}")
    print(f"Metrics: {result.metrics}")
    print(f"Reasons: {result.reasons}")
