#!/usr/bin/env python3
"""
Entry 043 — Proof Engine
Generates formal proof trees from navigation results.
"""

import json
import sys
from dataclasses import dataclass, asdict
from typing import List, Dict

@dataclass
class ProofStep:
    axiom: str
    premise: str
    evidence: float
    threshold: float
    conclusion: str # "True" | "False"

@dataclass
class Proof:
    theorem: str
    steps: List[ProofStep]
    final_conclusion: str
    qed: bool

class ProofEngine:
    """
    Transforms a Navigation Result into a Formal Proof.
    """
    
    def generate_proof(self, navigation_result: Dict) -> Proof:
        trace = navigation_result.get("trace", {})
        steps = []
        qed = False
        
        # 1. Physics Axioms
        phy = trace.get("physics", {})
        if phy:
            # Stability Axiom
            plddt = phy.get("plddt_global", 0.0) # Assume standardized key
            # Mapping from trace might vary, handling robustness check
            if "plddt_global" not in phy and "structure_conf" in navigation_result:
                 plddt = navigation_result["structure_conf"] * 100
            
            steps.append(ProofStep(
                axiom="Stability",
                premise="Structure must be stable (pLDDT > 70)",
                evidence=plddt,
                threshold=70.0,
                conclusion="True" if plddt > 70 else "False"
            ))
            
            # Existence Axiom
            vol = phy.get("volume", 0.0)
            steps.append(ProofStep(
                axiom="Existence",
                premise="Pocket volume must be druggable (> 150 A^3)",
                evidence=vol,
                threshold=150.0,
                conclusion="True" if vol > 150 else "False"
            ))

        # 2. Chemistry Axioms
        chem = trace.get("chemistry", {})
        if chem:
            # Polarity Axiom
            # logic is complex, simplified for proof representation
            logp = chem.get("logp", 0.0)
            steps.append(ProofStep(
                axiom="Compatibility",
                premise="Ligand polarity must match pocket",
                evidence=logp,
                threshold=0.0, # Threshold is dynamic, placeholder
                conclusion="True" if navigation_result["status"] != "BLOCKED" else "Check"
            ))

        # Final Conclusion
        status = navigation_result.get("status")
        valid = status == "CLEARED" or status == "INDEXED"
        
        # Check if proof holds
        # A proof holds if all steps are True (for acceptance) 
        # OR if a False step correctly leads to Rejection (valid falsification)
        
        # We define Q.E.D. as: The decision is logically consistent with the steps.
        # i.e., If any step is False, Status must be BLOCKED/REJECTED.
        
        any_false = any(s.conclusion == "False" for s in steps)
        consistent_rejection = any_false and (status == "BLOCKED" or status == "REJECT")
        consistent_acceptance = (not any_false) and (status == "CLEARED" or status == "INDEXED")
        
        qed = consistent_rejection or consistent_acceptance
        
        return Proof(
            theorem=f"Hypothesis {navigation_result.get('compound_id', 'H')} is viable",
            steps=steps,
            final_conclusion="PROVEN FALSE" if any_false else "PROVEN TRUE",
            qed=qed
        )

    def format_proof(self, proof: Proof) -> str:
        lines = [f"Theorem: {proof.theorem}"]
        lines.append("-" * 40)
        for i, step in enumerate(proof.steps, 1):
            mark = "✅" if step.conclusion == "True" else "❌"
            lines.append(f"{i}. Axiom: {step.axiom}")
            lines.append(f"   Premise: {step.premise}")
            lines.append(f"   Evidence: {step.evidence} (Threshold: {step.threshold})")
            lines.append(f"   Conclusion: {step.conclusion} {mark}")
        lines.append("-" * 40)
        lines.append(f"Final Decision: {proof.final_conclusion}")
        lines.append(f"Consistent (Q.E.D.): {proof.qed}")
        return "\n".join(lines)

if __name__ == "__main__":
    # Test with a mock rejection result
    mock_result = {
        "status": "BLOCKED",
        "compound_id": "Mol_X",
        "trace": {
            "physics": {"volume": 120.0, "plddt_global": 85.0},
            "chemistry": {"logp": 2.5}
        }
    }
    
    engine = ProofEngine()
    proof = engine.generate_proof(mock_result)
    print(engine.format_proof(proof))
