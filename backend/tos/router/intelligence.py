import os, yaml
from tos.ingestion.structure_object import StructureObject

class IntelligenceRouter:
    def decide(self, structure: StructureObject, tier1_verdict: str) -> dict:
        if tier1_verdict == "VETO":
            return {"state": "NON-REAL", "banner": "❌ Physical Impossibility Detected", "recommendation": "Structure violates physics."}
        if structure.confidence.is_low:
            return {"state": "UNCERTAIN", "banner": f"🟡 Physically Valid · Low Confidence ({structure.confidence.mean_plddt:.1f}%)", "recommendation": "Escalation recommended."}
        return {"state": "VALID", "banner": f"🟢 Physically Real · High Confidence ({structure.confidence.mean_plddt:.1f}%)", "recommendation": "Safe for docking."}
