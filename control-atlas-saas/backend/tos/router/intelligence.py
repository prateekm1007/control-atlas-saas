import yaml
import os
from ingestion.structure_object import StructureObject

class IntelligenceRouter:
    def __init__(self):
        config_path = os.path.join(os.path.dirname(__file__), "config.yaml")
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)

    def decide(self, structure: StructureObject, tier1_verdict: str) -> dict:
        t = self.config['thresholds']
        
        if tier1_verdict == "VETO":
            return {
                "state": "NON-REAL",
                "color": "red",
                "banner": "❌ Physical Impossibility Detected",
                "recommendation": "Structure violates invariant physical laws. Reject immediately."
            }
        
        # Routing logic: Only uses Protein pLDDT (Ligand confidence ignored per policy)
        if structure.confidence.mean_plddt < t['low_plddt']:
            return {
                "state": "REAL_BUT_UNCERTAIN",
                "color": "orange",
                "banner": f"🟡 Physically Valid · Low Confidence ({structure.confidence.mean_plddt:.1f}%)",
                "recommendation": f"Valid geometry but high ML uncertainty. Escalate to AF3 or Chai-1."
            }
            
        return {
            "state": "REAL_AND_VALID",
            "color": "green",
            "banner": f"🟢 Physically Real · High Confidence ({structure.confidence.mean_plddt:.1f}%)",
            "recommendation": "Structure is physically sound and high-confidence. Safe for engineering."
        }
