import json
from pathlib import Path

class GovernanceEnforcer:
    def __init__(self):
        root = Path(__file__).resolve().parent
        with open(root / "interaction_graph.json") as f: self.graph = json.load(f)
        with open(root / "write_permissions.json") as f: self.perms = json.load(f)
        with open(root / "architecture_rules.json") as f: self.arch_rules = json.load(f)

    def evaluate_interactions(self, all_results):
        triggered = []
        status_map = {l['law_id']: l['status'] for l in all_results}
        for rule in self.graph["interactions"]:
            if all(status_map.get(c.split(":")[0]) == c.split(":")[1] for c in rule["if"]):
                triggered.append(rule)
        return triggered

    def evaluate_architecture_eligibility(self, intent, all_results):
        penalties = []
        status_map = {l['law_id']: l['status'] for l in all_results}
        rule = self.arch_rules["rules"].get(intent)
        if rule:
            for req in rule["requires"]:
                law_id, required_status = req.split(":")
                if status_map.get(law_id) != required_status:
                    penalties.append({
                        "type": "ARCHITECTURE_INCOHERENCE",
                        "reason": rule["reason"],
                        "penalty": rule["penalty_if_missing"]
                    })
        return penalties
