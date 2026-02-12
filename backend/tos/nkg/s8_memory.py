import json, os
from datetime import datetime
from pathlib import Path
from .similarity import SimilarityEngine
from ..governance.constants import FAILURE_SEVERITY, S8_SIMILARITY_THRESHOLD

class S8NegativeMemory:
    """PILLAR 16: Institutional Memory."""
    def __init__(self, storage_path="/app/nkg/piu_moat.jsonl"):
        self.path = Path(storage_path)
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def ingest_wet_lab_failure(self, artifact_hash, failure_type, description, fingerprint):
        """
        S8 Ingestion: Only accepts if hash is validated as SEALED.
        (Validation check performed in main.py gateway)
        """
        record = {
            "timestamp": datetime.now().isoformat(),
            "artifact_hash": artifact_hash,
            "failure_type": failure_type,
            "description": description,
            "fingerprint": fingerprint,
            "governance_status": "SEALED_LINK"
        }
        with open(self.path, "a") as f:
            f.write(json.dumps(record) + "\n")
        return {"status": "S8_RECORDED", "artifact_seal": artifact_hash}

    def get_s8_signal(self, current_fingerprint: dict) -> dict:
        if not self.path.exists(): return {"match_count": 0}
        matches = []
        with open(self.path, "r") as f:
            for line in f:
                try:
                    failure = json.loads(line)
                    s = SimilarityEngine.compute_similarity(current_fingerprint, failure['fingerprint'])
                    if s >= S8_SIMILARITY_THRESHOLD:
                        matches.append({"s": s, "type": failure['failure_type']})
                except: continue
        
        if not matches: return {"match_count": 0}
        best_match = max(matches, key=lambda x: x['s'])
        return {
            "match_count": len(matches),
            "failure_rate": min(1.0, len(matches) / 10.0),
            "similarity_score": best_match['s'],
            "failure_class": best_match['type'],
            "severity_coeff": FAILURE_SEVERITY.get(best_match['type'], 0.50)
        }

_s8 = S8NegativeMemory()
def get_s8_memory(): return _s8
