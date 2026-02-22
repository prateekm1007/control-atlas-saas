import json, os, fcntl, logging
from datetime import datetime, timezone
from pathlib import Path
from ..governance.constants import AuditConstants

logger = logging.getLogger("toscanini.nkg")

class NKGManager:
    def __init__(self, path="/app/nkg/piu_moat.jsonl"):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def record_audit(self, payload):
        verdict = payload.get("verdict", {})
        governance = payload.get("governance", {})
        provenance = payload.get("provenance", {})
        tier3 = payload.get("tier3", {})
        architecture = payload.get("architecture", {})
        v = verdict.get("binary", "VETO")
        phys = verdict.get("physical_score", 0)
        conf = verdict.get("confidence_score", 0)
        
        record = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "artifact_hash": provenance.get("hash", "UNKNOWN")[:AuditConstants.HASH_LEN],
            "program_id": governance.get("program_id", "UNKNOWN"),
            "program_state": governance.get("state", "UNKNOWN"),
            "outcome_category": "SUCCESS" if v == "PASS" else "VETO",
            "failure_class": "NONE" if v == "PASS" else ("PHYSICAL_VETO" if phys < 100 else "CONFIDENCE_VETO"),
            "physical_score": phys,
            "confidence_score": conf,
            "epi_index": tier3.get("probability", 0),
            "architecture": architecture.get("authoritative_category", "UNKNOWN"),
            "s8_tax_applied": AuditConstants.S8_VETO_TAX if v == "VETO" else 0.0,
            "description": f"{v}: Phys={phys}%, Conf={conf}%, EPI={tier3.get('probability', 0)}%"
        }
        
        line = json.dumps(record) + "\n"
        with open(self.path, "a") as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            try:
                f.write(line); f.flush(); os.fsync(f.fileno())
                logger.info(f"NKG: {record['outcome_category']} — {record['artifact_hash']}")
            finally:
                fcntl.flock(f, fcntl.LOCK_UN)
        return True

    def read_records(self, limit=20):
        if not self.path.exists():
            return {"vetoes": [], "successes": []}
        records = []
        with open(self.path, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    try: records.append(json.loads(line))
                    except json.JSONDecodeError: continue
        return {
            "vetoes": [r for r in records if r.get("outcome_category") == "VETO"][-limit:],
            "successes": [r for r in records if r.get("outcome_category") == "SUCCESS"][-limit:]
        }

    def log_usage(self, record: dict):
        """
        Append a usage telemetry record.
        Separate from record_audit to avoid contaminating the existing NKG format.
        Writes to a dedicated usage log file.
        Failure is swallowed — telemetry must never break adjudication.
        """
        usage_path = self.path.parent / "usage_log.jsonl"
        try:
            line = json.dumps(record) + "\n"
            with open(usage_path, "a") as f:
                fcntl.flock(f, fcntl.LOCK_EX)
                try:
                    f.write(line)
                    f.flush()
                    os.fsync(f.fileno())
                finally:
                    fcntl.flock(f, fcntl.LOCK_UN)
        except Exception as e:
            logger.warning(f"Usage log write failed (non-fatal): {e}")

    def read_usage(self, limit=50):
        """Read recent usage records."""
        usage_path = self.path.parent / "usage_log.jsonl"
        if not usage_path.exists():
            return []
        records = []
        with open(usage_path, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        records.append(json.loads(line))
                    except json.JSONDecodeError:
                        continue
        return records[-limit:]


_nkg = None
def get_nkg():
    global _nkg
    if _nkg is None: _nkg = NKGManager()
    return _nkg
