import json, os, logging
from datetime import datetime, timezone
from pathlib import Path

import sys

try:
    import fcntl
except ModuleNotFoundError:
    fcntl = None

try:
    import msvcrt
except ModuleNotFoundError:
    msvcrt = None

from ..governance.constants import AuditConstants

logger = logging.getLogger("toscanini.nkg")

# WP-6: only these verdicts may enter the ledger.
VALID_VERDICTS = {"PASS", "VETO", "INDETERMINATE"}
# WP-6: SHA-256 of empty input — proves the "structure" never had content.
EMPTY_SOURCE_SHA256 = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"


class NKGManager:
    def __init__(self, path="/app/nkg/piu_moat.jsonl"):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def _already_recorded(self, artifact_hash: str, outcome_category: str) -> bool:
        """WP-6 duplicate policy (idempotent): a record with the same
        artifact hash AND outcome category is never written twice."""
        if not self.path.exists():
            return False
        with open(self.path, "r") as f:
            for line in f:
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if rec.get("artifact_hash") == artifact_hash and \
                   rec.get("outcome_category") == outcome_category:
                    return True
        return False

    def record_audit(self, payload):
        """Append one audit record to the NKG ledger.

        WP-6 validation — raises ValueError (record NOT written) when:
          - the verdict is missing or not one of PASS / VETO / INDETERMINATE
            (previously a missing verdict silently defaulted to VETO);
          - the structure hash is missing or equals the SHA-256 of empty
            input;
          - the provenance byte count is zero.

        Duplicate policy (idempotent): a record with the same artifact hash
        and outcome category is never written twice — returns False instead.
        """
        verdict = payload.get("verdict", {})
        governance = payload.get("governance", {})
        provenance = payload.get("provenance", {})
        tier3 = payload.get("tier3", {})
        architecture = payload.get("architecture", {})

        v = verdict.get("binary")
        if v not in VALID_VERDICTS:
            raise ValueError(
                f"NKG rejected: verdict {v!r} is not one of {sorted(VALID_VERDICTS)}"
            )

        source_hash = provenance.get("hash", "")
        if not source_hash or source_hash == EMPTY_SOURCE_SHA256:
            raise ValueError("NKG rejected: empty or blank-input structure hash")
        byte_count = provenance.get("byte_count")
        if byte_count == 0:
            raise ValueError("NKG rejected: zero source bytes")

        artifact_hash = source_hash[:AuditConstants.HASH_LEN]
        outcome_category = "SUCCESS" if v == "PASS" else "VETO"
        if self._already_recorded(artifact_hash, outcome_category):
            logger.info(f"NKG: duplicate suppressed for {artifact_hash} ({outcome_category})")
            return False

        phys = verdict.get("physical_score", 0)
        conf = verdict.get("confidence_score")
        conf_display = "UNAVAILABLE" if conf is None else f"{conf}%"
        
        record = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "artifact_hash": artifact_hash,
            "program_id": governance.get("program_id", "UNKNOWN"),
            "program_state": governance.get("state", "UNKNOWN"),
            "outcome_category": "SUCCESS" if v == "PASS" else "VETO",
            "failure_class": "NONE" if v == "PASS" else ("PHYSICAL_VETO" if phys < 100 else "CONFIDENCE_VETO"),
            "physical_score": phys,
            "confidence_score": conf,
            "epi_index": tier3.get("probability", 0),
            "architecture": architecture.get("authoritative_category", "UNKNOWN"),
            "s8_tax_applied": AuditConstants.S8_VETO_TAX if v == "VETO" else 0.0,
            "description": f"{v}: Phys={phys}%, Conf={conf_display}, EPI={tier3.get('probability', 0)}%"
        }
        
        line = json.dumps(record) + "\n"
        with open(self.path, "a") as f:
            _acquire_file_lock(f)
            try:
                f.write(line); f.flush(); os.fsync(f.fileno())
            finally:
                _release_file_lock(f)
            logger.info(f"NKG: {record['outcome_category']} — {record['artifact_hash']}")
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
                _acquire_file_lock(f)
                try:
                    f.write(line)
                    f.flush()
                    os.fsync(f.fileno())
                finally:
                    _release_file_lock(f)
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


def _acquire_file_lock(f):
    if fcntl is not None:
        fcntl.flock(f, fcntl.LOCK_EX)
    elif msvcrt is not None:
        msvcrt.locking(f.fileno(), msvcrt.LK_LOCK, 1)


def _release_file_lock(f):
    if fcntl is not None:
        fcntl.flock(f, fcntl.LOCK_UN)
    elif msvcrt is not None:
        try:
            msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, 1)
        except (OSError, IOError):
            pass


_nkg = None
def get_nkg():
    global _nkg
    if _nkg is None:
        # WP-6: path override for host/test environments. Production default
        # stays /app/nkg (nkg_data volume) unless TOSCANINI_NKG_PATH is set.
        default_path = os.environ.get("TOSCANINI_NKG_PATH", "/app/nkg/piu_moat.jsonl")
        _nkg = NKGManager(path=default_path)
    return _nkg
