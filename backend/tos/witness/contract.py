import re
from typing import Set
AUTHORITY_DISCLAIMER = "Authority Notice: Non-authoritative interpretive witness. Deterministic engines hold exclusive authority."
FORBIDDEN_PHRASES: Set[str] = {"will work", "is likely to succeed", "promising drug", "therapeutic potential", "guaranteed"}
MAX_WITNESS_TOKENS = 1000
class WitnessViolation(Exception): pass
class WitnessContractEnforcer:
    def __init__(self):
        self.violation_count = 0
        self.max_violations = 2
        self.enabled = True
    def validate_output(self, text: str) -> str:
        if not self.enabled: raise WitnessViolation("Witness disabled")
        if AUTHORITY_DISCLAIMER not in text: self._record(); raise WitnessViolation("Missing disclaimer")
        for f in FORBIDDEN_PHRASES:
            if f in text.lower(): self._record(); raise WitnessViolation(f"Forbidden phrase: {f}")
        return text
    def _record(self):
        self.violation_count += 1
        if self.violation_count >= self.max_violations: self.enabled = False
