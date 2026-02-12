import uuid, time, hashlib
from enum import Enum
from dataclasses import dataclass, field
from typing import Tuple, Optional, Dict
from ..utils.type_guards import force_bytes

class ProgramState(Enum):
    INSTANTIATED = 1
    ACQUIRING = 2
    AUDITED = 3
    DECIDED = 4
    SEALED = 5

@dataclass(frozen=True)
class ArchitectureDecision:
    category: str
    intent: str
    rationale: str
    prev_hash: str
    timestamp: float = field(default_factory=time.time)
    def compute_hash(self) -> str:
        return hashlib.sha256(force_bytes(f"{self.prev_hash}|{self.category}|{self.timestamp}")).hexdigest()[:16]

@dataclass(frozen=True)
class DiscoveryProgram:
    program_id: str = field(default_factory=lambda: "PRG-" + str(uuid.uuid4())[:8])
    state: ProgramState = ProgramState.INSTANTIATED
    decisions: Tuple[ArchitectureDecision, ...] = field(default_factory=tuple)
    artifact_hash: Optional[str] = None
    def transition_to(self, ns): return DiscoveryProgram(self.program_id, ns, self.decisions, self.artifact_hash)
    def seal(self, h): return DiscoveryProgram(self.program_id, ProgramState.SEALED, self.decisions, h)

class ProgramRepository:
    _storage: Dict[str, DiscoveryProgram] = {}
    @classmethod
    def get_or_create(cls, pid=None):
        if pid and pid in cls._storage: return cls._storage[pid]
        p = DiscoveryProgram(program_id=pid) if pid else DiscoveryProgram()
        cls._storage[p.program_id] = p; return p
    @classmethod
    def persist(cls, p): cls._storage[p.program_id] = p
