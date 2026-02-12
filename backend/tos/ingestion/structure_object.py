from dataclasses import dataclass, field
from typing import Tuple, Dict
import hashlib
from ..utils.type_guards import force_bytes

@dataclass(frozen=True)
class Atom:
    res_name: str; res_seq: int; atom_name: str; pos: Tuple[float, float, float]; element: str = "C"

@dataclass(frozen=True)
class ConfidenceSidecar:
    mean_plddt: float
    per_residue_plddt: Dict[int, float] = field(default_factory=dict)

@dataclass(frozen=True)
class StructureObject:
    audit_id: str; atoms: Tuple[Atom, ...]; confidence: ConfidenceSidecar; label: str = "Unknown"
    def get_coordinate_hash(self) -> str:
        coord_string = "|".join([f"{a.atom_name}:{a.pos}" for a in self.atoms])
        return hashlib.sha256(force_bytes(coord_string)).hexdigest()
