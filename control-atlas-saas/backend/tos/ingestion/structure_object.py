from dataclasses import dataclass, field
from typing import Tuple, Dict, Any, Optional
import numpy as np
import hashlib

@dataclass(frozen=True)
class Atom:
    res_name: str
    res_seq: int
    chain: str
    atom_name: str
    element: str
    pos: Tuple[float, float, float]
    plddt: Optional[float] = None

@dataclass(frozen=True)
class ConfidenceSidecar:
    mean_plddt: float
    status: str      # MEASURED, ABSENT, or PLACEHOLDER
    source: str      # BFACTOR_COLUMN, MA_QA_METRIC, etc.
    is_low: bool

@dataclass(frozen=True)
class StructureObject:
    audit_id: str
    source_model: str
    source_format: str
    atoms: Tuple[Atom, ...]
    ligands: Tuple[Atom, ...]
    confidence: ConfidenceSidecar
    metadata: Dict[str, Any] = field(default_factory=dict)

    def get_coordinate_hash(self) -> str:
        all_content = sorted(list(self.atoms) + list(self.ligands), 
                             key=lambda a: (a.chain, a.res_seq, a.atom_name))
        coord_string = "".join([f"{a.chain}{a.res_seq}{a.atom_name}{a.pos[0]:.3f}{a.pos[1]:.3f}{a.pos[2]:.3f}" 
                               for a in all_content])
        return hashlib.sha256(coord_string.encode()).hexdigest()
