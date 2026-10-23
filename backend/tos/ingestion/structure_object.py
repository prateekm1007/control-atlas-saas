from dataclasses import dataclass, field
from typing import Tuple, Dict

from ..engine.identity import compute_coordinate_hash


@dataclass(frozen=True)
class Atom:
    res_name: str
    res_seq: int
    atom_name: str
    pos: Tuple[float, float, float]
    element: str = "C"
    chain_id: str = ""
    insertion_code: str = ""
    b_iso: float = 0.0
    altloc: str = ""
    occupancy: float = 1.0
    is_hetatm: bool = False


@dataclass(frozen=True)
class ConfidenceSidecar:
    mean_plddt: float
    per_residue_plddt: Dict[int, float] = field(default_factory=dict)


@dataclass(frozen=True)
class StructureObject:
    audit_id: str
    atoms: Tuple[Atom, ...]
    confidence: ConfidenceSidecar
    label: str = "Unknown"

    def get_coordinate_hash(self) -> str:
        return compute_coordinate_hash(self.atoms)
