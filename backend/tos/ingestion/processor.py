import numpy as np
from dataclasses import dataclass
from typing import List, Optional
import uuid
import re

@dataclass
class Atom:
    """Industrial-grade Atom representation for forensic audit."""
    atom_name: str
    element: str
    pos: np.ndarray
    res_name: str
    res_seq: int
    chain_id: str
    insertion_code: str
    b_iso: float

class ConfidenceProfile:
    """Discriminates method-specific behaviors and carries structural metadata."""
    def __init__(self, values: List[float], is_experimental: bool, 
                 method: str = "predicted", resolution: Optional[float] = None):
        self.mean_plddt = float(np.mean(values)) if values else 0.0
        self.is_experimental = is_experimental
        self.method = method # xray, nmr, cryo_em, predicted
        self.resolution = resolution
        self._values = values

    @property
    def source_type(self) -> str:
        return "B-factor" if self.is_experimental else "pLDDT"

@dataclass
class Structure:
    """Canonical structure container used by the physics engine."""
    atoms: List[Atom]
    audit_id: str
    confidence: ConfidenceProfile

def _classify_method(lines: List[str]) -> str:
    """Classifies acquisition method. Default is predicted."""
    for line in lines:
        if line.startswith("EXPDTA"):
            l = line.upper()
            if "X-RAY" in l: return "xray"
            if "NMR" in l: return "nmr"
            if "ELECTRON MICROSCOPY" in l: return "cryo_em"
        if line.startswith("REMARK") and "ALPHAFOLD" in line.upper():
            return "predicted"
    return "predicted"

def _extract_resolution(lines: List[str]) -> Optional[float]:
    for line in lines:
        if line.startswith("REMARK   2 RESOLUTION."):
            match = re.search(r'RESOLUTION\.\s+(\d+\.\d+)\s+ANGSTROM', line)
            if match: return float(match.group(1))
    return None

class IngestionProcessor:
    @staticmethod
    def run(content: bytes, filename: str, task: str, mode: str):
        """
        Parses raw PDB bytes with Modality Filters.
        1. MODEL 1 Selection: Stops at ENDMDL (NMR Safety).
        2. Hydrogen Exclusion: Prevents LAW-130 inflation.
        3. Method Classification: Categorizes data regime.
        """
        lines = content.decode(errors="ignore").splitlines()
        atoms = []
        b_values = []
        
        # Pre-classify for filtering logic
        method = _classify_method(lines)
        is_experimental = (method != "predicted")
        resolution = _extract_resolution(lines) if is_experimental else None

        for line in lines:
            # 🛡️ MODEL 1 ISOLATION: Stop after first model
            if line.startswith("ENDMDL"):
                break
            
            if not line.startswith("ATOM  "):
                continue

            try:
                atom_name = line[12:16].strip()
                element = line[76:78].strip().upper()
                
                # 🛡️ HYDROGEN FILTERING: NMR/AF3 models often include H/D
                if element in ("H", "D"):
                    continue

                res_name = line[17:20].strip()
                chain_id = line[21:22].strip() or "A"
                res_seq = int(line[22:26].strip())
                icode = line[26:27].strip()

                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                b_iso = float(line[60:66].strip()) if line[60:66].strip() else 0.0

                atom_obj = Atom(
                    atom_name=atom_name, element=element, pos=np.array([x, y, z]),
                    res_name=res_name, res_seq=res_seq, chain_id=chain_id,
                    insertion_code=icode, b_iso=b_iso
                )
                atoms.append(atom_obj)
                b_values.append(b_iso)
            except (ValueError, IndexError):
                continue

        audit_id = str(uuid.uuid4())[:8]
        return Structure(
            atoms=atoms,
            audit_id=f"AUD-{audit_id}",
            confidence=ConfidenceProfile(b_values, is_experimental, method, resolution)
        )
