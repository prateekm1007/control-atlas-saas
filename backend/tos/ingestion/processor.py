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
    """Discriminates pLDDT (AlphaFold) from B-factor (experimental).
    Carries resolution for sigma scaling."""

    def __init__(self, values: List[float], is_experimental: bool,
                 resolution: Optional[float] = None):
        self.mean_plddt = float(np.mean(values)) if values else 0.0
        self.is_experimental = is_experimental
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


def _detect_experimental(lines: List[str], b_values: List[float]) -> bool:
    """Priority cascade: EXPDTA > AlphaFold REMARK > B-factor heuristic."""
    for line in lines:
        if line.startswith("EXPDTA"):
            return True
        if line.startswith("REMARK") and "ALPHAFOLD" in line.upper():
            return False
        if line.startswith("REMARK") and "DBREF" in line and "AF-" in line:
            return False

    if not b_values:
        return False
    if float(np.max(b_values)) > 100.0:
        return True
    return False


def _extract_resolution(lines: List[str]) -> Optional[float]:
    """
    Extract crystallographic resolution from PDB header.

    Canonical format (wwPDB standard):
        REMARK   2 RESOLUTION.    1.74 ANGSTROMS.

    Returns resolution in Angstroms, or None if not found.
    """
    for line in lines:
        if line.startswith("REMARK   2 RESOLUTION."):
            match = re.search(r'RESOLUTION\.\s+(\d+\.\d+)\s+ANGSTROM', line)
            if match:
                return float(match.group(1))
    return None


class IngestionProcessor:
    @staticmethod
    def run(content: bytes, filename: str, task: str, mode: str):
        """
        Parses raw PDB bytes using strict column offsets.
        ATOM records only. Extracts resolution for experimental structures.
        """
        lines = content.decode(errors="ignore").splitlines()
        atoms = []
        b_values = []

        for line in lines:
            if not line.startswith("ATOM  "):
                continue

            try:
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22].strip() or "A"
                res_seq = int(line[22:26].strip())
                icode = line[26:27].strip()

                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                b_iso_str = line[60:66].strip()
                b_iso = float(b_iso_str) if b_iso_str else 0.0

                element = line[76:78].strip() if len(line) > 77 else ""

                atom_obj = Atom(
                    atom_name=atom_name,
                    element=element,
                    pos=np.array([x, y, z]),
                    res_name=res_name,
                    res_seq=res_seq,
                    chain_id=chain_id,
                    insertion_code=icode,
                    b_iso=b_iso,
                )
                atoms.append(atom_obj)
                b_values.append(b_iso)
            except (ValueError, IndexError):
                continue

        is_experimental = _detect_experimental(lines, b_values)
        resolution = _extract_resolution(lines) if is_experimental else None
        audit_id = str(uuid.uuid4())[:8]

        structure_obj = Structure(
            atoms=atoms,
            audit_id=f"AUD-{audit_id}",
            confidence=ConfidenceProfile(b_values, is_experimental, resolution),
        )

        return structure_obj
