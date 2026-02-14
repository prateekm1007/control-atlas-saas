import numpy as np
from dataclasses import dataclass
from typing import List
import uuid

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
    b_iso: float # This holds the critical pLDDT / B-factor data

@dataclass
class Structure:
    """Canonical structure container used by the physics engine."""
    atoms: List[Atom]
    audit_id: str
    confidence: any

class IngestionProcessor:
    @staticmethod
    def run(content: bytes, filename: str, task: str, mode: str):
        """
        Parses raw PDB bytes using strict column offsets.
        Ensures AlphaFold pLDDT is captured from the B-factor column.
        """
        lines = content.decode(errors="ignore").splitlines()
        atoms = []
        plddt_values = []
        
        for line in lines:
            # Only process standard ATOM and HETATM records
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                try:
                    # PDB Format Standard Column Offsets
                    # Reference: https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    chain_id = line[21:22].strip() or "A"
                    res_seq = int(line[22:26].strip())
                    icode = line[26:27].strip() # Insertion Code
                    
                    # Coordinates
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    
                    # Column 61-66: Temperature Factor / B-Factor (holds pLDDT in AF models)
                    b_iso_str = line[60:66].strip()
                    b_iso = float(b_iso_str) if b_iso_str else 0.0
                    
                    # Column 77-78: Element symbol
                    element = line[76:78].strip()
                    
                    atom_obj = Atom(
                        atom_name=atom_name,
                        element=element,
                        pos=np.array([x, y, z]),
                        res_name=res_name,
                        res_seq=res_seq,
                        chain_id=chain_id,
                        insertion_code=icode,
                        b_iso=b_iso
                    )
                    atoms.append(atom_obj)
                    plddt_values.append(b_iso)
                except (ValueError, IndexError):
                    # Skip malformed lines, but don't crash the whole ingestion
                    continue
        
        # Calculate mean pLDDT for the strategic confidence object
        class ConfidenceProfile:
            def __init__(self, values):
                self.mean_plddt = np.mean(values) if values else 0.0
        
        audit_id = str(uuid.uuid4())[:8]
        structure_obj = Structure(
            atoms=atoms,
            audit_id=f"AUD-{audit_id}",
            confidence=ConfidenceProfile(plddt_values)
        )
        
        return structure_obj
