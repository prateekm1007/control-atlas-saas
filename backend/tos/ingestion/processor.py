import numpy as np
from dataclasses import dataclass
from typing import List
import uuid

@dataclass
class Atom:
    atom_name: str
    element: str
    pos: np.ndarray
    res_name: str
    res_seq: int
    chain_id: str
    insertion_code: str
    b_iso: float

@dataclass
class Structure:
    atoms: List[Atom]
    audit_id: str
    confidence: any

class IngestionProcessor:
    @staticmethod
    def run(content: bytes, filename: str, task: str, mode: str):
        lines = content.decode(errors="ignore").splitlines()
        atoms = []
        plddts = []
        
        for line in lines:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                try:
                    # PRO-GRADE COLUMN EXTRACTION
                    name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    chain = line[21:22].strip() or "A"
                    res_seq = int(line[22:26].strip())
                    icode = line[26:27].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    # pLDDT is in B-factor col (61-66)
                    b_iso = float(line[60:66].strip())
                    elem = line[76:78].strip()
                    
                    atoms.append(Atom(
                        atom_name=name, element=elem, pos=np.array([x,y,z]),
                        res_name=res_name, res_seq=res_seq, chain_id=chain,
                        insertion_code=icode, b_iso=b_iso
                    ))
                    plddts.append(b_iso)
                except: continue
        
        class Conf: pass
        c = Conf()
        c.mean_plddt = np.mean(plddts) if plddts else 0.0
        return Structure(atoms=atoms, audit_id=str(uuid.uuid4())[:8], confidence=c)
