import numpy as np
from dataclasses import dataclass
from typing import List, Optional
import uuid
import re

# P0 item 5 (2026-10-19 audit): HETATM records that are polymer
# constituents are RETAINED and mapped to their parent amino acid.
# Selenomethionine is the canonical case: it is a normal Met side chain
# with Se replacing S, deposited as HETATM MSE — dropping it silently
# deletes residues from the chain and fabricates backbone gaps.
# Waters, ligands, ions, and nucleotides have no amino-acid parent and
# remain excluded. D-amino-acid variants (DAL, DLE, ...) are
# deliberately NOT mapped: LAW-145 must be able to flag them.
MODIFIED_RESIDUE_PARENTS = {
    "MSE": "MET",  # selenomethionine
    "SEC": "CYS",  # selenocysteine
    "PYL": "LYS",  # pyrrolysine
    "HYP": "PRO",  # 4-hydroxyproline (collagen)
    "CME": "CYS", "CSO": "CYS", "CSD": "CYS",  # modified cysteines
    "OCS": "CYS", "CSX": "CYS", "CAS": "CYS", "CSM": "CYS",
    "MLY": "LYS", "M3L": "LYS", "KCX": "LYS",  # modified lysines
    "SEP": "SER", "TPO": "THR", "PTR": "TYR",  # phosphorylated
    "PCA": "GLU", "NLE": "LEU",  # pyroglutamate, norleucine
}

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
    # P0 item 4: mmCIF label_seq_id — the authoritative sequence position
    # for residue ordering when the source is mmCIF (None for PDB sources,
    # where file order is authoritative). See coord_math.get_sequential_residues.
    label_seq_id: Optional[int] = None
    # P0 item 6: alternate-conformer identifier (PDB column 17) and
    # occupancy (columns 55-60). One conformer per atom is retained —
    # highest occupancy, then first in file — everywhere in the engine.
    altloc: str = " "
    occupancy: float = 1.0

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
        4. HETATM Retention (P0 item 5): polymer-associated modified
           residues (MSE, SEC, ...) are kept and mapped to their parent
           amino acid; waters/ligands/ions stay excluded.
        """
        lines = content.decode(errors="ignore").splitlines()
        atoms = []
        best_by_key = {}

        # Pre-classify for filtering logic
        method = _classify_method(lines)
        is_experimental = (method != "predicted")
        resolution = _extract_resolution(lines) if is_experimental else None

        for line in lines:
            # 🛡️ MODEL 1 ISOLATION: Stop after first model
            if line.startswith("ENDMDL"):
                break

            record_type = line[:6]
            if record_type == "ATOM  ":
                res_name = line[17:20].strip()
            elif record_type == "HETATM":
                # P0 item 5: keep only polymer-associated modified
                # residues, mapped to their parent amino acid.
                het_name = line[17:20].strip()
                if het_name not in MODIFIED_RESIDUE_PARENTS:
                    continue
                res_name = MODIFIED_RESIDUE_PARENTS[het_name]
            else:
                continue

            try:
                atom_name = line[12:16].strip()
                element = line[76:78].strip().upper()

                # 🛡️ HYDROGEN FILTERING: NMR/AF3 models often include H/D
                if element in ("H", "D"):
                    continue

                altloc = line[16:17] or " "
                chain_id = line[21:22].strip() or "A"
                res_seq = int(line[22:26].strip())
                icode = line[26:27].strip()

                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                occ_field = line[54:60].strip()
                occupancy = float(occ_field) if occ_field else 1.0
                b_iso = float(line[60:66].strip()) if line[60:66].strip() else 0.0

                atom_obj = Atom(
                    atom_name=atom_name, element=element, pos=np.array([x, y, z]),
                    res_name=res_name, res_seq=res_seq, chain_id=chain_id,
                    insertion_code=icode, b_iso=b_iso,
                    altloc=altloc, occupancy=occupancy
                )

                # P0 item 6: exactly one conformer per atom — highest
                # occupancy wins, ties keep the FIRST record in file order.
                # Selection happens at parse time so that every consumer
                # of structure.atoms (LAW-130/155/200 included) sees one
                # atom per (chain, residue, atom_name).
                key = (chain_id, res_seq, icode, atom_name)
                incumbent = best_by_key.get(key)
                if incumbent is None:
                    best_by_key[key] = len(atoms)
                    atoms.append(atom_obj)
                elif atom_obj.occupancy > atoms[incumbent].occupancy:
                    atoms[incumbent] = atom_obj
            except (ValueError, IndexError):
                continue

        # Confidence reflects the RETAINED conformers only.
        b_values = [a.b_iso for a in atoms]
        audit_id = str(uuid.uuid4())[:8]
        return Structure(
            atoms=atoms,
            audit_id=f"AUD-{audit_id}",
            confidence=ConfidenceProfile(b_values, is_experimental, method, resolution)
        )
