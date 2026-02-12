import gemmi, hashlib, logging
from .structure_object import StructureObject, Atom, ConfidenceSidecar
from ..utils.type_guards import force_str

logger = logging.getLogger("toscanini.ingestion")

class IngestionProcessor:
    @staticmethod
    def run(content: bytes, filename: str, label: str, mode: str) -> StructureObject:
        audit_id = hashlib.sha256(content).hexdigest()[:12]
        atoms = []
        try:
            raw_str = force_str(content)
            if filename.lower().endswith(('.cif', '.mmcif')):
                structure = gemmi.make_structure(gemmi.cif.read_string(raw_str).sole_block())
            else:
                structure = gemmi.read_pdb_string(raw_str)
            
            plddts = {}
            for model in structure:
                for chain in model:
                    for residue in chain:
                        res_idx = int(residue.seqid.num)
                        for atom in residue:
                            atoms.append(Atom(residue.name, res_idx, atom.name, 
                                            (atom.pos.x, atom.pos.y, atom.pos.z), atom.element.name))
                            if atom.name == "CA":
                                plddts[res_idx] = atom.b_iso
            
            mean_plddt = sum(plddts.values()) / len(plddts) if plddts else 50.0
            return StructureObject(audit_id, tuple(atoms), ConfidenceSidecar(round(mean_plddt, 1), plddts), label)
        except Exception as e:
            logger.error(f"Ingestion Fault: {e}")
            return StructureObject(audit_id, (Atom("UNK", 1, "CA", (0.0, 0.0, 0.0), "C"),), 
                                 ConfidenceSidecar(0.0), label)
