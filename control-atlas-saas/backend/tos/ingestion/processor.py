import hashlib, os, tempfile, numpy as np
from Bio.PDB import PDBParser, MMCIFParser
from tos.ingestion.structure_object import StructureObject, Atom, ConfidenceSidecar

class IngestionProcessor:
    @staticmethod
    def run(content: bytes, filename: str, gen: str) -> StructureObject:
        ext = filename.split(".")[-1].lower()
        if ext not in ["pdb", "cif"]: ext = "pdb"
        with tempfile.NamedTemporaryFile(suffix=f".{ext}", mode='wb', delete=False) as tf:
            tf.write(content); t_path = tf.name
        parser = PDBParser(QUIET=True) if ext=="pdb" else MMCIFParser(QUIET=True)
        struct = parser.get_structure("CLAIM", t_path); os.remove(t_path)
        atoms, ligands, plddts = [], [], []
        model = struct[0]
        for chain in model:
            for residue in chain:
                is_lig = residue.get_id()[0].strip() != ""
                for atom in residue:
                    if atom.get_altloc() not in [" ", "A"]: continue
                    if atom.element == "H": continue
                    val = float(atom.get_bfactor())
                    plddts.append(val)
                    a_obj = Atom(residue.get_resname(), residue.get_id()[1], chain.id, atom.get_name(), atom.element, tuple(float(x) for x in atom.get_coord()), val)
                    if is_lig: ligands.append(a_obj)
                    else: atoms.append(a_obj)
        p_arr = np.array(plddts)
        mean_p = float(np.mean(p_arr)) if len(p_arr) > 0 else 0.0
        is_placeholder = len(p_arr) > 0 and (np.all(p_arr == 0) or np.all(p_arr == 100.0) or np.all(p_arr == 1.0))
        if len(p_arr) == 0: sidecar = ConfidenceSidecar(0.0, "ABSENT", "NONE", False)
        elif is_placeholder: sidecar = ConfidenceSidecar(mean_p, "PLACEHOLDER", "BFACTOR_STATIC", False)
        else:
            if 0 < mean_p <= 1.0: mean_p *= 100
            sidecar = ConfidenceSidecar(mean_p, "MEASURED", "BFACTOR_COLUMN", (mean_p < 70.0))
        return StructureObject(hashlib.sha256(content).hexdigest()[:16], gen, ext.upper(), tuple(atoms), tuple(ligands), sidecar)
