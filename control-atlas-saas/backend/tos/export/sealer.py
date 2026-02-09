import hashlib

class ForensicSealer:
    @staticmethod
    def canonical_serialize(atoms_list) -> str:
        sorted_atoms = sorted(atoms_list, key=lambda a: (a.chain, a.res_seq, a.atom_name))
        return "".join([f"{a.chain}{a.res_seq}{a.atom_name}{a.pos[0]:.3f}{a.pos[1]:.3f}{a.pos[2]:.3f}" 
                        for a in sorted_atoms])

    @staticmethod
    def generate_hash(canonical_str: str) -> str:
        return hashlib.sha256(canonical_str.encode()).hexdigest()

    @staticmethod
    def seal_structure(content_bytes: bytes, audit_id: str, verdict: str, coord_hash: str, ext: str) -> str:
        """Sovereign Sealing: Handles PDB (REMARK) and CIF (Comment) safely."""
        # Safe decoding for metadata injection
        body = content_bytes.decode('utf-8', errors='ignore')
        
        prefix = "REMARK" if ext.lower() == "pdb" else "#"
        header = [
            f"{prefix} 900 TOSCANINI FORENSIC SEAL v14.2",
            f"{prefix} 900 AUDIT_ID: {audit_id}",
            f"{prefix} 901 VERDICT: {verdict}",
            f"{prefix} 902 COORD_HASH: {coord_hash}",
            f"{prefix} 903 STATUS: NOTARIZED_SOVEREIGN",
            f"{prefix} 903 " + ("="*30)
        ]
        return "\n".join(header) + "\n" + body
