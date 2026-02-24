"""
TOSCANINI Coordinate Math Library (Module 1.1)

Pure geometric primitives for structural measurements.
No governance logic. No thresholds. No verdicts. No side effects.

All functions operate on numpy arrays of shape (3,) representing
3D atomic positions in Angstroms.

References:
- Engh & Huber (1991), Acta Cryst. A47, 392-400
- Ramachandran & Sasisekharan (1968), Adv. Protein Chem. 23
"""
import numpy as np
from typing import List, Optional, Tuple, Set


def distance(a: np.ndarray, b: np.ndarray) -> float:
    """
    Euclidean distance between two 3D points.

    Args:
        a: position vector (3,)
        b: position vector (3,)

    Returns:
        Distance in Angstroms.
    """
    return float(np.linalg.norm(a - b))


def angle_deg(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """
    Angle at vertex B formed by vectors BA and BC, in degrees.

    Args:
        a: position of atom A
        b: position of atom B (vertex)
        c: position of atom C

    Returns:
        Angle in degrees [0, 180].
        Returns 0.0 if any vector has zero length (degenerate geometry).
    """
    ba = a - b
    bc = c - b

    norm_ba = np.linalg.norm(ba)
    norm_bc = np.linalg.norm(bc)

    if norm_ba < 1e-10 or norm_bc < 1e-10:
        return 0.0

    cos_angle = np.dot(ba, bc) / (norm_ba * norm_bc)
    # Clamp to [-1, 1] to handle floating point edge cases
    cos_angle = np.clip(cos_angle, -1.0, 1.0)

    return float(np.degrees(np.arccos(cos_angle)))


def dihedral_deg(a: np.ndarray, b: np.ndarray,
                 c: np.ndarray, d: np.ndarray) -> float:
    """
    Dihedral (torsion) angle defined by four points, in degrees.

    Computes the angle between the planes (A,B,C) and (B,C,D).

    Used for:
    - Phi:   C(i-1) - N(i) - CA(i) - C(i)
    - Psi:   N(i) - CA(i) - C(i) - N(i+1)
    - Omega: CA(i) - C(i) - N(i+1) - CA(i+1)
    - Chi1:  N - CA - CB - CG

    Args:
        a, b, c, d: position vectors of the four atoms

    Returns:
        Dihedral angle in degrees [-180, 180].
        Returns 0.0 if geometry is degenerate (collinear atoms).
    """
    b1 = b - a
    b2 = c - b
    b3 = d - c

    # Normal vectors to the two planes
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    norm_n1 = np.linalg.norm(n1)
    norm_n2 = np.linalg.norm(n2)

    # Degenerate case: collinear atoms produce zero-length cross products
    if norm_n1 < 1e-10 or norm_n2 < 1e-10:
        return 0.0

    n1 = n1 / norm_n1
    n2 = n2 / norm_n2

    # Unit vector along b2
    b2_norm = np.linalg.norm(b2)
    if b2_norm < 1e-10:
        return 0.0
    m1 = np.cross(n1, b2 / b2_norm)

    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    return float(np.degrees(np.arctan2(-y, x)))


def get_atom_pos(atoms, chain_id: str, res_seq: int,
                 atom_name: str, insertion_code: str = "") -> Optional[np.ndarray]:
    """
    Locate a specific atom and return its position.

    Args:
        atoms: list of Atom objects (from IngestionProcessor)
        chain_id: chain identifier (e.g., "A")
        res_seq: residue sequence number
        atom_name: PDB atom name (e.g., "CA", "N", "C", "CB", "SG")
        insertion_code: PDB insertion code (default "")

    Returns:
        numpy array of shape (3,) with atom position, or None if not found.
    """
    for atom in atoms:
        if (atom.chain_id == chain_id and
                atom.res_seq == res_seq and
                atom.atom_name == atom_name and
                atom.insertion_code == insertion_code):
            return atom.pos
    return None


def get_chains(atoms) -> Set[str]:
    """
    Return the set of unique chain IDs in the structure.

    Args:
        atoms: list of Atom objects

    Returns:
        Set of chain ID strings.
    """
    return {atom.chain_id for atom in atoms}


def get_sequential_residues(atoms, chain_id: str) -> List[Tuple[int, str, str]]:
    """
    Return an ordered list of residues in a chain.

    Each entry is (res_seq, insertion_code, res_name).
    Sorted by (res_seq, insertion_code) to handle insertion codes correctly.

    Args:
        atoms: list of Atom objects
        chain_id: chain to extract

    Returns:
        Ordered list of (res_seq, insertion_code, res_name) tuples.
        Duplicates removed (multiple atoms per residue collapsed to one entry).
    """
    seen = set()
    residues = []
    for atom in atoms:
        if atom.chain_id != chain_id:
            continue
        key = (atom.res_seq, atom.insertion_code)
        if key not in seen:
            seen.add(key)
            residues.append((atom.res_seq, atom.insertion_code, atom.res_name))

    # Sort by sequence number, then insertion code
    residues.sort(key=lambda r: (r[0], r[1]))
    return residues


def get_residue_atoms(atoms, chain_id: str, res_seq: int,
                      insertion_code: str = "") -> dict:
    """
    Return a dictionary of atom_name -> position for all atoms in a residue.

    Args:
        atoms: list of Atom objects
        chain_id: chain identifier
        res_seq: residue sequence number
        insertion_code: PDB insertion code

    Returns:
        dict mapping atom name strings to numpy position arrays.
        Example: {"N": array([1,2,3]), "CA": array([4,5,6]), ...}
    """
    result = {}
    for atom in atoms:
        if (atom.chain_id == chain_id and
                atom.res_seq == res_seq and
                atom.insertion_code == insertion_code):
            result[atom.atom_name] = atom.pos
    return result


# ═══════════════════════════════════════════════════════════════
# CHI1 ATOM DEFINITIONS
# Maps residue name to the four atoms defining the chi1 dihedral.
# Residues without chi1 (GLY, ALA) are excluded.
# ═══════════════════════════════════════════════════════════════

CHI1_ATOMS = {
    "ARG": ("N", "CA", "CB", "CG"),
    "ASN": ("N", "CA", "CB", "CG"),
    "ASP": ("N", "CA", "CB", "CG"),
    "CYS": ("N", "CA", "CB", "SG"),
    "GLN": ("N", "CA", "CB", "CG"),
    "GLU": ("N", "CA", "CB", "CG"),
    "HIS": ("N", "CA", "CB", "CG"),
    "ILE": ("N", "CA", "CB", "CG1"),
    "LEU": ("N", "CA", "CB", "CG"),
    "LYS": ("N", "CA", "CB", "CG"),
    "MET": ("N", "CA", "CB", "CG"),
    "PHE": ("N", "CA", "CB", "CG"),
    "PRO": ("N", "CA", "CB", "CG"),
    "SER": ("N", "CA", "CB", "OG"),
    "THR": ("N", "CA", "CB", "OG1"),
    "TRP": ("N", "CA", "CB", "CG"),
    "TYR": ("N", "CA", "CB", "CG"),
    "VAL": ("N", "CA", "CB", "CG1"),
}

# Residues that do not have a chi1 angle
NO_CHI1_RESIDUES = {"GLY", "ALA"}
