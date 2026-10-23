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

from ..ingestion.atom_selection import (
    atom_identity,
    choose_atom,
    normalize_insertion_code,
)


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


def angle_deg_valid(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> tuple[float, bool]:
    """Return `(angle_deg, valid)` where degenerate vectors are not valid."""
    angle = angle_deg(a, b, c)
    ba = a - b
    bc = c - b
    valid = bool(np.linalg.norm(ba) >= 1e-10 and np.linalg.norm(bc) >= 1e-10)
    return angle, valid


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


def dihedral_deg_valid(a: np.ndarray, b: np.ndarray,
                       c: np.ndarray, d: np.ndarray) -> tuple[float, bool]:
    """Return `(dihedral_deg, valid)` where degenerate geometry is not valid."""
    angle = dihedral_deg(a, b, c, d)
    b1 = b - a
    b2 = c - b
    b3 = d - c
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    valid = bool(
        np.linalg.norm(n1) >= 1e-10
        and np.linalg.norm(n2) >= 1e-10
        and np.linalg.norm(b2) >= 1e-10
    )
    return angle, valid


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
    best = None
    insertion_code = normalize_insertion_code(insertion_code)
    for atom in atoms:
        if (atom.chain_id == chain_id and
                atom.res_seq == res_seq and
                atom.atom_name == atom_name and
                normalize_insertion_code(atom.insertion_code) == insertion_code):
            best = choose_atom(best, atom)
    return best.pos if best is not None else None


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

    Order follows the authoritative sequence position (P1 ordering contract):
    a valid mmCIF label_seq_id wins when present; otherwise the numeric residue
    sequence number is used. The sort is STABLE, so residues sharing a primary
    position (e.g. 184A before 184 in 1PPE) keep their source-file order, which
    is the depositor's chemical order. Alphabetical insertion-code sorting is
    intentionally NOT used, because it can re-pair peptide geometry (P0 item 4).

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
        key = atom_identity(atom)[:4]
        if key not in seen:
            seen.add(key)
            ls = getattr(atom, "label_seq_id", None)
            try:
                ls_int = int(ls) if ls is not None else None
            except (TypeError, ValueError):
                ls_int = None
            residues.append(
                (atom.res_seq, normalize_insertion_code(atom.insertion_code),
                 atom.res_name, ls_int)
            )

    residues.sort(key=lambda r: (r[3] if r[3] is not None else r[0],))
    return [(r[0], r[1], r[2]) for r in residues]


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
    selected = {}
    insertion_code = normalize_insertion_code(insertion_code)
    for atom in atoms:
        if (atom.chain_id == chain_id and
                atom.res_seq == res_seq and
                normalize_insertion_code(atom.insertion_code) == insertion_code):
            selected[atom.atom_name] = choose_atom(selected.get(atom.atom_name), atom)
    return {name: atom.pos for name, atom in selected.items()}


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
