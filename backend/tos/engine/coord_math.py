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

# Maximum C(i)-N(i+1) distance for a genuine peptide bond (Engh & Huber
# ideal 1.329 A; the 2.0 A bound matches the LAW-110 backbone-gap
# criterion). Used to validate that a file-order residue pair is truly
# bonded before measuring peptide-plane geometry.
PEPTIDE_BOND_MAX_ANGSTROM = 2.0


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

    P0 item 6: when several records match (alternate conformers), the
    retained atom is the one with the HIGHEST occupancy, ties resolved to
    the FIRST record in file order — the same rule the parser and
    _extract apply. Previously this returned the first match blindly.

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
    for atom in atoms:
        if (atom.chain_id == chain_id and
                atom.res_seq == res_seq and
                atom.atom_name == atom_name and
                atom.insertion_code == insertion_code):
            if best is None or atom.occupancy > best.occupancy:
                best = atom
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


def residues_are_adjacent(chain_a: str, seq_a: int, icode_a: str,
                          chain_b: str, seq_b: int, icode_b: str) -> bool:
    """
    No-missing-span predicate for file-order-consecutive residues
    (P0 item 4, 2026-10-19 audit).

    Callers iterate residues in FILE ORDER within one chain (see
    get_sequential_residues / Tier1Measurements._extract), which is the
    depositor's chemical bond order. This predicate answers a narrower
    question: is there a MISSING RESIDUE SPAN between two
    file-order-consecutive residues, i.e. numbering that jumps by more
    than one base number (45 -> 67)?

    Alphabetical insertion-code comparison is intentionally REMOVED.
    Insertion codes cannot be ordered by string comparison because real
    deposits list insertion-coded residues either side of their base
    number: 1PPE chain E file order is 183, 184A, 184, 185 (sha256
    fe51b428..., evidence/item-4/). Only the file order is ground truth;
    pair-level chemical validation is done with C-N bond distances
    (peptide_bonded).

    Decision table:
      - chain mismatch            -> False (never join chains)
      - seq_b - seq_a > 1         -> False (missing span: absent density,
                                    not a physical chain break; caller
                                    skips instead of measuring across)
      - anything else (seq_b - seq_a in {..., -1, 0, 1}, any icodes)
                                  -> True  (no missing span; measure the
                                    pair — geometry decides pass/break)

    Args:
        chain_a, seq_a, icode_a: chain ID, residue sequence number and
            insertion code of the first residue (icode "" when absent).
        chain_b, seq_b, icode_b: same for the second residue.

    Returns:
        True when the pair is file-order-consecutive with no missing
        residue span between them.
    """
    if chain_a != chain_b:
        return False
    if seq_b - seq_a > 1:
        return False
    return True


def peptide_bonded(atoms_r1: dict, atoms_r2: dict) -> bool:
    """
    Chemical validation that residue 2 follows residue 1 (P0 item 4).

    True when the C(i)-N(i+1) distance is within the peptide-bond
    tolerance (PEPTIDE_BOND_MAX_ANGSTROM). Used by the peptide-plane
    measurements (LAW-120 O-C-N / CA-C-N, LAW-125 phi/psi, LAW-135
    omega) so that angles are only ever measured across a genuine bond,
    never across a discontinuity — the same class of defect as the
    repaired intra-residue pseudo-angle.

    Args:
        atoms_r1: atom_name -> position dict of residue i (needs "C").
        atoms_r2: atom_name -> position dict of residue i+1 (needs "N").

    Returns:
        True when C(i)-N(i+1) is a peptide bond. Missing C or N atoms
        return False (cannot validate -> do not measure).
    """
    if "C" not in atoms_r1 or "N" not in atoms_r2:
        return False
    return distance(atoms_r1["C"], atoms_r2["N"]) <= PEPTIDE_BOND_MAX_ANGSTROM


def get_sequential_residues(atoms, chain_id: str) -> List[Tuple[int, str, str]]:
    """
    Return an ordered list of residues in a chain, in FILE ORDER.

    Each entry is (res_seq, insertion_code, res_name).

    P0 item 4: residues are returned in order of first appearance in the
    source file, NOT sorted by (res_seq, insertion_code). Real deposits
    list insertion-coded residues before their base number (1PPE chain E:
    183, 184A, 184, 185), so any (res_seq, icode) sort silently re-pairs
    the peptide bonds and fabricates backbone gaps. mmCIF inputs ordered
    by label_seq_id set Atom.label_seq_id, which takes precedence when
    present.

    Args:
        atoms: list of Atom objects
        chain_id: chain to extract

    Returns:
        Ordered list of (res_seq, insertion_code, res_name) tuples.
        Duplicates removed (multiple atoms per residue collapsed to one entry).
    """
    seen = set()
    residues = []
    pending_label_seq: dict = {}
    for atom in atoms:
        if atom.chain_id != chain_id:
            continue
        key = (atom.res_seq, atom.insertion_code)
        if key not in seen:
            seen.add(key)
            residues.append((atom.res_seq, atom.insertion_code, atom.res_name))
            label = getattr(atom, "label_seq_id", None)
            if label is not None:
                pending_label_seq[key] = label

    if pending_label_seq:
        # mmCIF: label_seq_id is the authoritative sequence position.
        residues.sort(key=lambda r: pending_label_seq.get(
            (r[0], r[1]), float("inf")))
    return residues


def get_residue_atoms(atoms, chain_id: str, res_seq: int,
                      insertion_code: str = "") -> dict:
    """
    Return a dictionary of atom_name -> position for all atoms in a residue.

    P0 item 6: alternate conformers are collapsed to one atom per name —
    highest occupancy, then first in file order. Previously this kept the
    LAST matching record, disagreeing with get_atom_pos (first match).

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
    for atom in atoms:
        if (atom.chain_id == chain_id and
                atom.res_seq == res_seq and
                atom.insertion_code == insertion_code):
            incumbent = selected.get(atom.atom_name)
            if incumbent is None or atom.occupancy > incumbent.occupancy:
                selected[atom.atom_name] = atom
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
