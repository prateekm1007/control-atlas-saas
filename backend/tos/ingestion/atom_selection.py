"""Deterministic atom identity and alternate-conformer selection."""

from __future__ import annotations

import math
from typing import Iterable, Optional


def normalize_altloc(value) -> str:
    text = str(value or "").strip()
    if text in {"", "?", "."} or text == "\x00":
        return ""
    return text


def normalize_insertion_code(value) -> str:
    text = str(value or "").strip()
    if text in {"?", "."} or text == "\x00":
        return ""
    return text


def finite_float(value, default: Optional[float] = None) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return number if math.isfinite(number) else default


def atom_altloc(atom) -> str:
    return normalize_altloc(getattr(atom, "altloc", ""))


def atom_occupancy(atom) -> float:
    value = finite_float(getattr(atom, "occupancy", 1.0), None)
    return value if value is not None else float("-inf")


def atom_is_hetatm(atom) -> bool:
    return bool(getattr(atom, "is_hetatm", False))


def atom_position(atom) -> tuple[float, float, float]:
    pos = getattr(atom, "pos", None)
    return (
        round(float(pos[0]), 3),
        round(float(pos[1]), 3),
        round(float(pos[2]), 3),
    )


def atom_quality_key(atom):
    altloc = atom_altloc(atom)
    return (
        -atom_occupancy(atom),
        0 if altloc == "" else 1,
        altloc,
        0 if not atom_is_hetatm(atom) else 1,
        atom_position(atom),
    )


def atom_identity(atom) -> tuple:
    return (
        str(getattr(atom, "chain_id", "")),
        int(getattr(atom, "res_seq", 0)),
        normalize_insertion_code(getattr(atom, "insertion_code", "")),
        str(getattr(atom, "res_name", "")).strip().upper(),
        str(getattr(atom, "atom_name", "")).strip(),
    )


def residue_identity(atom) -> tuple:
    return atom_identity(atom)[:4]


def atom_sort_key(atom):
    return atom_identity(atom) + (atom_altloc(atom), 0 if not atom_is_hetatm(atom) else 1)


def choose_atom(current, candidate):
    if current is None:
        return candidate
    return min((current, candidate), key=atom_quality_key)


def canonicalize_atoms(atoms: Iterable) -> list:
    # Group atoms by residue identity so altloc selection is residue-level.
    # This prevents constructing an artificial residue by mixing atoms from
    # incompatible alternate conformations.
    residues = {}
    for atom in atoms:
        res_key = residue_identity(atom)
        residues.setdefault(res_key, []).append(atom)

    selected = {}
    for res_key, res_atoms in residues.items():
        # Determine the dominant specific altloc (non-blank) for this residue
        # by total occupancy.  Blank altloc atoms are universal (present in
        # all conformers) and are not counted toward any conformer's total.
        specific_altloc_occ = {}
        for atom in res_atoms:
            alt = atom_altloc(atom)
            if alt == "":
                continue
            occ = atom_occupancy(atom)
            if occ == float("-inf"):
                occ = 0.0
            specific_altloc_occ[alt] = specific_altloc_occ.get(alt, 0.0) + occ

        # Determine dominant conformer among EXPLICIT (non-blank) altlocs
        # only.  Blank-altloc atoms are shared/universal and never compete
        # as a conformer.  If no explicit altlocs exist, dominant is blank.
        if not specific_altloc_occ:
            dominant = ""
        else:
            max_specific = max(specific_altloc_occ.values())
            dominant_candidates = [
                alt for alt, occ in specific_altloc_occ.items()
                if occ == max_specific
            ]
            dominant = sorted(dominant_candidates)[0]

        # Select atoms: prefer the dominant conformer when it has the atom;
        # fall back to blank (universal) atoms otherwise.
        atom_names = {a.atom_name for a in res_atoms}
        for name in atom_names:
            dominant_atoms = [
                a for a in res_atoms
                if a.atom_name == name and atom_altloc(a) == dominant
            ]
            blank_atoms = [
                a for a in res_atoms
                if a.atom_name == name and atom_altloc(a) == ""
            ]
            if dominant and dominant_atoms:
                # Use dominant conformer atom with highest occupancy
                best = max(dominant_atoms, key=atom_occupancy)
            elif blank_atoms:
                # Fall back to shared (blank) atom
                best = max(blank_atoms, key=atom_occupancy)
            else:
                # Neither dominant conformer nor shared atom has this atom.
                # Do NOT substitute from a non-dominant conformer — that would
                # mix conformations.  Skip this atom entirely.
                best = None
            if best is not None:
                key = atom_identity(best)
                selected[key] = choose_atom(selected.get(key), best)
    return list(selected.values())
