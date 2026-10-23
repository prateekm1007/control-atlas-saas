"""OD-07 canonical coordinate identity and audit identifiers."""

from __future__ import annotations

import hashlib
import json
import math

from ..ingestion.atom_selection import (
    atom_altloc,
    atom_identity,
    atom_is_hetatm,
    atom_position,
    normalize_insertion_code,
)

COORDINATE_IDENTITY_SCHEMA = "od07-coordinate-identity/1.0.0"
AUDIT_ID_SCHEMA = "od07-audit-id/1.0.0"


def _canonical_json(value) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")


def _coordinate_record(atom) -> dict:
    position = atom_position(atom)
    if not all(math.isfinite(value) for value in position):
        raise ValueError("non-finite coordinates cannot produce an audit identity")
    chain_id, res_seq, insertion_code, res_name, atom_name = atom_identity(atom)
    return {
        "atom_name": atom_name,
        "chain_id": chain_id,
        "element": str(getattr(atom, "element", "")).strip().upper(),
        "insertion_code": normalize_insertion_code(insertion_code),
        "is_hetatm": atom_is_hetatm(atom),
        "res_name": res_name,
        "res_seq": res_seq,
        "altloc": atom_altloc(atom),
        "x": position[0],
        "y": position[1],
        "z": position[2],
    }


def compute_coordinate_hash(structure_or_atoms) -> str:
    atoms = getattr(structure_or_atoms, "atoms", structure_or_atoms)
    records = [_coordinate_record(atom) for atom in atoms]
    records.sort(
        key=lambda record: (
            record["chain_id"],
            record["res_seq"],
            record["insertion_code"],
            record["res_name"],
            record["atom_name"],
            record["altloc"],
            record["is_hetatm"],
        )
    )
    payload = {
        "schema": COORDINATE_IDENTITY_SCHEMA,
        "atoms": records,
    }
    return hashlib.sha256(_canonical_json(payload)).hexdigest()


def compute_audit_id(coordinate_hash: str, canon_hash: str, matrix_hash: str) -> str:
    for name, value in (
        ("coordinate_hash", coordinate_hash),
        ("canon_hash", canon_hash),
        ("matrix_hash", matrix_hash),
    ):
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{name} must be a non-empty canonical digest")
    payload = {
        "schema": AUDIT_ID_SCHEMA,
        "coordinate_hash": coordinate_hash.strip().lower(),
        "canon_hash": canon_hash.strip().lower(),
        "matrix_hash": matrix_hash.strip().lower(),
    }
    return hashlib.sha256(_canonical_json(payload)).hexdigest()
