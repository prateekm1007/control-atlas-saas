from dataclasses import dataclass, field
from typing import Dict, List, Optional
import math
import re

import gemmi
import numpy as np

from .atom_selection import (
    canonicalize_atoms,
    finite_float,
    normalize_altloc,
    normalize_insertion_code,
    residue_identity,
)


class IngestionError(ValueError):
    """Raised when a source cannot produce a canonical structure."""


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
    altloc: str = ""
    occupancy: float = 1.0
    is_hetatm: bool = False
    label_asym_id: str = ""
    label_seq_id: Optional[int] = None
    model_index: int = 1
    # 0-based index of the source record (chain/residue iteration order) as
    # read from the PDB/mmCIF file. This is the authoritative depositor order
    # for peptide adjacency; alphabetical insertion-code sorting must not
    # override it (P0 item 4 / P1 ordering contract).
    source_order: int = -1


@dataclass
class IngestionMetadata:
    format: str = "unknown"
    model_policy: str = "first_model_only"
    altloc_policy: str = "highest_occupancy_per_atom; tie-break blank then lexicographic"
    polymer_hetatm_policy: str = (
        "include polymer HETATM and backbone-complete non-water/non-ligand HETATM; "
        "exclude waters, non-polymer ligands, and other HETATM"
    )
    source_atom_records: int = 0
    selected_atoms: int = 0
    total_residues: int = 0
    complete_backbone_residues: int = 0
    excluded_hydrogens: int = 0
    excluded_waters: int = 0
    excluded_ligands: int = 0
    excluded_other: int = 0
    excluded_models: int = 0
    malformed_records: int = 0
    altloc_conflicts: int = 0
    altloc_selected_atoms: int = 0
    missing_occupancy: int = 0
    missing_b_factor: int = 0
    # Source polymer-sequence coverage metadata (P2 source-backed coverage).
    # per chain id -> number of residues in the authoritative polymer sequence
    # (PDB SEQRES / mmCIF entity_poly_seq) for that chain's entity.
    source_polymer_residues_by_chain: Dict[str, int] = field(default_factory=dict)
    source_sequence_available: bool = False
    warnings: List[str] = field(default_factory=list)


class ConfidenceProfile:
    """Discriminates method-specific behaviors and carries structural metadata."""
    def __init__(self, values: List[float], is_experimental: bool,
                 method: str = "predicted", resolution: Optional[float] = None):
        finite_values = [
            number for number in (finite_float(value, None) for value in values)
            if number is not None
        ]
        self.is_experimental = is_experimental
        self.method = method
        self.resolution = resolution
        if is_experimental:
            self._values = []
            self._bfactors = finite_values
            self.mean_plddt = None
        else:
            self._values = finite_values
            self._bfactors = []
            self.mean_plddt = float(np.mean(finite_values)) if finite_values else 0.0

    @property
    def source_type(self) -> str:
        if self._values:
            return "pLDDT"
        if self._bfactors:
            return "B-factor"
        return "unavailable"


@dataclass
class Structure:
    """Canonical structure container used by the physics engine."""
    atoms: List[Atom]
    audit_id: str
    confidence: ConfidenceProfile
    metadata: IngestionMetadata = field(default_factory=IngestionMetadata)


def _source_text(content: bytes) -> str:
    return content.decode("utf-8-sig", errors="replace")


def _detect_format(content: bytes, filename: str) -> str:
    lower = (filename or "").lower()
    text = _source_text(content).lstrip()
    if lower.endswith((".cif", ".mmcif")):
        return "mmcif"
    if lower.endswith((".pdb", ".ent")):
        return "pdb"
    if text.startswith("data_") or "_atom_site." in text[:5000]:
        return "mmcif"
    if any(
        line.startswith(("ATOM  ", "HETATM", "HEADER", "MODEL", "REMARK"))
        for line in text.splitlines()[:20]
    ):
        return "pdb"
    raise IngestionError("unrecognized structure format; expected PDB or mmCIF")


def _read_gemmi(content: bytes, source_format: str):
    text = _source_text(content)
    try:
        if source_format == "mmcif":
            document = gemmi.cif.read_string(text)
            block = document.sole_block()
            structure = gemmi.make_structure_from_block(block)
        else:
            structure = gemmi.read_pdb_string(text)
    except Exception as exc:
        raise IngestionError(f"{source_format} parsing failed: {exc}") from exc
    if len(structure) == 0 or not structure[0]:
        raise IngestionError(f"{source_format} source contains no models or atoms")
    return structure


def _classify_method(info: dict, text: str) -> str:
    """Classify acquisition method.

    Returns one of: "xray", "nmr", "cryo_em", "predicted", "unknown".
    Unknown provenance must NOT silently default to "predicted" — that would
    cause a real experimental structure to be judged by the predicted cohort
    (pLDDT gate, pLDDT coverage basis) rather than the experimental one.
    """
    raw = str(info.get("_exptl.method", "") or "").upper()
    if "X-RAY" in raw:
        return "xray"
    if "NMR" in raw:
        return "nmr"
    if "ELECTRON MICROSCOPY" in raw or "CRYO" in raw:
        return "cryo_em"
    upper = text.upper()
    if "ALPHAFOLD" in upper or "PREDICTED" in raw or "COMPUTATIONAL" in raw:
        return "predicted"
    # No acquisition method was declared in the source.  The caller must
    # decide whether this structure is experimental or predicted; treating
    # it as predicted here would silently downgrade a real crystal structure.
    return "unknown"


def _extract_resolution(structure, text: str) -> Optional[float]:
    resolution = finite_float(getattr(structure, "resolution", None), None)
    if resolution is not None and resolution > 0.0:
        return resolution
    for line in text.splitlines():
        if line.startswith("REMARK   2 RESOLUTION."):
            match = re.search(r"RESOLUTION\.\s+(\d+\.\d+)\s+ANGSTROM", line)
            if match:
                return float(match.group(1))
    return None


def _infer_element(atom_name: str) -> str:
    letters = re.sub(r"[^A-Za-z]", "", atom_name or "").upper()
    return letters[:1] or "X"


def _is_polymer_like(residue) -> bool:
    if residue.is_water():
        return False
    if residue.entity_type == gemmi.EntityType.Polymer:
        return True
    if residue.het_flag == "A":
        return True
    names = {atom.name.strip().upper() for atom in residue}
    return {"N", "CA", "C"}.issubset(names)


def _normalize_chain(value) -> str:
    text = str(value or "").strip()
    return "" if text in {"?", "."} else text


def _source_polymer_positions(structure, model):
    """Return (by_chain, available).

    by_chain maps a normalized chain id to the number of residues in the
    authoritative polymer sequence for that chain's entity (PDB SEQRES or
    mmCIF entity_poly_seq). This count may exceed the number of residues
    observed in the coordinate records (unresolved/absent positions), which
    is exactly the source-backed coverage denominator (P2).

    available is False when a polymer chain's source sequence cannot be
    established (no entity, unknown polymer type, or empty full_sequence).
    """
    by_chain: Dict[str, int] = {}
    available = True
    for chain in model:
        entity = None
        # Prefer the chain's polymer span (ResidueSpan) for get_entity_of.
        try:
            polymer_span = chain.get_polymer()
        except Exception:
            polymer_span = None
        if polymer_span is not None and len(polymer_span) > 0:
            try:
                entity = structure.get_entity_of(polymer_span)
            except Exception:
                entity = None
            has_polymer = True
        else:
            has_polymer = False
        if not has_polymer or entity is None:
            continue
        seq = getattr(entity, "full_sequence", None)
        if not seq:
            available = False
            continue
        cid = _normalize_chain(chain.name)
        by_chain[cid] = len(seq)
    return by_chain, available


def _select_model(structure):
    return structure[0]


class IngestionProcessor:
    @staticmethod
    def run(content: bytes, filename: str, task: str, mode: str):
        source_format = _detect_format(content, filename)
        text = _source_text(content)
        structure = _read_gemmi(content, source_format)
        try:
            structure.setup_entities()
        except Exception:
            pass

        info = dict(structure.info.items())
        method = _classify_method(info, text)
        is_experimental = method != "predicted"
        resolution = _extract_resolution(structure, text) if is_experimental else None

        metadata = IngestionMetadata(
            format=source_format,
            excluded_models=max(0, len(structure) - 1),
        )
        model = _select_model(structure)
        selected = []
        identities = {}
        b_values = []
        source_order = 0

        for chain in model:
            chain_id = _normalize_chain(chain.name)
            for residue in chain:
                atom_records = list(residue)
                metadata.source_atom_records += len(atom_records)
                if residue.is_water():
                    metadata.excluded_waters += len(atom_records)
                    continue
                if residue.entity_type in (gemmi.EntityType.NonPolymer, gemmi.EntityType.Branched):
                    metadata.excluded_ligands += len(atom_records)
                    continue
                if not _is_polymer_like(residue):
                    metadata.excluded_other += len(atom_records)
                    continue

                res_name = residue.name.strip().upper()
                res_seq = int(residue.seqid.num)
                insertion_code = normalize_insertion_code(residue.seqid.icode)
                label_seq_id = residue.label_seq if residue.label_seq is not None else None
                is_hetatm = residue.het_flag == "H"

                for atom in atom_records:
                    atom_name = atom.name.strip()
                    element = atom.element.name.strip().upper()
                    if not element or element == "X":
                        element = _infer_element(atom_name)
                    if element in {"H", "D"}:
                        metadata.excluded_hydrogens += 1
                        continue

                    try:
                        pos = np.array(
                            [float(atom.pos.x), float(atom.pos.y), float(atom.pos.z)],
                            dtype=float,
                        )
                    except (TypeError, ValueError):
                        metadata.malformed_records += 1
                        continue
                    if pos.shape != (3,) or not np.all(np.isfinite(pos)):
                        metadata.malformed_records += 1
                        continue

                    occupancy = finite_float(atom.occ, None)
                    if occupancy is None:
                        metadata.missing_occupancy += 1
                        occupancy = float("nan")
                    b_iso = finite_float(atom.b_iso, None)
                    if b_iso is None:
                        metadata.missing_b_factor += 1
                        b_iso = float("nan")

                    altloc = normalize_altloc(atom.altloc)
                    atom_obj = Atom(
                        atom_name=atom_name,
                        element=element,
                        pos=pos,
                        res_name=res_name,
                        res_seq=res_seq,
                        chain_id=chain_id,
                        insertion_code=insertion_code,
                        b_iso=b_iso,
                        altloc=altloc,
                        occupancy=occupancy,
                        is_hetatm=is_hetatm,
                        label_asym_id=_normalize_chain(getattr(chain, "name", "")),
                        label_seq_id=int(label_seq_id) if label_seq_id is not None else None,
                        model_index=1,
                        source_order=source_order,
                    )
                    selected.append(atom_obj)
                    source_order += 1
                    if finite_float(b_iso, None) is not None:
                        b_values.append(float(b_iso))

        canonical_atoms = canonicalize_atoms(selected)
        identity_counts = {}
        identity_altlocs = {}
        for atom in selected:
            key = (
                atom.chain_id,
                atom.res_seq,
                atom.insertion_code,
                atom.res_name,
                atom.atom_name,
            )
            identity_counts[key] = identity_counts.get(key, 0) + 1
            identity_altlocs.setdefault(key, set()).add(atom.altloc)
        metadata.altloc_conflicts = sum(
            1 for key, count in identity_counts.items()
            if count > 1 or len(identity_altlocs[key]) > 1
        )
        metadata.altloc_selected_atoms = sum(
            1 for atom in canonical_atoms if atom.altloc
        )
        metadata.selected_atoms = len(canonical_atoms)

        residues = {}
        for atom in canonical_atoms:
            residues.setdefault(residue_identity(atom), set()).add(atom.atom_name)
        metadata.total_residues = len(residues)
        metadata.complete_backbone_residues = sum(
            1 for names in residues.values() if {"N", "CA", "C"}.issubset(names)
        )
        src_by_chain, src_available = _source_polymer_positions(structure, model)
        metadata.source_polymer_residues_by_chain = src_by_chain
        metadata.source_sequence_available = src_available
        if src_available and src_by_chain:
            metadata.warnings.append(
                "coverage denominator uses the source polymer sequence "
                f"({sum(src_by_chain.values())} positions across "
                f"{len(src_by_chain)} chain(s))"
            )
        elif not src_available:
            metadata.warnings.append(
                "source polymer sequence metadata unavailable; LAW-105 coverage "
                "is reported on observed coordinates only (P2)"
            )
        if metadata.malformed_records:
            metadata.warnings.append(
                f"{metadata.malformed_records} atom record(s) had non-finite coordinates and were discarded"
            )
        if metadata.missing_b_factor:
            metadata.warnings.append(
                f"{metadata.missing_b_factor} atom record(s) had no finite B-factor"
            )
        if not canonical_atoms:
            raise IngestionError(
                f"{source_format} source produced no canonical polymer atoms"
            )

        return Structure(
            atoms=canonical_atoms,
            audit_id="",
            confidence=ConfidenceProfile(b_values, is_experimental, method, resolution),
            metadata=metadata,
        )
