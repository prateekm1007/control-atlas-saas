import numpy as np
from dataclasses import dataclass
from typing import List, Optional
import re

from ..generation.afdb import detect_structure_format

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

# N8 (Round-11, v2 — audit-corrected default): alternate-conformer
# selection POLICY.
#   per_residue  — DEFAULT (Round-11 audit direction: "per-residue
#                  default"): choose ONE altloc ID per residue (highest
#                  TOTAL occupancy, first-seen on ties), keep blank-altloc
#                  (shared) atoms plus the chosen conformer's atoms. This
#                  prevents chimera residues that mix atom sets from
#                  different conformers.
#   per_atom     — opt-in legacy behavior (the original P0 item 6 rule):
#                  for each (chain, residue, atom_name) keep the
#                  highest-occupancy record; ties keep the first in file
#                  order. Still available explicitly via
#                  altloc_policy=ALTLOC_POLICY_PER_ATOM for byte-level
#                  comparison with pre-N8 deliveries.
# v1 of this patch shipped per_residue as opt-in with per_atom kept as
# the engine default; the Round-11 audit rejected that (N8 default
# error). v2 flips the module default to per_residue as directed. The
# flip is an audit-directed semantic change and is recorded as such in
# the delivery register (not a silent change).
ALTLOC_POLICY_PER_ATOM = "per_atom"
ALTLOC_POLICY_PER_RESIDUE = "per_residue"
ALTLOC_POLICY = ALTLOC_POLICY_PER_RESIDUE


def _select_conformers(records: List["Atom"], policy: str) -> List["Atom"]:
    """Apply the alternate-conformer selection policy to parsed records.

    per_residue (N8 v2, DEFAULT per the Round-11 audit direction): first
    choose one altloc ID per residue — the ID whose atoms carry the
    highest TOTAL occupancy, ties broken by first appearance in file
    order — then keep blank-altloc atoms (shared atoms) and the chosen
    ID's atoms, applying the per-atom rule within the survivors as a
    safety net.

    per_atom (opt-in legacy, the original P0 item 6 rule): for each
    (chain, res_seq, icode, atom_name) keep the highest-occupancy
    record; ties keep the FIRST record in file order. Byte-for-byte the
    pre-N8 behavior, available explicitly for comparison.

    Both policies guarantee exactly one atom per (chain, residue,
    atom_name).
    """
    if policy == ALTLOC_POLICY_PER_RESIDUE:
        totals: dict = {}            # (residue key, altloc) -> total occupancy
        for rec in records:
            alt = rec.altloc
            if alt in ("", " "):
                continue  # shared atoms belong to every conformer
            rk = (rec.chain_id, rec.res_seq, rec.insertion_code)
            k = (rk, alt)
            totals[k] = totals.get(k, 0.0) + rec.occupancy
        chosen: dict = {}            # residue key -> (altloc, total)
        for (rk, alt), total in totals.items():
            # dict preserves first-seen order: overwrite only on a
            # strictly greater total, so ties keep the FIRST-seen ID.
            if rk not in chosen or total > chosen[rk][1]:
                chosen[rk] = (alt, total)
        records = [
            r for r in records
            if r.altloc in ("", " ") or (
                (r.chain_id, r.res_seq, r.insertion_code) in chosen
                and r.altloc == chosen[(r.chain_id, r.res_seq,
                                        r.insertion_code)][0])
        ]
    elif policy != ALTLOC_POLICY_PER_ATOM:
        raise ValueError(f"unknown altloc policy: {policy!r}")

    atoms: List[Atom] = []
    best_by_key: dict = {}
    for rec in records:
        key = (rec.chain_id, rec.res_seq, rec.insertion_code, rec.atom_name)
        incumbent = best_by_key.get(key)
        if incumbent is None:
            best_by_key[key] = len(atoms)
            atoms.append(rec)
        elif rec.occupancy > atoms[incumbent].occupancy:
            atoms[incumbent] = rec
    return atoms

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
    # OD-07 (auditor, delegated by owner, 2026-10-20): the audit identity
    # is the full 64-hex SHA-256 of the canonical input computed at the
    # API layer (main.py) and combined with the canon/matrix hashes. The
    # ingestion layer carries NO identifier — the previous per-parse
    # uuid4()[:8] value was random, collision-prone, and could masquerade
    # as an audit identity, so it was removed (empty string here).
    audit_id: str
    confidence: ConfidenceProfile
    # P0 item 9: number of distinct polymer residues present in the SOURCE
    # file (ATOM records + retained HETATM polymer records), counted from
    # the raw records before per-atom parse filtering. The coverage
    # denominator uses this so that unparseable coordinate lines or
    # dropped records cannot silently inflate LAW-105. None when the
    # source carried no polymer residues at all.
    source_residue_count: Optional[int] = None

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
    def run(content: bytes, filename: str, task: str, mode: str,
            altloc_policy: str = ALTLOC_POLICY):
        """
        Parses raw structure bytes with Modality Filters. The format is
        detected from the CONTENT (P0 item 13) — mmCIF input is parsed
        with gemmi (maintained parser, not in-house), PDB input by the
        line parser below.
        1. MODEL 1 Selection: Stops at ENDMDL (NMR Safety).
        2. Hydrogen Exclusion: Prevents LAW-130 inflation.
        3. Method Classification: Categorizes data regime.
        4. HETATM Retention (P0 item 5): polymer-associated modified
           residues (MSE, SEC, ...) are kept and mapped to their parent
           amino acid; waters/ligands/ions stay excluded.
        5. Conformer selection (N8 v2, audit-corrected default): ONE
           conformer per RESIDUE by default (highest total occupancy,
           first-seen ID on ties — chimera prevention). The legacy
           per-atom rule (P0 item 6) remains available explicitly via
           altloc_policy="per_atom" for byte-level comparison with
           pre-N8 deliveries.
        """
        try:
            fmt = detect_structure_format(content)
        except ValueError:
            fmt = "pdb"  # legacy behavior: let the line parser reject it
        if fmt == "cif":
            return IngestionProcessor._run_cif(content, altloc_policy)
        return IngestionProcessor._run_pdb(content, altloc_policy)

    @staticmethod
    def _run_cif(content: bytes,
                 altloc_policy: str = ALTLOC_POLICY) -> "Structure":
        """mmCIF ingestion via gemmi (P0 item 13: no in-house CIF parser).

        Reads the ``_atom_site`` table through gemmi's CIF API and maps it
        onto the same Atom/Structure semantics as the PDB path:
          - group_PDB ATOM rows kept; HETATM rows kept only when in
            MODIFIED_RESIDUE_PARENTS and mapped to the parent amino acid;
          - hydrogens/deuteriums excluded;
          - one conformer per atom: highest occupancy, then first;
          - label_seq_id carried on every atom so residue ordering uses
            the authoritative mmCIF label sequence (P0 item 4);
          - auth_* identifiers preferred, label_* fallback;
          - source_residue_count counts polymer + retained-het residues;
          - classification from _exptl.method when present (AFDB files
            carry none -> predicted).
        """
        import gemmi

        doc = gemmi.cif.read_string(content.decode("utf-8", errors="ignore"))
        block = doc.sole_block()

        full_tags = ["group_PDB", "type_symbol", "label_atom_id",
                     "label_alt_id", "label_comp_id", "label_asym_id",
                     "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x",
                     "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv",
                     "auth_seq_id", "auth_asym_id", "pdbx_PDB_model_num"]
        min_tags = ["group_PDB", "label_atom_id", "label_comp_id",
                    "label_asym_id", "label_seq_id", "Cartn_x", "Cartn_y",
                    "Cartn_z", "B_iso_or_equiv"]
        idx: dict = {}
        tab = None
        for tags in (full_tags, min_tags):
            candidate = block.find("_atom_site.", tags)
            if candidate.width() == len(tags) and len(candidate) > 0:
                tab, idx = candidate, {t: i for i, t in enumerate(tags)}
                break
        if tab is None:
            return Structure(
                atoms=[], audit_id="AUD-empty",
                confidence=ConfidenceProfile([], False, "predicted"),
                source_residue_count=None)

        def val(row, tag):
            i = idx.get(tag)
            if i is None or not row.has(i):
                return None
            v = row[i]
            return None if v in (".", "?", None) else v

        # classification from _exptl.method (pair or loop form)
        method = "predicted"
        is_experimental = False
        exptl = block.find_value("_exptl.method")
        if exptl in (None, ".", "?"):
            t = block.find(["_exptl.method"])
            if t.width() == 1 and len(t) > 0:
                exptl = t[0][0]
        if exptl:
            m = str(exptl).upper()
            if "X-RAY" in m:
                method, is_experimental = "xray", True
            elif "NMR" in m:
                method, is_experimental = "nmr", True
            elif "ELECTRON MICROSCOPY" in m or "CRYO-EM" in m:
                method, is_experimental = "cryo_em", True
        resolution = None
        if is_experimental:
            res_val = block.find_value("_refine.ls_d_res_high")
            try:
                resolution = float(res_val) if res_val not in (None, ".", "?") else None
            except (TypeError, ValueError):
                resolution = None

        records: List[Atom] = []
        source_residue_keys = set()

        for row in tab:
            model = val(row, "pdbx_PDB_model_num")
            if model is not None:
                try:
                    if int(model) != 1:
                        continue  # MODEL 1 isolation, as in the PDB path
                except ValueError:
                    continue
            group = (val(row, "group_PDB") or "ATOM").strip().upper()
            comp = (val(row, "label_comp_id") or "").strip().upper()
            if group == "ATOM":
                res_name = comp
            elif group == "HETATM":
                if comp not in MODIFIED_RESIDUE_PARENTS:
                    continue
                res_name = MODIFIED_RESIDUE_PARENTS[comp]
            else:
                continue
            atom_name = (val(row, "label_atom_id") or "").strip()
            if not atom_name:
                continue
            chain_id = (val(row, "auth_asym_id") or val(row, "label_asym_id")
                        or "A").strip() or "A"
            try:
                res_seq = int(val(row, "auth_seq_id") or val(row, "label_seq_id"))
            except (TypeError, ValueError):
                continue
            label_seq_raw = val(row, "label_seq_id")
            try:
                label_seq = int(label_seq_raw) if label_seq_raw is not None else None
            except ValueError:
                label_seq = None
            icode = (val(row, "pdbx_PDB_ins_code") or "").strip()
            try:
                pos = np.array([float(val(row, "Cartn_x")),
                                float(val(row, "Cartn_y")),
                                float(val(row, "Cartn_z"))])
            except (TypeError, ValueError):
                continue
            occ_raw = val(row, "occupancy")
            try:
                occupancy = float(occ_raw) if occ_raw is not None else 1.0
            except ValueError:
                occupancy = 1.0
            try:
                b_iso = float(val(row, "B_iso_or_equiv") or 0.0)
            except ValueError:
                b_iso = 0.0
            element = (val(row, "type_symbol") or "").strip().upper()
            if not element:
                element = gemmi.Element(atom_name).name.upper()
            if element in ("H", "D"):
                continue
            altloc = (val(row, "label_alt_id") or " ").strip() or " "

            source_residue_keys.add((chain_id, res_seq, icode))
            atom_obj = Atom(
                atom_name=atom_name, element=element, pos=pos,
                res_name=res_name, res_seq=res_seq, chain_id=chain_id,
                insertion_code=icode, b_iso=b_iso,
                label_seq_id=label_seq, altloc=altloc, occupancy=occupancy,
            )
            records.append(atom_obj)

        # P0 item 6 / N8: conformer selection per policy (single shared
        # implementation for both parse paths).
        atoms = _select_conformers(records, altloc_policy)
        b_values = [a.b_iso for a in atoms]
        # OD-07: no ingestion-level identifier (see Structure.audit_id).
        return Structure(
            atoms=atoms,
            audit_id="",
            confidence=ConfidenceProfile(b_values, is_experimental, method,
                                         resolution),
            source_residue_count=len(source_residue_keys) if source_residue_keys else None,
        )

    @staticmethod
    def _run_pdb(content: bytes,
                 altloc_policy: str = ALTLOC_POLICY) -> "Structure":
        lines = content.decode(errors="ignore").splitlines()
        records: List[Atom] = []
        source_residue_keys = set()

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

            # P0 item 9: count SOURCE residues from the raw record keys,
            # independent of coordinate parse success, so the coverage
            # denominator cannot be silently shrunk by unparseable lines.
            try:
                _src_chain = line[21:22].strip() or "A"
                _src_seq = int(line[22:26].strip())
                _src_icode = line[26:27].strip()
                source_residue_keys.add((_src_chain, _src_seq, _src_icode))
            except (ValueError, IndexError):
                pass

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

                records.append(atom_obj)
            except (ValueError, IndexError):
                continue

        # P0 item 6 / N8: conformer selection per policy — exactly one
        # atom per (chain, residue, atom_name) under BOTH policies, so
        # every consumer of structure.atoms (LAW-130/155/200 included)
        # keeps seeing one atom per key.
        atoms = _select_conformers(records, altloc_policy)

        # Confidence reflects the RETAINED conformers only.
        b_values = [a.b_iso for a in atoms]
        # OD-07: no ingestion-level identifier (see Structure.audit_id).
        return Structure(
            atoms=atoms,
            audit_id="",
            confidence=ConfidenceProfile(b_values, is_experimental, method, resolution),
            source_residue_count=len(source_residue_keys) if source_residue_keys else None
        )
