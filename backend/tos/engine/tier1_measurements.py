"""
TOSCANINI Tier-1 Physics Engine (v22.5.3 — Module 1.2.1)

Canonical math migration complete. All geometric calculations
delegate to coord_math.py. No local math functions.

Changes from pre-migration:
  - Removed: _dist, _angle, _dihedral (local duplicates)
  - LAW-100: uses coord_math.distance
  - LAW-110: uses coord_math.distance
  - LAW-125: uses coord_math.dihedral_deg
  - LAW-120: IMPLEMENTED (Bond Angle RMSD via coord_math.angle_deg)
  - LAW-160: sample field fixed to report actual CA pairs measured
  - LAW-155: explicit result inserted (was missing entirely)
  - Zero-observation laws report NOT_MEASURED or NOT_APPLICABLE, never PASS.
  - LAW-100 uses the canonical percentage-of-outliers threshold in every modality.
"""
from .coord_math import (
    distance, angle_deg_valid, dihedral_deg_valid,
    get_atom_pos, get_chains, get_sequential_residues,
    CHI1_ATOMS, NO_CHI1_RESIDUES
)
from ..ingestion.atom_selection import (
    canonicalize_atoms,
    finite_float,
    normalize_insertion_code,
)
import math
import numpy as np
import logging
from ..governance.station_sop import (
    LAW_CANON, SIGMA_TABLE, IDEAL_TABLE,
    STANDARD_RESIDUES, HYDROPHOBIC_RESIDUES
)

# Phase A.5: Residue-level granularity constants
GRANULARITY_RESIDUE      = "residue"
GRANULARITY_RESIDUE_PAIR = "residue_pair"
GRANULARITY_GLOBAL       = "global"

logger = logging.getLogger("toscanini.tier1")

MEASURED_STATUSES = frozenset({"PASS", "FAIL", "VETO"})
MISSING_STATUSES = frozenset({"NOT_MEASURED", "MEASUREMENT_FAILED", "NOT_APPLICABLE"})
VALID_MEASUREMENT_STATUSES = MEASURED_STATUSES | MISSING_STATUSES

INSERTION_ALPHABET = (
    "",
    *tuple("ABCDEFGHIJKLMNOPQRSTUVWXYZ"),
    *tuple("abcdefghijklmnopqrstuvwxyz"),
)
INSERTION_RANK = {code: index for index, code in enumerate(INSERTION_ALPHABET)}
REQUIRED_BACKBONE_ATOMS = frozenset({"N", "CA", "C"})
PEPTIDE_BOND_MIN = 0.5
PEPTIDE_BOND_MAX = 2.0


def _residue_label_sort_key(label: str):
    chain_id, _, residue = label.partition(":")
    digits = "".join(character for character in residue if character.isdigit())
    return chain_id, int(digits) if digits else 0, residue


class Tier1Measurements:
    @staticmethod
    def _not_measured(reason):
        return {
            "observed": None,
            "status": "NOT_MEASURED",
            "sample": 0,
            "reason": reason,
        }

    @staticmethod
    def _not_applicable(reason):
        return {
            "observed": None,
            "status": "NOT_APPLICABLE",
            "sample": 0,
            "reason": reason,
        }

    @staticmethod
    def _numeric_adjacent(r1, r2):
        """Numbering compatibility (PDB residue number + insertion code).

        Used as a guard on top of source-order consecutiveness so that two
        residues separated by a genuinely missing position are not paired.
        Insertions (same residue number, distinct insertion codes) are adjacent
        when consecutive in source order; a step across a numbered position
        requires the canonical blank/insertion transition.
        """
        s1, s2 = int(r1["_seq"]), int(r2["_seq"])
        i1 = normalize_insertion_code(r1.get("_icode", ""))
        i2 = normalize_insertion_code(r2.get("_icode", ""))
        if s2 == s1:
            return i1 != i2
        if s2 == s1 + 1:
            return (not i1 and not i2) or (bool(i1) and not i2)
        return False

    @staticmethod
    def _is_sequential(r1, r2):
        """Adjacency defined by the canonical residue-order contract.

        Two residues are sequential iff they share a chain, ``r2`` is the next
        canonical position after ``r1`` (``_pos`` + 1), and — for coordinate-
        derived bases — they are numbering-compatible. For the sequence-position
        basis (label_seq_id) adjacency is fully encoded by the positions
        themselves. This keeps ordering and adjacency on one contract while
        never pairing residues separated by a genuinely missing position.
        """
        if r1["_chain"] != r2["_chain"]:
            return False
        p1 = r1.get("_pos")
        p2 = r2.get("_pos")
        if p1 is None or p2 is None:
            return False
        if p2 != p1 + 1:
            return False
        if r1.get("_basis") == "label_seq_id":
            return True
        return Tier1Measurements._numeric_adjacent(r1, r2)

    @staticmethod
    def _has_peptide_bond(r1, r2):
        if not Tier1Measurements._is_sequential(r1, r2):
            return False
        if "C" not in r1["_atoms"] or "N" not in r2["_atoms"]:
            return False
        length = distance(r1["_atoms"]["C"], r2["_atoms"]["N"])
        return math.isfinite(length) and PEPTIDE_BOND_MIN <= length <= PEPTIDE_BOND_MAX

    @staticmethod
    def _is_rama_outlier(res_name, phi, psi):
        if phi is None or psi is None:
            return False
        if res_name == "GLY":
            return (-20 < phi < 20) and (-20 < psi < 20)
        if res_name == "PRO":
            return not (-100 < phi < -30 and -50 < psi < 180)
        alpha = (-180 < phi < 0) and (-90 < psi < 50)
        beta = (-180 < phi < -20) and (20 < psi < 180)
        left = (20 < phi < 120) and (-60 < psi < 80)
        bridge = (-180 < phi < -20) and (-10 < psi < 70)
        return not (alpha or beta or left or bridge)

    @staticmethod
    def _extract(structure):
        res_map = {}
        for a in canonicalize_atoms(getattr(structure, "atoms", []) or []):
            k = (
                a.chain_id,
                int(a.res_seq),
                normalize_insertion_code(a.insertion_code),
                str(a.res_name).strip().upper(),
            )
            label_seq_id = getattr(a, "label_seq_id", None)
            label_seq_id = int(label_seq_id) if isinstance(label_seq_id, int) and label_seq_id >= 0 else None
            source_order = getattr(a, "source_order", None)
            source_order = int(source_order) if isinstance(source_order, int) and source_order >= 0 else None
            if k not in res_map:
                res_map[k] = {
                    "_name": str(a.res_name).strip().upper(), "_chain": a.chain_id,
                    "_seq": int(a.res_seq), "_icode": normalize_insertion_code(a.insertion_code),
                    "_lseq": label_seq_id, "_src_raw": source_order,
                    "_conf_acc": [], "_atoms": {}
                }
            else:
                if res_map[k]["_lseq"] is None:
                    res_map[k]["_lseq"] = label_seq_id
                if res_map[k]["_src_raw"] is None:
                    res_map[k]["_src_raw"] = source_order
            res_map[k]["_atoms"][a.atom_name] = a.pos
            b_value = finite_float(getattr(a, "b_iso", None), None)
            if b_value is not None:
                res_map[k]["_conf_acc"].append(b_value)
        for r in res_map.values():
            r["_conf"] = float(np.mean(r["_conf_acc"])) if r["_conf_acc"] else None
        return Tier1Measurements._order_residues(list(res_map.values()))

    @staticmethod
    def _order_residues(residues):
        """Order residues per chain under ONE consistent evidence basis.

        Basis precedence per chain: valid ``label_seq_id`` (gap-sensitive),
        else ``source_order`` (deposited file order), else residue number.
        Both ordering and adjacency (via ``_is_sequential``) derive from the
        resulting ``_pos``, so the order contract and the adjacency contract
        cannot diverge.
        """
        by_chain = {}
        for r in residues:
            by_chain.setdefault(r["_chain"], []).append(r)
        out = []
        for chain in sorted(by_chain):
            members = by_chain[chain]
            if all(r["_lseq"] is not None for r in members):
                ordered = sorted(members, key=lambda r: r["_lseq"])
                for r in ordered:
                    r["_pos"] = r["_lseq"]
                    r["_basis"] = "label_seq_id"
                    r["_gap_known"] = True
            elif all(r["_src_raw"] is not None for r in members):
                ordered = sorted(members, key=lambda r: r["_src_raw"])
                for index, r in enumerate(ordered):
                    r["_pos"] = index
                    r["_basis"] = "source_order"
                    r["_gap_known"] = False
            else:
                ordered = sorted(members, key=lambda r: (r["_seq"], r["_icode"]))
                for index, r in enumerate(ordered):
                    r["_pos"] = index
                    r["_basis"] = "residue_number"
                    r["_gap_known"] = False
            out.extend(ordered)
        return out

    @staticmethod
    def _plddt_advisory(residues, is_experimental):
        if is_experimental or not residues:
            return None
        scores = [r["_conf"] for r in residues if r["_conf"] is not None]
        if not scores:
            return None
        above = sum(1 for s in scores if s >= 70.0)
        return {
            "mean_plddt": round(float(np.mean(scores)), 2),
            "fraction_at_or_above_70": round(above / len(scores) * 100, 2),
            "observation_count": len(scores),
            "policy": "advisory only; does not affect status or verdict",
        }

    @staticmethod
    def run_full_audit(structure, user_intent="NONE"):
        method = getattr(structure.confidence, "method", "predicted")
        is_experimental = (method != "predicted")

        for atom in getattr(structure, "atoms", []) or []:
            try:
                pos = np.asarray(getattr(atom, "pos", None), dtype=float)
            except (TypeError, ValueError):
                pos = None
            if pos is None or pos.shape != (3,) or not np.all(np.isfinite(pos)):
                return {
                    "LAW-105": {
                        "observed": None,
                        "status": "MEASUREMENT_FAILED",
                        "sample": 0,
                        "reason": (
                            "non-finite atomic coordinates cannot be measured; "
                            "no law may pass on malformed geometry"
                        ),
                    }
                }, 0.0, is_experimental

        residues = Tier1Measurements._extract(structure)
        if not residues:
            return {}, 0, is_experimental

        complete = [
            r for r in residues
            if REQUIRED_BACKBONE_ATOMS.issubset(r["_atoms"].keys())
        ]
        if is_experimental:
            core = complete
        else:
            core = [
                r for r in complete
                if r["_conf"] is not None and r["_conf"] >= 70.0
            ]

        metadata = getattr(structure, "metadata", None)

        # Source-backed coverage: when an authoritative source polymer sequence
        # is available (label_seq_id), the denominator is the total number of
        # source polymer positions and the numerator is the set of eligible core
        # residues that map to a source position — so positions absent from the
        # coordinates count against coverage. When the source sequence cannot be
        # established, coverage is computed from observed residues but is
        # explicitly marked non-authoritative (source_backed=False) and is never
        # presented as source-backed.
        source_seq = getattr(metadata, "source_sequences", {}) if metadata is not None else {}
        source_available = bool(source_seq) and all(
            (source_seq.get(r["_chain"]) or {}).get("positions") is not None
            for r in residues
        )
        if source_available:
            coverage_basis = "source_polymer_positions"
            source_backed = True
            denominator_positions = sum(
                entry["positions"] for entry in source_seq.values()
                if entry.get("positions") is not None)
            numerator_positions = sum(
                1 for r in core if r.get("_lseq") is not None)
            coverage = (
                (numerator_positions / denominator_positions * 100)
                if denominator_positions else 0.0
            )
        else:
            if is_experimental:
                coverage_basis = "complete_backbone_residues"
            else:
                coverage_basis = "plddt_core_complete_backbone_residues"
            source_backed = False
            denominator_positions = len(residues)
            numerator_positions = len(core)
            coverage = (
                (numerator_positions / denominator_positions * 100)
                if denominator_positions else 0.0
            )

        results = {}

        # ═══════════════════════════════════════════════════════════
        # LAW-100: Bond Integrity (uses coord_math.distance)
        # Measures z-scores of covalent bond lengths against Engh-Huber ideals.
        # Canonical metric: percentage of bonds with |z| > 4σ. Applies to every
        # modality; the former experimental-only RMSZ comparison was unregistered.
        # ─────────────────────────────────────────────────────────────────────
        z_scores, outliers = [], 0
        for r in core:
            for b in ["N-CA", "CA-C", "C-O", "CA-CB"]:
                if b not in IDEAL_TABLE:
                    continue
                a1_name, a2_name = b.split('-')
                if a1_name in r["_atoms"] and a2_name in r["_atoms"]:
                    d = distance(r["_atoms"][a1_name], r["_atoms"][a2_name])
                    z = abs(d - IDEAL_TABLE[b]) / SIGMA_TABLE.get(b, 0.02)
                    z_scores.append(z)
                    if z > 4.0:
                        outliers += 1

        if not z_scores:
            results["LAW-100"] = Tier1Measurements._not_measured(
                "no measurable covalent bonds in the eligible residues"
            )
        else:
            bond_pct = outliers / len(z_scores) * 100
            threshold = LAW_CANON["LAW-100"]["threshold"]
            results["LAW-100"] = {
                "observed": round(bond_pct, 2),
                "status": "PASS" if bond_pct <= threshold else "VETO",
                "sample": len(z_scores),
            }

        # ═══════════════════════════════════════════════════════════
        # LAW-105: Reliability Coverage
        # ═══════════════════════════════════════════════════════════
        metadata = getattr(structure, "metadata", None)
        law105 = {
            "observed": round(coverage, 2),
            "status": "PASS" if coverage >= 70.0 else "FAIL",
            "sample": len(residues),
            "coverage_basis": coverage_basis,
            "source_backed": source_backed,
            "coverage_numerator": numerator_positions,
            "coverage_denominator": denominator_positions,
            "required_atoms": "N,CA,C",
            "complete_backbone_residues": len(complete),
            "source_atom_records": getattr(metadata, "source_atom_records", None),
            "selected_atoms": getattr(metadata, "selected_atoms", None),
            "malformed_atom_records": getattr(metadata, "malformed_records", None),
        }
        if not source_backed:
            law105["note"] = (
                "source polymer sequence unavailable; coverage is computed from "
                "observed residues and is not source-backed / not authoritative"
            )
        results["LAW-105"] = law105

        # ═══════════════════════════════════════════════════════════
        # LAW-110: Backbone Gap (uses coord_math.distance)
        # Measures C(i)-N(i+1) peptide bond distance for sequential pairs.
        # ═══════════════════════════════════════════════════════════
        gaps = 0
        measured_pairs = 0
        unknown_gap_found = False
        for i in range(len(residues) - 1):
            r1, r2 = residues[i], residues[i + 1]
            if r1["_chain"] != r2["_chain"]:
                continue
            gap_known = bool(r1.get("_gap_known"))
            pos_delta = r2.get("_pos", 0) - r1.get("_pos", 0)
            if gap_known and pos_delta > 1:
                # A genuine source polymer position lies between observed residues;
                # peptide-bond accounting cannot be assumed complete.
                unknown_gap_found = True
                continue
            if not Tier1Measurements._is_sequential(r1, r2):
                continue
            if "C" in r1["_atoms"] and "N" in r2["_atoms"]:
                measured_pairs += 1
                if distance(r1["_atoms"]["C"], r2["_atoms"]["N"]) > 2.0:
                    gaps += 1
        if unknown_gap_found:
            results["LAW-110"] = Tier1Measurements._not_measured(
                "missing polymer positions prevent complete peptide-bond accounting"
            )
        elif measured_pairs == 0:
            results["LAW-110"] = Tier1Measurements._not_measured(
                "no eligible sequential C-N pairs with both atoms present"
            )
        else:
            results["LAW-110"] = {
                "observed": gaps,
                "status": "PASS" if gaps == 0 else "VETO",
                "sample": measured_pairs
            }

        # ═══════════════════════════════════════════════════════════
        # LAW-120: Bond Angle RMSD (Module 1.2.1 — uses coord_math.angle_deg)
        # Measures N-CA-C, O(i)-C(i)-N(i+1), and CA-C-N angles against
        # Engh-Huber ideals. Reports RMSD of angular deviations in degrees.
        # Threshold: 18.0° from station_sop.py
        # ─────────────────────────────────────────────────────────────────────
        angle_sq_diffs = []
        for i in range(len(core)):
            r = core[i]
            atoms = r["_atoms"]

            # Intra-residue: N-CA-C
            if all(a in atoms for a in ["N", "CA", "C"]):
                measured, angle_valid = angle_deg_valid(atoms["N"], atoms["CA"], atoms["C"])
                if not angle_valid:
                    continue
                ideal = IDEAL_TABLE["N-CA-C"]
                angle_sq_diffs.append((measured - ideal) ** 2)

        # Inter-residue: O(i)-C(i)-N(i+1) and CA-C-N across sequential pairs
        for i in range(len(core) - 1):
            r1, r2 = core[i], core[i + 1]
            if not Tier1Measurements._has_peptide_bond(r1, r2):
                continue
            if all(a in r1["_atoms"] for a in ["O", "C"]) and "N" in r2["_atoms"]:
                measured, angle_valid = angle_deg_valid(
                    r1["_atoms"]["O"], r1["_atoms"]["C"], r2["_atoms"]["N"]
                )
                if angle_valid:
                    ideal = IDEAL_TABLE["O-C-N"]
                    angle_sq_diffs.append((measured - ideal) ** 2)
            if "CA" in r1["_atoms"] and "C" in r1["_atoms"] and "N" in r2["_atoms"]:
                measured, angle_valid = angle_deg_valid(
                    r1["_atoms"]["CA"], r1["_atoms"]["C"], r2["_atoms"]["N"]
                )
                if angle_valid:
                    ideal = IDEAL_TABLE["CA-C-N"]
                    angle_sq_diffs.append((measured - ideal) ** 2)

        if not angle_sq_diffs:
            results["LAW-120"] = Tier1Measurements._not_measured(
                "no measurable backbone bond angles in the eligible residues"
            )
        else:
            angle_rmsd = float(np.sqrt(np.mean(angle_sq_diffs)))
            results["LAW-120"] = {
                "observed": round(angle_rmsd, 2),
                "status": "PASS" if angle_rmsd <= LAW_CANON["LAW-120"]["threshold"] else "VETO",
                "sample": len(angle_sq_diffs)
            }

        # ═══════════════════════════════════════════════════════════
        # LAW-125: Ramachandran (uses coord_math.dihedral_deg)
        # Computes phi/psi dihedrals and flags outliers.
        # ═══════════════════════════════════════════════════════════
        rama_out, rama_tot = 0, 0
        rama_outliers = []  # Phase A.5: inline residue capture
        for i in range(1, len(core) - 1):
            p, c, n = core[i - 1], core[i], core[i + 1]
            if not (Tier1Measurements._has_peptide_bond(p, c) and
                    Tier1Measurements._has_peptide_bond(c, n)):
                continue
            if ("C" not in p["_atoms"] or
                    not all(k in c["_atoms"] for k in ["N", "CA", "C"]) or
                    "N" not in n["_atoms"]):
                continue

            phi, phi_valid = dihedral_deg_valid(
                p["_atoms"]["C"], c["_atoms"]["N"],
                c["_atoms"]["CA"], c["_atoms"]["C"]
            )
            psi, psi_valid = dihedral_deg_valid(
                c["_atoms"]["N"], c["_atoms"]["CA"],
                c["_atoms"]["C"], n["_atoms"]["N"]
            )

            if not (phi_valid and psi_valid):
                continue
            rama_tot += 1
            if Tier1Measurements._is_rama_outlier(c["_name"], phi, psi):
                rama_out += 1
                rama_outliers.append(f"{c['_chain']}:{c['_seq']}")  # chain:seq format

        if rama_tot == 0:
            results["LAW-125"] = Tier1Measurements._not_measured(
                "no eligible sequential residue triples with measurable phi/psi"
            )
            results["LAW-125"]["granularity"] = GRANULARITY_RESIDUE
            results["LAW-125"]["failing_residues"] = []
        else:
            rama_pct = rama_out / rama_tot * 100
            results["LAW-125"] = {
                "observed": round(rama_pct, 2),
                "status": "PASS" if rama_pct <= LAW_CANON["LAW-125"]["threshold"] else "VETO",
                "sample": rama_tot,
                "granularity": GRANULARITY_RESIDUE,
                "failing_residues": sorted(rama_outliers, key=lambda x: (x.split(":")[0], int(x.split(":")[1])))
            }

        # ═══════════════════════════════════════════════════════════
        # LAW-130: Clashscore (Neighbor Exclusion)
        # Uses scipy KDTree for efficiency. Geometry via numpy directly
        # (pairwise distance already handled by KDTree radius query).
        # ═══════════════════════════════════════════════════════════
        clash = 0
        clash_pairs = []  # Phase A.5: inline residue-pair capture
        canonical_atoms = canonicalize_atoms(getattr(structure, "atoms", []) or [])
        all_coords = np.array([a.pos for a in canonical_atoms], dtype=float)
        neighbor_search_error = None
        residue_pos = {}
        for r in residues:
            residue_pos[(r["_chain"], r["_seq"], r["_icode"])] = r["_pos"]

        def _is_peptide_neighbor(ai, aj):
            ka = (ai.chain_id, ai.res_seq, ai.insertion_code)
            kb = (aj.chain_id, aj.res_seq, aj.insertion_code)
            pa = residue_pos.get(ka)
            pb = residue_pos.get(kb)
            return pa is not None and pb is not None and pb == pa + 1
        if len(all_coords) > 1:
            try:
                from scipy.spatial import cKDTree
                tree = cKDTree(all_coords)
                pairs = tree.query_pairs(r=1.5)
                for ii, jj in pairs:
                    ai, aj = canonical_atoms[ii], canonical_atoms[jj]
                    same_residue = (
                        ai.chain_id,
                        ai.res_seq,
                        ai.insertion_code,
                    ) == (
                        aj.chain_id,
                        aj.res_seq,
                        aj.insertion_code,
                    )
                    if same_residue:
                        continue
                    if (
                        _is_peptide_neighbor(ai, aj)
                        and ((ai.atom_name, aj.atom_name) == ("C", "N")
                             or (ai.atom_name, aj.atom_name) == ("N", "C"))
                    ):
                        continue
                    clash += 1
                    r1 = f"{ai.chain_id}:{ai.res_seq}{normalize_insertion_code(ai.insertion_code)}"
                    r2 = f"{aj.chain_id}:{aj.res_seq}{normalize_insertion_code(aj.insertion_code)}"
                    pair = tuple(sorted([r1, r2], key=_residue_label_sort_key))
                    clash_pairs.append(pair)
            except Exception as exc:
                neighbor_search_error = f"{type(exc).__name__}: {exc}"
        # Deduplicate and sort pairs
        unique_pairs = sorted(
            set(clash_pairs),
            key=lambda pair: (
                _residue_label_sort_key(pair[0]),
                _residue_label_sort_key(pair[1]),
            ),
        )
        atom_input_count = len(all_coords)
        if len(all_coords) < 2:
            results["LAW-130"] = {
                "observed": None, "status": "NOT_MEASURED", "sample": 0,
                "atom_input_count": atom_input_count,
                "granularity": GRANULARITY_RESIDUE_PAIR,
                "failing_residue_pairs": [],
                "reason": "fewer than two atomic coordinates available",
                "plddt_advisory": Tier1Measurements._plddt_advisory(residues, is_experimental),
            }
        elif neighbor_search_error is not None:
            results["LAW-130"] = {
                "observed": None, "status": "MEASUREMENT_FAILED", "sample": 0,
                "atom_input_count": atom_input_count,
                "granularity": GRANULARITY_RESIDUE_PAIR,
                "failing_residue_pairs": [],
                "reason": f"neighbor search failed: {neighbor_search_error}",
                "plddt_advisory": Tier1Measurements._plddt_advisory(residues, is_experimental),
            }
        else:
            threshold = LAW_CANON["LAW-130"]["threshold"]
            # Precision contract: the verdict is computed from the SAME rounded
            # value that is reported, so the displayed value and status never
            # contradict one another.
            reported = round((clash / atom_input_count) * 1000, 2)
            results["LAW-130"] = {
                "observed": reported,
                "status": "PASS" if reported <= threshold else "VETO",
                "sample": atom_input_count,
                "atom_input_count": atom_input_count,
                "clash_count": clash,
                "precision": "verdict evaluates the reported rounded observed value",
                "granularity": GRANULARITY_RESIDUE_PAIR,
                "failing_residue_pairs": unique_pairs,
                "plddt_advisory": Tier1Measurements._plddt_advisory(residues, is_experimental),
            }

        # ═══════════════════════════════════════════════════════════
        # LAW-170: Residue Identity
        # ═══════════════════════════════════════════════════════════
        non_std = sum(1 for r in residues if r["_name"] not in STANDARD_RESIDUES)
        results["LAW-170"] = {
            "observed": non_std,
            "status": "PASS" if non_std == 0 else "VETO",
            "sample": len(residues)
        }

        # ═══════════════════════════════════════════════════════════
        # Implemented laws follow.
        # ═══════════════════════════════════════════════════════════
        # LAW-145: Chirality (Module 1.2.3 — uses coord_math.dihedral_deg)
        # Checks Cα chirality via improper dihedral N-CA-C-CB.
        # L-amino acids: improper dihedral is negative (~-34°).
        # D-amino acids or errors: improper dihedral is positive.
        # GLY is excluded (no CB atom, no chirality center).
        # Threshold: 0 violations (any violation = VETO).
        # Reference: IUPAC-IUB (1970), Biochemistry 9, 3471-3479
        # ═══════════════════════════════════════════════════════════
        chiral_violations = 0
        chiral_checked = 0
        chiral_candidates = 0
        chiral_unmeasured = 0
        chiral_outlier_residues = []  # Phase A.5: inline residue capture
        for r in core:
            if r["_name"] == "GLY":
                continue
            chiral_candidates += 1
            atoms = r["_atoms"]
            if not all(a in atoms for a in ["N", "CA", "C", "CB"]):
                chiral_unmeasured += 1
                continue
            improper, improper_valid = dihedral_deg_valid(
                atoms["N"], atoms["CA"], atoms["C"], atoms["CB"]
            )
            if not improper_valid:
                chiral_unmeasured += 1
                continue
            chiral_checked += 1
            # Improper dihedral: N-CA-C-CB
            # L-amino acid → negative, D-amino acid → positive
            # Genuine D-amino acid: improper between 0° and 90°
            # Values near ±180° indicate degenerate/coplanar geometry
            # (CB in the N-CA-C plane) — refinement artifact, not chirality error.
            # Reference: 4HHB Chain D Res 78 LEU has improper=179.7° with
            # CB only 0.009Å from the N-CA-C plane.
            if 0.0 < improper < 90.0:
                chiral_violations += 1
                chiral_outlier_residues.append(f"{r['_chain']}:{r['_seq']}")

        if chiral_candidates == 0:
            results["LAW-145"] = Tier1Measurements._not_applicable(
                "no chiral residue centres (glycine-only eligible residues)"
            )
        elif chiral_checked == 0 or chiral_unmeasured > 0:
            results["LAW-145"] = Tier1Measurements._not_measured(
                "chirality is not measurable for every eligible residue "
                f"({chiral_unmeasured} lacked valid N-CA-C-CB geometry)"
            )
        else:
            results["LAW-145"] = {
                "observed": chiral_violations,
                "status": "PASS" if chiral_violations == 0 else "VETO",
                "sample": chiral_checked,
            }
        results["LAW-145"]["granularity"] = GRANULARITY_RESIDUE
        results["LAW-145"]["failing_residues"] = sorted(
            chiral_outlier_residues, key=lambda x: (x.split(":")[0], int(x.split(":")[1]))
        )
        # ═══════════════════════════════════════════════════════════
        # LAW-135: Omega Planarity (Module 1.2.2 — uses coord_math.dihedral_deg)
        # Measures omega dihedral CA(i)-C(i)-N(i+1)-CA(i+1) for sequential pairs.
        # Trans peptide: omega ≈ ±180°. Cis peptide: omega ≈ 0° (typically Pro).
        # Outlier: deviation > 30° from nearest planar ideal (0° or 180°).
        # Reports percentage of outlier omegas.
        # Threshold: 3.0% from station_sop.py
        # Reference: Ramachandran & Sasisekharan (1968)
        # ═══════════════════════════════════════════════════════════
        omega_outliers, omega_total = 0, 0
        omega_outlier_residues = []  # Phase A.5: inline residue capture
        for i in range(len(core) - 1):
            r1, r2 = core[i], core[i + 1]
            if not Tier1Measurements._has_peptide_bond(r1, r2):
                continue
            if not all(a in r1["_atoms"] for a in ["CA", "C"]):
                continue
            if not all(a in r2["_atoms"] for a in ["N", "CA"]):
                continue

            omega, omega_valid = dihedral_deg_valid(
                r1["_atoms"]["CA"], r1["_atoms"]["C"],
                r2["_atoms"]["N"], r2["_atoms"]["CA"]
            )
            if not omega_valid:
                continue

            omega_total += 1
            # Deviation from nearest planar ideal (0° cis or ±180° trans)
            dev_from_trans = 180.0 - abs(omega)
            dev_from_cis = abs(omega)
            min_dev = min(dev_from_trans, dev_from_cis)
            if min_dev > 30.0:
                omega_outliers += 1
                # Capture the peptide bond as r1-r2 pair (report r2 as the affected residue)
                omega_outlier_residues.append(f"{r2['_chain']}:{r2['_seq']}")

        if omega_total == 0:
            results["LAW-135"] = Tier1Measurements._not_measured(
                "no measurable sequential peptide bonds"
            )
            results["LAW-135"]["failing_residues"] = []
        else:
            omega_pct = omega_outliers / omega_total * 100
            results["LAW-135"] = {
                "observed": round(omega_pct, 2),
                "status": "PASS" if omega_pct <= LAW_CANON["LAW-135"]["threshold"] else "VETO",
                "sample": omega_total,
                "failing_residues": sorted(
                    omega_outlier_residues, key=lambda x: (x.split(":")[0], int(x.split(":")[1]))
                )
            }
        results["LAW-135"]["granularity"] = GRANULARITY_RESIDUE
        # ═══════════════════════════════════════════════════════════
        # LAW-150: Rotamer Audit (Module 1.2.5 — uses coord_math.dihedral_deg)
        # Measures chi1 dihedral for non-GLY/ALA core residues.
        # Expected rotamer wells (±30° tolerance):
        #   gauche+ (g+):  60° ± 30°  → [30°, 90°]
        #   gauche- (g-): -60° ± 30°  → [-90°, -30°]
        #   trans   (t):  180° ± 30°  → [150°, 180°] or [-180°, -150°]
        # Outlier: chi1 outside all three wells.
        # Reports percentage of outlier rotamers.
        # Threshold: 20.0% from station_sop.py
        # ═══════════════════════════════════════════════════════════
        rot_outliers, rot_total = 0, 0
        chi1_candidates = 0
        chi1_unmeasured = 0
        rotamer_outliers = []  # Phase A.5: inline residue capture
        for r in core:
            res_name = r["_name"]
            if res_name in NO_CHI1_RESIDUES:
                continue
            if res_name not in CHI1_ATOMS:
                continue
            chi1_candidates += 1
            chi1_def = CHI1_ATOMS[res_name]
            atoms = r["_atoms"]
            if not all(a in atoms for a in chi1_def):
                chi1_unmeasured += 1
                continue

            chi1, chi1_valid = dihedral_deg_valid(
                atoms[chi1_def[0]], atoms[chi1_def[1]],
                atoms[chi1_def[2]], atoms[chi1_def[3]]
            )
            if not chi1_valid:
                chi1_unmeasured += 1
                continue

            rot_total += 1
            # Check if chi1 falls in any expected rotamer well
            abs_chi1 = abs(chi1)
            in_gauche_plus = (30.0 <= chi1 <= 90.0)
            in_gauche_minus = (-90.0 <= chi1 <= -30.0)
            in_trans = (abs_chi1 >= 150.0)  # covers both 150-180 and -180 to -150
            if not (in_gauche_plus or in_gauche_minus or in_trans):
                rot_outliers += 1
                rotamer_outliers.append(f"{r['_chain']}:{r['_seq']}")  # inline capture

        if chi1_candidates == 0:
            results["LAW-150"] = Tier1Measurements._not_applicable(
                "no eligible residue types with a chi1 rotamer"
            )
        elif rot_total == 0 or chi1_unmeasured > 0:
            results["LAW-150"] = Tier1Measurements._not_measured(
                "chi1 is not measurable for every eligible residue "
                f"({chi1_unmeasured} lacked valid chi1 geometry)"
            )
        else:
            rot_pct = rot_outliers / rot_total * 100
            results["LAW-150"] = {
                "observed": round(rot_pct, 2),
                "status": "PASS" if rot_pct <= LAW_CANON["LAW-150"]["threshold"] else "VETO",
                "sample": rot_total,
            }
        results["LAW-150"]["granularity"] = GRANULARITY_RESIDUE
        results["LAW-150"]["failing_residues"] = sorted(
            rotamer_outliers, key=lambda x: (x.split(":")[0], int(x.split(":")[1]))
        )

        # ═══════════════════════════════════════════════════════════
        # LAW-160: Chain Integrity (uses coord_math.distance/get_atom_pos)
        # Measures maximum sequential Cα-Cα distance across all chains.
        # Threshold: 4.5Å from station_sop.py
        # ═══════════════════════════════════════════════════════════
        max_gap = 0.0
        ca_pairs_checked = 0
        missing_ca_pairs = 0
        by_chain = {}
        for residue in residues:
            by_chain.setdefault(residue["_chain"], []).append(residue)
        for chain_id in sorted(by_chain):
            chain_residues = by_chain[chain_id]
            for i in range(len(chain_residues) - 1):
                r1, r2 = chain_residues[i], chain_residues[i + 1]
                if not Tier1Measurements._is_sequential(r1, r2):
                    continue
                if "CA" not in r1["_atoms"] or "CA" not in r2["_atoms"]:
                    missing_ca_pairs += 1
                    continue
                d = distance(r1["_atoms"]["CA"], r2["_atoms"]["CA"])
                ca_pairs_checked += 1
                if d > max_gap:
                    max_gap = d

        if ca_pairs_checked == 0:
            reason = "no eligible sequential CA pairs with both coordinates present"
            if missing_ca_pairs:
                reason = (
                    f"{missing_ca_pairs} sequential pair(s) lack CA coordinates; "
                    "no measurable CA-CA distance"
                )
            results["LAW-160"] = Tier1Measurements._not_measured(reason)
        else:
            results["LAW-160"] = {
                "observed": round(max_gap, 2),
                "status": "PASS" if max_gap <= LAW_CANON["LAW-160"]["threshold"] else "VETO",
                "sample": ca_pairs_checked
            }

        # ═══════════════════════════════════════════════════════════
        # LAW-195: Disulfide Geometry (Module 1.2.4 — uses coord_math.distance)
        # Finds CYS pairs with SG-SG distance < 3.0Å (putative disulfides).
        # Measures deviation from ideal S-S bond length (2.033Å).
        # Reports maximum deviation across all disulfide pairs.
        # Threshold: 0.20Å from station_sop.py
        # Reference: Engh & Huber (1991), Acta Cryst. A47, 392-400
        # ═══════════════════════════════════════════════════════════
        ss_ideal = IDEAL_TABLE["S-S"]
        cys_residues = [r for r in residues if r["_name"] == "CYS"]
        cys_sg = [r for r in cys_residues if "SG" in r["_atoms"]]

        max_ss_dev = 0.0
        ss_pairs_found = 0
        for i in range(len(cys_sg)):
            for j in range(i + 1, len(cys_sg)):
                sg1 = cys_sg[i]["_atoms"]["SG"]
                sg2 = cys_sg[j]["_atoms"]["SG"]
                d = distance(sg1, sg2)
                if d < 3.0:
                    ss_pairs_found += 1
                    dev = abs(d - ss_ideal)
                    if dev > max_ss_dev:
                        max_ss_dev = dev

        if not cys_residues:
            results["LAW-195"] = Tier1Measurements._not_applicable(
                "no cysteine residues are present"
            )
        elif len(cys_sg) < len(cys_residues):
            results["LAW-195"] = Tier1Measurements._not_measured(
                "one or more cysteine SG coordinates are unavailable"
            )
        elif ss_pairs_found == 0:
            results["LAW-195"] = Tier1Measurements._not_applicable(
                "no disulfide pairs were detected"
            )
        else:
            results["LAW-195"] = {
                "observed": round(max_ss_dev, 3),
                "status": "PASS" if max_ss_dev <= LAW_CANON["LAW-195"]["threshold"] else "VETO",
                "sample": ss_pairs_found
            }
        # ═══════════════════════════════════════════════════════════
        # LAW-182: Hydrophobic Burial (Module 1.2.7 — heuristic)
        # Per-chain evaluation of hydrophobic core compactness.
        # For each chain: compute centroid of hydrophobic CA atoms,
        # count fraction within 15Å radius of centroid.
        # Reports minimum ratio across all chains (worst chain governs).
        # Well-folded proteins: ratio ~0.8-1.0.
        # Unfolded/disordered: ratio << 0.3.
        # Threshold: 0.3 from station_sop.py (operator >=)
        # ═══════════════════════════════════════════════════════════
        burial_radius = 15.0
        min_burial_ratio = 1.0
        total_hydrophobic = 0
        hydrophobic_residues = 0
        chains_missing_hydrophobic_ca = []

        for chain_id in get_chains(structure.atoms):
            # Collect CA positions of hydrophobic residues in this chain
            hydro_cas = []
            chain_hydrophobic = 0
            chain_res = get_sequential_residues(structure.atoms, chain_id)
            for res_seq, icode, res_name in chain_res:
                if res_name in HYDROPHOBIC_RESIDUES:
                    hydrophobic_residues += 1
                    chain_hydrophobic += 1
                    ca_pos = get_atom_pos(structure.atoms, chain_id, res_seq, "CA", icode)
                    if ca_pos is not None:
                        hydro_cas.append(ca_pos)

            if chain_hydrophobic > 0 and len(hydro_cas) == 0:
                chains_missing_hydrophobic_ca.append(chain_id)
                continue
            if len(hydro_cas) == 0:
                continue

            total_hydrophobic += len(hydro_cas)
            hydro_array = np.array(hydro_cas)
            centroid = np.mean(hydro_array, axis=0)

            within = sum(1 for pos in hydro_cas if distance(pos, centroid) <= burial_radius)
            ratio = within / len(hydro_cas)
            if ratio < min_burial_ratio:
                min_burial_ratio = ratio

        if hydrophobic_residues == 0:
            results["LAW-182"] = Tier1Measurements._not_applicable(
                "no hydrophobic residues are present"
            )
        elif chains_missing_hydrophobic_ca:
            results["LAW-182"] = Tier1Measurements._not_measured(
                "hydrophobic residues lack CA coordinates on chain(s): "
                + ", ".join(sorted(chains_missing_hydrophobic_ca))
            )
        elif total_hydrophobic != hydrophobic_residues:
            results["LAW-182"] = Tier1Measurements._not_measured(
                "hydrophobic residue CA coverage is incomplete"
            )
        else:
            results["LAW-182"] = {
                "observed": round(min_burial_ratio, 3),
                "status": "PASS" if min_burial_ratio >= LAW_CANON["LAW-182"]["threshold"] else "VETO",
                "sample": total_hydrophobic
            }
        # ═══════════════════════════════════════════════════════════
        # LAW-155: Voxel Occupancy (Module 1.2.8 — heuristic)
        # Divides bounding box into 2Å voxels. Counts voxels containing
        # at least one atom. Reports occupied voxels per residue.
        # Well-resolved structures: ~5-15 voxels/residue.
        # Sparse/incomplete structures: < 2.0 voxels/residue.
        # Threshold: 2.0 from station_sop.py (operator >=)
        # ═══════════════════════════════════════════════════════════
        voxel_size = 2.0
        n_residues = len(residues)
        if len(all_coords) > 0 and n_residues > 0:
            bbox_min = np.min(all_coords, axis=0)
            # Quantize each atom position into a voxel index
            voxel_indices = set()
            for pos in all_coords:
                ix = int(np.floor((pos[0] - bbox_min[0]) / voxel_size))
                iy = int(np.floor((pos[1] - bbox_min[1]) / voxel_size))
                iz = int(np.floor((pos[2] - bbox_min[2]) / voxel_size))
                voxel_indices.add((ix, iy, iz))
            voxels_per_residue = len(voxel_indices) / n_residues
            results["LAW-155"] = {
                "observed": round(voxels_per_residue, 2),
                "status": "PASS" if voxels_per_residue >= LAW_CANON["LAW-155"]["threshold"] else "VETO",
                "sample": len(all_coords)
            }
        else:
            results["LAW-155"] = {
                "observed": None, "status": "NOT_MEASURED", "sample": 0,
                "reason": "no atomic coordinates available"
            }
        # ═══════════════════════════════════════════════════════════
        # LAW-200: Packing Quality (Module 1.2.6 — heuristic)
        # Computes bounding box volume / number of atoms.
        # Well-packed proteins: ~10-20 Å³/atom.
        # Exploded/sparse structures: >> 300 Å³/atom.
        # Threshold: 300.0 Å³/atom from station_sop.py
        # ═══════════════════════════════════════════════════════════
        if len(all_coords) >= 2:
            bbox_min = np.min(all_coords, axis=0)
            bbox_max = np.max(all_coords, axis=0)
            bbox_dims = bbox_max - bbox_min
            # Clamp dimensions to minimum 1.0Å to avoid zero volume for planar/linear structures
            bbox_dims = np.maximum(bbox_dims, 1.0)
            bbox_volume = float(bbox_dims[0] * bbox_dims[1] * bbox_dims[2])
            vol_per_atom = bbox_volume / len(all_coords)
            results["LAW-200"] = {
                "observed": round(vol_per_atom, 2),
                "status": "PASS" if vol_per_atom <= LAW_CANON["LAW-200"]["threshold"] else "VETO",
                "sample": len(all_coords),
                "plddt_advisory": Tier1Measurements._plddt_advisory(residues, is_experimental),
            }
        else:
            results["LAW-200"] = {
                "observed": None, "status": "NOT_MEASURED", "sample": 0,
                "reason": "fewer than two atomic coordinates available",
                "plddt_advisory": Tier1Measurements._plddt_advisory(residues, is_experimental),
            }

        return results, coverage, is_experimental

    @staticmethod
    def compute_structural_characterization(structure):
        atoms = canonicalize_atoms(getattr(structure, "atoms", []) or [])
        res = Tier1Measurements._extract(structure)
        return {
            "total_atoms": len(atoms),
            "total_residues": len(res),
            "source_type": getattr(structure.confidence, "method", "predicted"),
            "resolution": getattr(structure.confidence, "resolution", None),
            "method": "Standard"
        }

    @staticmethod
    def _plddt_observations(structure):
        profile = getattr(structure, "confidence", None)
        values = list(getattr(profile, "_values", []) or [])
        return [
            number for number in (finite_float(value, None) for value in values)
            if number is not None
        ]

    @staticmethod
    def _b_factor_observations(structure):
        profile = getattr(structure, "confidence", None)
        values = list(getattr(profile, "_bfactors", []) or [])
        if values:
            return [
                number for number in (finite_float(value, None) for value in values)
                if number is not None
            ]
        if getattr(profile, "_values", []):
            return []
        return [
            number for number in (
                finite_float(getattr(atom, "b_iso", None), None)
                for atom in getattr(structure, "atoms", [])
            )
            if number is not None
        ]

    @staticmethod
    def detect_confidence_source(structure):
        profile = getattr(structure, "confidence", None)
        method = getattr(profile, "method", "unavailable")
        plddt_values = Tier1Measurements._plddt_observations(structure)
        source = "pLDDT" if plddt_values else getattr(profile, "source_type", "unavailable")
        mean = round(float(np.mean(plddt_values)), 2) if plddt_values else None
        return (source, mean, method)

    @staticmethod
    def decompose_confidence(structure):
        plddt_values = Tier1Measurements._plddt_observations(structure)
        b_values = Tier1Measurements._b_factor_observations(structure)
        mean = round(float(np.mean(plddt_values)), 2) if plddt_values else None
        b_mean = round(float(np.mean(b_values)), 2) if b_values else None
        meta = {
            "schema_version": "confidence-meta/1.0.0",
            "mean": mean,
            "data_available": bool(plddt_values),
            "observation_count": len(plddt_values),
        }
        if b_mean is not None:
            meta["b_factor_mean"] = b_mean
            meta["b_factor_unit"] = "B-factor"
            meta["b_factor_observation_count"] = len(b_values)
        return meta
