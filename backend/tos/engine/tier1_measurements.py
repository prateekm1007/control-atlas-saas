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

WP-SCI-1 measurement-status semantics (directive §2, Constitution §4):
  Every law distinguishes measured PASS / measured VETO / measured FAIL /
  NOT_MEASURED (empty sample) / NOT_APPLICABLE (law does not apply) /
  MEASUREMENT_FAILED (acquisition or tooling failure).
  An empty sample is never a zero-valued measurement and never a
  favorable PASS. Missing data is never converted into a fabricated
  physical violation.

WP-SCI-1 measurement-definition repairs:
  - LAW-120: peptide O(i)-C(i)-N(i+1) now uses the NEXT residue's
    nitrogen (was the same residue's N, which is not bonded to C(i)).
  - LAW-160: CA-CA measured only for immediately sequential residue
    pairs (coord_math.residues_are_adjacent); numbering gaps are not
    physical chain breaks.
  - LAW-195: no applicable disulfide pair -> NOT_APPLICABLE.
"""
from .coord_math import (
    distance, angle_deg, dihedral_deg,
    get_atom_pos, get_chains, get_sequential_residues,
    residues_are_adjacent, peptide_bonded,
    CHI1_ATOMS, NO_CHI1_RESIDUES
)
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


class Tier1Measurements:
    @staticmethod
    def _is_sequential(r1, r2):
        """Immediate sequence adjacency within one chain.

        WP-SCI-1: delegates to the shared predicate in coord_math so that
        every inter-residue law uses one adjacency definition (insertion
        codes honored, chains never joined, numbering gaps never treated
        as physical continuity). Decision semantics are unchanged from
        the historical inline implementation."""
        return residues_are_adjacent(
            r1["_chain"], r1["_seq"], r1.get("_icode", ""),
            r2["_chain"], r2["_seq"], r2.get("_icode", ""),
        )

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
        chain_first_seen = {}
        for a in structure.atoms:
            k = (a.chain_id, a.res_seq, a.insertion_code)
            if k not in res_map:
                res_map[k] = {
                    "_name": a.res_name, "_chain": a.chain_id,
                    "_seq": a.res_seq, "_icode": a.insertion_code,
                    "_conf_acc": [], "_atoms": {}
                }
                chain_first_seen.setdefault(a.chain_id, len(chain_first_seen))
            res_map[k]["_atoms"][a.atom_name] = a.pos
            res_map[k]["_conf_acc"].append(a.b_iso)
        for r in res_map.values():
            r["_conf"] = float(np.mean(r["_conf_acc"]))
        # P0 item 4: residues are returned in FILE ORDER (order of first
        # appearance), grouped by chain in first-appearance order — NOT
        # sorted by (chain, res_seq, icode). Deposits exist that list an
        # insertion-coded residue before its base number (1PPE chain E:
        # 183, 184A, 184, 185); any (res_seq, icode) sort re-pairs the
        # peptide bonds and fabricates backbone gaps. mmCIF structures
        # whose atoms carry label_seq_id are ordered by it instead (see
        # coord_math.get_sequential_residues).
        return sorted(res_map.values(),
                      key=lambda r: chain_first_seen[r["_chain"]])

    @staticmethod
    def run_full_audit(structure, user_intent="NONE"):
        residues = Tier1Measurements._extract(structure)
        if not residues:
            return {}, 0, False

        method = getattr(structure.confidence, "method", "predicted")
        is_experimental = (method != "predicted")
        core = residues if is_experimental else [r for r in residues if r["_conf"] >= 70]
        coverage = 100.0 if is_experimental else (len(core) / len(residues) * 100)

        results = {}

        # ═══════════════════════════════════════════════════════════
        # LAW-100: Bond Integrity (uses coord_math.distance)
        # Measures z-scores of covalent bond lengths against Engh-Huber ideals.
        # Reports % outliers (z > 4σ) for predicted, RMSZ for experimental.
        # ═══════════════════════════════════════════════════════════
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

        # WP-SCI-1: empty sample -> NOT_MEASURED, not a favorable 0.0.
        if z_scores:
            rmsz = float(np.sqrt(np.mean(np.array(z_scores)**2)))
            bond_pct = (outliers / len(z_scores) * 100)
            results["LAW-100"] = {
                "observed": round(rmsz if is_experimental else bond_pct, 2),
                "status": "PASS" if (rmsz <= 8.0 if is_experimental else bond_pct <= 5.0) else "VETO",
                "sample": len(z_scores)
            }
        else:
            results["LAW-100"] = {
                "observed": None, "status": "NOT_MEASURED", "sample": 0,
                "reason": "no covalent bond pairs available in core residues"
            }

        # ═══════════════════════════════════════════════════════════
        # LAW-105: Reliability Coverage
        # ═══════════════════════════════════════════════════════════
        results["LAW-105"] = {
            "observed": round(coverage, 2),
            "status": "PASS" if coverage >= 70.0 else "FAIL",
            "sample": len(residues)
        }

        # ═══════════════════════════════════════════════════════════
        # LAW-110: Backbone Gap (uses coord_math.distance)
        # Measures C(i)-N(i+1) peptide bond distance for sequential pairs.
        # ═══════════════════════════════════════════════════════════
        gaps = 0
        peptide_pairs_checked = 0
        for i in range(len(residues) - 1):
            r1, r2 = residues[i], residues[i + 1]
            if r1["_chain"] == r2["_chain"] and Tier1Measurements._is_sequential(r1, r2):
                if "C" in r1["_atoms"] and "N" in r2["_atoms"]:
                    peptide_pairs_checked += 1
                    if distance(r1["_atoms"]["C"], r2["_atoms"]["N"]) > 2.0:
                        gaps += 1
        # WP-SCI-1: sample now reports the actual number of peptide pairs
        # measured (not len(residues)-1); zero measurable pairs is
        # NOT_MEASURED, not a vacuous PASS.
        if peptide_pairs_checked > 0:
            results["LAW-110"] = {
                "observed": gaps,
                "status": "PASS" if gaps == 0 else "VETO",
                "sample": peptide_pairs_checked
            }
        else:
            results["LAW-110"] = {
                "observed": None, "status": "NOT_MEASURED", "sample": 0,
                "reason": "no sequentially adjacent C-N pairs available"
            }

        # ═══════════════════════════════════════════════════════════
        # LAW-120: Bond Angle RMSD (Module 1.2.1 — uses coord_math.angle_deg)
        # WP-SCI-1 measurement-definition repair (directive §1–§2):
        #   - N(i)-CA(i)-C(i): intra-residue, unchanged.
        #   - O(i)-C(i)-N(i+1): peptide carbonyl angle at C(i), using the
        #     NEXT residue's nitrogen over valid sequential pairs. The
        #     previous code used N(i) — the same residue's amide nitrogen,
        #     which is NOT bonded to C(i) — so the recorded RMSD compared a
        #     non-bonded pseudo-angle against the Engh-Huber peptide ideal.
        #   - CA(i)-C(i)-N(i+1): distinct peptide angle, unchanged.
        # Ideals (IDEAL_TABLE, unchanged): N-CA-C 111.0°, O-C-N 123.0°,
        # CA-C-N 116.6°. Units: degrees. Threshold: LAW_CANON (18.0°).
        # Empty sample -> NOT_MEASURED (never a favorable zero).
        # ═══════════════════════════════════════════════════════════
        angle_sq_diffs = []
        for r in core:
            atoms = r["_atoms"]

            # Intra-residue: N-CA-C
            if all(a in atoms for a in ["N", "CA", "C"]):
                measured = angle_deg(atoms["N"], atoms["CA"], atoms["C"])
                ideal = IDEAL_TABLE["N-CA-C"]
                angle_sq_diffs.append((measured - ideal) ** 2)

        # Peptide angles across valid sequential pairs only: same chain,
        # no missing span (file order), AND chemically bonded — C(i)-N(i+1)
        # within the peptide-bond tolerance (P0 item 4). Angles are never
        # measured across a discontinuity: a pseudo-angle against the
        # Engh-Huber peptide ideal is exactly the defect class repaired in
        # WP-SCI-1.
        for i in range(len(core) - 1):
            r1, r2 = core[i], core[i + 1]
            if not Tier1Measurements._is_sequential(r1, r2):
                continue
            if not peptide_bonded(r1["_atoms"], r2["_atoms"]):
                continue
            if "C" in r1["_atoms"] and "N" in r2["_atoms"]:
                if "O" in r1["_atoms"]:
                    measured = angle_deg(r1["_atoms"]["O"], r1["_atoms"]["C"], r2["_atoms"]["N"])
                    ideal = IDEAL_TABLE["O-C-N"]
                    angle_sq_diffs.append((measured - ideal) ** 2)
                if "CA" in r1["_atoms"]:
                    measured = angle_deg(r1["_atoms"]["CA"], r1["_atoms"]["C"], r2["_atoms"]["N"])
                    ideal = IDEAL_TABLE["CA-C-N"]
                    angle_sq_diffs.append((measured - ideal) ** 2)

        if angle_sq_diffs:
            angle_rmsd = float(np.sqrt(np.mean(angle_sq_diffs)))
            results["LAW-120"] = {
                "observed": round(angle_rmsd, 2),
                "status": "PASS" if angle_rmsd <= LAW_CANON["LAW-120"]["threshold"] else "VETO",
                "sample": len(angle_sq_diffs)
            }
        else:
            # WP-SCI-1: an empty sample is not a zero-valued measurement.
            results["LAW-120"] = {
                "observed": None, "status": "NOT_MEASURED", "sample": 0,
                "reason": "no backbone angle samples (no residues with N/CA/C and no valid sequential pairs)"
            }

        # ═══════════════════════════════════════════════════════════
        # LAW-125: Ramachandran (uses coord_math.dihedral_deg)
        # Computes phi/psi dihedrals and flags outliers.
        # ═══════════════════════════════════════════════════════════
        rama_out, rama_tot = 0, 0
        rama_outliers = []  # Phase A.5: inline residue capture
        for i in range(1, len(core) - 1):
            p, c, n = core[i - 1], core[i], core[i + 1]
            if not (Tier1Measurements._is_sequential(p, c) and
                    Tier1Measurements._is_sequential(c, n)):
                continue
            # P0 item 4: phi/psi only across genuinely bonded pairs.
            if not (peptide_bonded(p["_atoms"], c["_atoms"]) and
                    peptide_bonded(c["_atoms"], n["_atoms"])):
                continue
            if ("C" not in p["_atoms"] or
                    not all(k in c["_atoms"] for k in ["N", "CA", "C"]) or
                    "N" not in n["_atoms"]):
                continue

            phi = dihedral_deg(
                p["_atoms"]["C"], c["_atoms"]["N"],
                c["_atoms"]["CA"], c["_atoms"]["C"]
            )
            psi = dihedral_deg(
                c["_atoms"]["N"], c["_atoms"]["CA"],
                c["_atoms"]["C"], n["_atoms"]["N"]
            )

            # dihedral_deg returns float (never None), but 0.0 on degenerate geometry.
            # Only count if both are non-degenerate (non-zero or genuinely zero).
            # In practice, real backbone dihedrals are never exactly 0.0 for both.
            if phi == 0.0 and psi == 0.0:
                continue
            rama_tot += 1
            if Tier1Measurements._is_rama_outlier(c["_name"], phi, psi):
                rama_out += 1
                rama_outliers.append(f"{c['_chain']}:{c['_seq']}")  # chain:seq format

        # WP-SCI-1: zero measurable phi/psi triples is NOT_MEASURED.
        if rama_tot > 0:
            rama_pct = (rama_out / rama_tot * 100)
            results["LAW-125"] = {
                "observed": round(rama_pct, 2),
                "status": "PASS" if rama_pct <= LAW_CANON["LAW-125"]["threshold"] else "VETO",
                "sample": rama_tot,
                "granularity": GRANULARITY_RESIDUE,
                "failing_residues": sorted(rama_outliers, key=lambda x: (x.split(":")[0], int(x.split(":")[1])))
            }
        else:
            results["LAW-125"] = {
                "observed": None, "status": "NOT_MEASURED", "sample": 0,
                "granularity": GRANULARITY_RESIDUE, "failing_residues": [],
                "reason": "no consecutive residue triples with complete phi/psi atoms"
            }

        # ═══════════════════════════════════════════════════════════
        # LAW-130: Clashscore (Neighbor Exclusion)
        # Uses scipy KDTree for efficiency. Geometry via numpy directly
        # (pairwise distance already handled by KDTree radius query).
        # ═══════════════════════════════════════════════════════════
        clash = 0
        clash_pairs = []  # Phase A.5: inline residue-pair capture
        all_coords = np.array([a.pos for a in structure.atoms])
        clash_search_failed = False  # WP-SCI-1: distinguish tooling failure
        if len(all_coords) > 1:
            try:
                from scipy.spatial import cKDTree
                tree = cKDTree(all_coords)
                pairs = tree.query_pairs(r=1.5)
                for ii, jj in pairs:
                    ai, aj = structure.atoms[ii], structure.atoms[jj]
                    if ai.res_seq == aj.res_seq and ai.chain_id == aj.chain_id:
                        continue
                    if ai.chain_id == aj.chain_id and abs(ai.res_seq - aj.res_seq) <= 1:
                        continue
                    clash += 1
                    # Capture as residue pair, deterministic order (lower chain:seq first)
                    r1 = f"{ai.chain_id}:{ai.res_seq}"
                    r2 = f"{aj.chain_id}:{aj.res_seq}"
                    pair = tuple(sorted([r1, r2], key=lambda x: (x.split(":")[0], int(x.split(":")[1]))))
                    clash_pairs.append(pair)
            except Exception:
                # WP-SCI-1: a failed clash search must not masquerade as a
                # measured, clash-free structure.
                clash_search_failed = True
        # Deduplicate and sort pairs
        unique_pairs = sorted(set(clash_pairs), key=lambda p: (p[0].split(":")[0], int(p[0].split(":")[1]), p[1].split(":")[0], int(p[1].split(":")[1])))
        # WP-SCI-1: explicit non-success states for LAW-130.
        if clash_search_failed:
            results["LAW-130"] = {
                "observed": None, "status": "MEASUREMENT_FAILED", "sample": 0,
                "granularity": GRANULARITY_RESIDUE_PAIR, "failing_residue_pairs": [],
                "reason": "clash search unavailable (neighbor-search tooling failure)"
            }
        elif len(all_coords) == 0:
            results["LAW-130"] = {
                "observed": None, "status": "NOT_MEASURED", "sample": 0,
                "granularity": GRANULARITY_RESIDUE_PAIR, "failing_residue_pairs": [],
                "reason": "no atomic coordinates available"
            }
        else:
            clash_score = (clash / len(all_coords)) * 1000
            results["LAW-130"] = {
                "observed": round(clash_score, 2),
                "status": "PASS" if clash_score < LAW_CANON["LAW-130"]["threshold"] else "VETO",
                "sample": len(all_coords),
                "granularity": GRANULARITY_RESIDUE_PAIR,
                "failing_residue_pairs": unique_pairs
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
        # STUBS: Laws awaiting Phase 1 implementation
        # These return hardcoded PASS values. Each will be replaced
        # in subsequent modules with real coord_math calculations.
        # ═══════════════════════════════════════════════════════════
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
        chiral_outlier_residues = []  # Phase A.5: inline residue capture
        for r in core:
            if r["_name"] == "GLY":
                continue
            atoms = r["_atoms"]
            if not all(a in atoms for a in ["N", "CA", "C", "CB"]):
                continue
            chiral_checked += 1
            # Improper dihedral: N-CA-C-CB
            # L-amino acid → negative, D-amino acid → positive
            improper = dihedral_deg(
                atoms["N"], atoms["CA"], atoms["C"], atoms["CB"]
            )
            # Genuine D-amino acid: improper between 0° and 90°
            # Values near ±180° indicate degenerate/coplanar geometry
            # (CB in the N-CA-C plane) — refinement artifact, not chirality error.
            # Reference: 4HHB Chain D Res 78 LEU has improper=179.7° with
            # CB only 0.009Å from the N-CA-C plane.
            if 0.0 < improper < 90.0:
                chiral_violations += 1
                chiral_outlier_residues.append(f"{r['_chain']}:{r['_seq']}")

        # WP-SCI-1: no chiral center checked (all GLY / missing CB or N/CA/C)
        # is NOT_MEASURED, not a vacuous zero-violation PASS.
        if chiral_checked > 0:
            results["LAW-145"] = {
                "observed": chiral_violations,
                "status": "PASS" if chiral_violations == 0 else "VETO",
                "sample": chiral_checked,
                "granularity": GRANULARITY_RESIDUE,
                "failing_residues": sorted(chiral_outlier_residues, key=lambda x: (x.split(":")[0], int(x.split(":")[1])))
            }
        else:
            results["LAW-145"] = {
                "observed": None, "status": "NOT_MEASURED", "sample": 0,
                "granularity": GRANULARITY_RESIDUE, "failing_residues": [],
                "reason": "no residue with N/CA/C/CB available to check chirality"
            }
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
            if not Tier1Measurements._is_sequential(r1, r2):
                continue
            # P0 item 4: omega only across a genuinely bonded pair.
            if not peptide_bonded(r1["_atoms"], r2["_atoms"]):
                continue
            if not all(a in r1["_atoms"] for a in ["CA", "C"]):
                continue
            if not all(a in r2["_atoms"] for a in ["N", "CA"]):
                continue

            omega = dihedral_deg(
                r1["_atoms"]["CA"], r1["_atoms"]["C"],
                r2["_atoms"]["N"], r2["_atoms"]["CA"]
            )

            omega_total += 1
            # Deviation from nearest planar ideal (0° cis or ±180° trans)
            dev_from_trans = 180.0 - abs(omega)
            dev_from_cis = abs(omega)
            min_dev = min(dev_from_trans, dev_from_cis)
            if min_dev > 30.0:
                omega_outliers += 1
                # Capture the peptide bond as r1-r2 pair (report r2 as the affected residue)
                omega_outlier_residues.append(f"{r2['_chain']}:{r2['_seq']}")

        # WP-SCI-1: zero measurable omega pairs is NOT_MEASURED.
        if omega_total > 0:
            omega_pct = (omega_outliers / omega_total * 100)
            results["LAW-135"] = {
                "observed": round(omega_pct, 2),
                "status": "PASS" if omega_pct <= LAW_CANON["LAW-135"]["threshold"] else "VETO",
                "sample": omega_total,
                "granularity": GRANULARITY_RESIDUE,
                "failing_residues": sorted(omega_outlier_residues, key=lambda x: (x.split(":")[0], int(x.split(":")[1])))
            }
        else:
            results["LAW-135"] = {
                "observed": None, "status": "NOT_MEASURED", "sample": 0,
                "granularity": GRANULARITY_RESIDUE, "failing_residues": [],
                "reason": "no sequentially adjacent CA-C-N-CA omega pairs available"
            }
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
        rotamer_outliers = []  # Phase A.5: inline residue capture
        for r in core:
            res_name = r["_name"]
            if res_name in NO_CHI1_RESIDUES:
                continue
            if res_name not in CHI1_ATOMS:
                continue
            chi1_def = CHI1_ATOMS[res_name]
            atoms = r["_atoms"]
            if not all(a in atoms for a in chi1_def):
                continue

            rot_total += 1
            chi1 = dihedral_deg(
                atoms[chi1_def[0]], atoms[chi1_def[1]],
                atoms[chi1_def[2]], atoms[chi1_def[3]]
            )

            # Check if chi1 falls in any expected rotamer well
            abs_chi1 = abs(chi1)
            in_gauche_plus = (30.0 <= chi1 <= 90.0)
            in_gauche_minus = (-90.0 <= chi1 <= -30.0)
            in_trans = (abs_chi1 >= 150.0)  # covers both 150-180 and -180 to -150
            if not (in_gauche_plus or in_gauche_minus or in_trans):
                rot_outliers += 1
                rotamer_outliers.append(f"{r['_chain']}:{r['_seq']}")  # inline capture

        # WP-SCI-1: zero rotatable chi1 side chains checked is NOT_MEASURED.
        if rot_total > 0:
            rot_pct = (rot_outliers / rot_total * 100)
            results["LAW-150"] = {
                "observed": round(rot_pct, 2),
                "status": "PASS" if rot_pct <= LAW_CANON["LAW-150"]["threshold"] else "VETO",
                "sample": rot_total,
                "granularity": GRANULARITY_RESIDUE,
                "failing_residues": sorted(rotamer_outliers, key=lambda x: (x.split(":")[0], int(x.split(":")[1])))
            }
        else:
            results["LAW-150"] = {
                "observed": None, "status": "NOT_MEASURED", "sample": 0,
                "granularity": GRANULARITY_RESIDUE, "failing_residues": [],
                "reason": "no side chain with a measurable chi1 dihedral available"
            }

        # ═══════════════════════════════════════════════════════════
        # LAW-160: Chain Integrity (uses coord_math.distance/get_atom_pos)
        # WP-SCI-1: CA(i)-CA(i+1) is measured only for immediately
        # sequential residue pairs (coord_math.residues_are_adjacent).
        # A missing residue span (numbering gap, e.g. 45 -> 67) is absent
        # density, not a physical chain break, and is skipped. Chains are
        # never joined; insertion codes are honored; residues with missing
        # CA atoms are skipped. Genuine coordinate discontinuities between
        # adjacent-numbered residues still VETO. Empty sample ->
        # NOT_MEASURED. Threshold: 4.5Å from station_sop.py.
        # ═══════════════════════════════════════════════════════════
        max_gap = None
        ca_pairs_checked = 0
        for chain_id in get_chains(structure.atoms):
            chain_residues = get_sequential_residues(structure.atoms, chain_id)
            for i in range(len(chain_residues) - 1):
                curr_res = chain_residues[i]
                next_res = chain_residues[i + 1]
                if not residues_are_adjacent(chain_id, curr_res[0], curr_res[1],
                                             chain_id, next_res[0], next_res[1]):
                    continue
                p1 = get_atom_pos(structure.atoms, chain_id, curr_res[0], "CA", curr_res[1])
                p2 = get_atom_pos(structure.atoms, chain_id, next_res[0], "CA", next_res[1])
                if p1 is not None and p2 is not None:
                    d = distance(p1, p2)
                    ca_pairs_checked += 1
                    if max_gap is None or d > max_gap:
                        max_gap = d

        if ca_pairs_checked > 0:
            results["LAW-160"] = {
                "observed": round(max_gap, 2),
                "status": "PASS" if max_gap <= LAW_CANON["LAW-160"]["threshold"] else "VETO",
                "sample": ca_pairs_checked
            }
        else:
            results["LAW-160"] = {
                "observed": None, "status": "NOT_MEASURED", "sample": 0,
                "reason": "no sequentially adjacent CA pairs available"
            }

        # ═══════════════════════════════════════════════════════════
        # STUBS: Remaining laws
        # ═══════════════════════════════════════════════════════════
        # ═══════════════════════════════════════════════════════════
        # LAW-195: Disulfide Geometry (Module 1.2.4 — uses coord_math.distance)
        # Finds CYS pairs with SG-SG distance < 3.0Å (putative disulfides).
        # Measures deviation from ideal S-S bond length (2.033Å).
        # Reports maximum deviation across all disulfide pairs.
        # Threshold: 0.20Å from station_sop.py
        # Reference: Engh & Huber (1991), Acta Cryst. A47, 392-400
        # ═══════════════════════════════════════════════════════════
        ss_ideal = IDEAL_TABLE["S-S"]
        cys_sg = []
        for r in residues:
            if r["_name"] == "CYS" and "SG" in r["_atoms"]:
                cys_sg.append(r)

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

        # WP-SCI-1 (directive §2): a structure with no applicable disulfide
        # pair has NOT passed this law — the law simply does not apply.
        # Recording it as PASS claimed a measurement that never happened.
        if ss_pairs_found > 0:
            results["LAW-195"] = {
                "observed": round(max_ss_dev, 3),
                "status": "PASS" if max_ss_dev <= LAW_CANON["LAW-195"]["threshold"] else "VETO",
                "sample": ss_pairs_found
            }
        else:
            results["LAW-195"] = {
                "observed": None, "status": "NOT_APPLICABLE", "sample": 0,
                "reason": "no CYS-SG pair within 3.0 Å — no applicable disulfide to measure"
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
        # WP-SCI-1: distinguish "no hydrophobic residues at all" (law does
        # not apply) from "hydrophobic residues present but not measurable"
        # (missing CA atoms).
        hydro_named = sum(1 for r in residues if r["_name"] in HYDROPHOBIC_RESIDUES)

        for chain_id in get_chains(structure.atoms):
            # Collect CA positions of hydrophobic residues in this chain
            hydro_cas = []
            chain_res = get_sequential_residues(structure.atoms, chain_id)
            for res_seq, icode, res_name in chain_res:
                if res_name in HYDROPHOBIC_RESIDUES:
                    ca_pos = get_atom_pos(structure.atoms, chain_id, res_seq, "CA", icode)
                    if ca_pos is not None:
                        hydro_cas.append(ca_pos)

            if len(hydro_cas) == 0:
                continue

            total_hydrophobic += len(hydro_cas)
            hydro_array = np.array(hydro_cas)
            centroid = np.mean(hydro_array, axis=0)

            within = sum(1 for pos in hydro_cas if distance(pos, centroid) <= burial_radius)
            ratio = within / len(hydro_cas)
            if ratio < min_burial_ratio:
                min_burial_ratio = ratio

        # WP-SCI-1: the historical "vacuously true" 1.0 fallback is gone.
        if total_hydrophobic > 0:
            results["LAW-182"] = {
                "observed": round(min_burial_ratio, 3),
                "status": "PASS" if min_burial_ratio >= LAW_CANON["LAW-182"]["threshold"] else "VETO",
                "sample": total_hydrophobic
            }
        elif hydro_named > 0:
            results["LAW-182"] = {
                "observed": None, "status": "NOT_MEASURED", "sample": 0,
                "reason": "hydrophobic residues present but no CA positions measurable"
            }
        else:
            results["LAW-182"] = {
                "observed": None, "status": "NOT_APPLICABLE", "sample": 0,
                "reason": "no hydrophobic residues — no hydrophobic core to measure"
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
            # WP-SCI-1: missing coordinates must not become a fabricated
            # VETO (the historical 0.0 fallback read as a physical violation).
            results["LAW-155"] = {
                "observed": None, "status": "NOT_MEASURED", "sample": 0,
                "reason": "no atomic coordinates available to grid"
            }
        # ═══════════════════════════════════════════════════════════
        # LAW-200: Packing Quality (Module 1.2.6 — heuristic)
        # Computes bounding box volume / number of atoms.
        # Well-packed proteins: ~10-20 Å³/atom.
        # Exploded/sparse structures: >> 300 Å³/atom.
        # Threshold: 300.0 Å³/atom from station_sop.py
        # ═══════════════════════════════════════════════════════════
        if len(all_coords) > 0:
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
                "sample": len(all_coords)
            }
        else:
            # WP-SCI-1: missing coordinates must not become a favorable 0.0.
            results["LAW-200"] = {
                "observed": None, "status": "NOT_MEASURED", "sample": 0,
                "reason": "no atomic coordinates available"
            }

        return results, coverage, False

    @staticmethod
    def compute_structural_characterization(structure):
        res = Tier1Measurements._extract(structure)
        return {
            "total_atoms": len(structure.atoms),
            "total_residues": len(res),
            "source_type": getattr(structure.confidence, "method", "predicted"),
            "resolution": getattr(structure.confidence, "resolution", None),
            "method": "Standard"
        }

    @staticmethod
    def detect_confidence_source(structure):
        return ("pLDDT", 90.0, "static")

    @staticmethod
    def decompose_confidence(structure):
        return {"mean": 90.0, "data_available": True}
