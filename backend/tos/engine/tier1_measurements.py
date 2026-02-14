"""
TOSCANINI Tier-1 Physical Invariant Engine v22.4.3
All 15 laws implemented. Residue-specific physics, insertion-aware adjacency, stratified auditing.
"""
import numpy as np
import logging
from ..governance.station_sop import (
    LAW_CANON, SIGMA_TABLE, IDEAL_TABLE,
    STANDARD_RESIDUES, HYDROPHOBIC_RESIDUES, LAW_METHOD_CLASSIFICATIONS
)

logger = logging.getLogger("toscanini.tier1")


# ── Geometry Primitives ─────────────────────────────────────
def _dihedral(p1, p2, p3, p4):
    """Standard right-handed dihedral (BioPython/IUPAC)."""
    b0 = -1.0 * (p2 - p1)
    b1 = p3 - p2
    b2 = p4 - p3
    b1_n = np.linalg.norm(b1)
    if b1_n < 1e-8:
        return None
    b1 /= b1_n
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    vn, wn = np.linalg.norm(v), np.linalg.norm(w)
    if vn < 1e-8 or wn < 1e-8:
        return None
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def _angle(p1, p2, p3):
    """3-point angle in degrees."""
    v1, v2 = p1 - p2, p3 - p2
    norm = np.linalg.norm(v1) * np.linalg.norm(v2)
    if norm < 1e-8:
        return None
    return float(np.degrees(np.arccos(np.clip(np.dot(v1, v2) / norm, -1.0, 1.0))))


def _dist(p1, p2):
    """Euclidean distance."""
    return float(np.linalg.norm(p1 - p2))


# ── Core Engine ─────────────────────────────────────────────
class Tier1Measurements:

    @staticmethod
    def _is_sequential(r1, r2):
        """Insertion-aware adjacency: 45 -> 45A -> 46."""
        if r1["_chain"] != r2["_chain"]:
            return False
        s1, s2 = r1["_seq"], r2["_seq"]
        i1, i2 = r1.get("_icode", ""), r2.get("_icode", "")
        if s2 - s1 == 1 and not i1 and not i2:
            return True
        if s2 == s1 and i2 > i1:
            return True
        if s2 - s1 == 1 and i1 and not i2:
            return True
        return False

    @staticmethod
    def _is_rama_outlier(res_name, phi, psi):
        """Residue-specific Ramachandran (rectangular Lovell approximations)."""
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
        """Extract residue map with atoms, confidence, chain info."""
        res_map = {}
        for a in structure.atoms:
            k = (a.chain_id, a.res_seq, a.insertion_code)
            if k not in res_map:
                res_map[k] = {
                    "_name": a.res_name, "_chain": a.chain_id,
                    "_seq": a.res_seq, "_icode": a.insertion_code,
                    "_conf_acc": [], "_atoms": {},
                }
            res_map[k]["_atoms"][a.atom_name] = a.pos
            res_map[k]["_conf_acc"].append(a.b_iso)
        for r in res_map.values():
            r["_conf"] = float(np.mean(r["_conf_acc"]))
        return sorted(res_map.values(), key=lambda x: (x["_chain"], x["_seq"], x["_icode"]))

    # ── Main Audit ──────────────────────────────────────────
    @staticmethod
    def run_full_audit(structure, user_intent="NONE"):
        try:
            residues = Tier1Measurements._extract(structure)
        except Exception as e:
            return {lid: ("FAIL", f"Extraction error: {e}", "ERR", "error") for lid in LAW_CANON}, 0, False

        if not residues:
            return {lid: ("FAIL", "No residues", "ERR", "error") for lid in LAW_CANON}, 0, False

        # Epistemic stratification
        core = [r for r in residues if r["_conf"] >= 70]
        fringe = [r for r in residues if r["_conf"] < 70]
        coverage = (len(core) / len(residues) * 100) if residues else 0

        results = {}
        fatal_fringe = False

        # Collect all atom coordinates for spatial laws
        all_coords = []
        all_elements = []
        for a in structure.atoms:
            all_coords.append(a.pos)
            all_elements.append(a.element)
        all_coords = np.array(all_coords) if all_coords else np.empty((0, 3))

        # ════════════════════════════════════════════════════
        # LAW-100: Bond Integrity (4-sigma Engh-Huber)
        # ════════════════════════════════════════════════════
        def _count_bond_fails(r_list):
            cnt = 0
            for r in r_list:
                for b, ideal in IDEAL_TABLE.items():
                    parts = b.split('-')
                    if len(parts) != 2:
                        continue
                    a1, a2 = parts
                    if a1 in r["_atoms"] and a2 in r["_atoms"]:
                        d = _dist(r["_atoms"][a1], r["_atoms"][a2])
                        sigma = SIGMA_TABLE.get(b, 0.02)
                        z = abs(d - ideal) / sigma if sigma > 0 else 0
                        if z > 4.0:
                            cnt += 1
                            break
            return cnt

        c_bond = _count_bond_fails(core)
        f_bond = _count_bond_fails(fringe)
        if f_bond > 0:
            fatal_fringe = True
        results["LAW-100"] = (
            "PASS" if c_bond == 0 else "VETO",
            f"Core outliers: {c_bond}, Fringe: {f_bond}",
            "sigma:4.0", "deterministic",
        )

        # ════════════════════════════════════════════════════
        # LAW-105: Reliability Coverage
        # ════════════════════════════════════════════════════
        cov_thresh = 30.0
        results["LAW-105"] = (
            "PASS" if coverage >= cov_thresh else "FAIL",
            f"Coverage: {round(coverage, 1)}% (threshold: {cov_thresh}%)",
            f"THRESH:{cov_thresh}%", "deterministic",
        )

        # ════════════════════════════════════════════════════
        # LAW-110: Backbone Gap Detection
        # ════════════════════════════════════════════════════
        gaps = 0
        for i in range(len(residues) - 1):
            r1, r2 = residues[i], residues[i + 1]
            if not Tier1Measurements._is_sequential(r1, r2):
                continue
            if "C" in r1["_atoms"] and "N" in r2["_atoms"]:
                d = _dist(r1["_atoms"]["C"], r2["_atoms"]["N"])
                if d > 1.5:
                    gaps += 1
        results["LAW-110"] = (
            "PASS" if gaps == 0 else "VETO",
            f"{gaps} gaps > 1.5A",
            "THRESH:1.5A", "deterministic",
        )

        # ════════════════════════════════════════════════════
        # LAW-120: Bond Angle Deviation
        # ════════════════════════════════════════════════════
        angle_deviations = []
        for r in core:
            for b, ideal in IDEAL_TABLE.items():
                parts = b.split('-')
                if len(parts) != 3:
                    continue
                a1, a2, a3 = parts
                if a1 in r["_atoms"] and a2 in r["_atoms"] and a3 in r["_atoms"]:
                    ang = _angle(r["_atoms"][a1], r["_atoms"][a2], r["_atoms"][a3])
                    if ang is not None:
                        angle_deviations.append(abs(ang - ideal))
        mean_angle_dev = float(np.mean(angle_deviations)) if angle_deviations else 0.0
        results["LAW-120"] = (
            "PASS" if mean_angle_dev < 4.6 else "VETO",
            f"Mean deviation: {round(mean_angle_dev, 2)} deg",
            "THRESH:4.6deg", "deterministic",
        )

        # ════════════════════════════════════════════════════
        # LAW-125: Ramachandran Outliers
        # ════════════════════════════════════════════════════
        def _count_rama(r_list):
            out, tot = 0, 0
            for i in range(1, len(r_list) - 1):
                p, c, n = r_list[i - 1], r_list[i], r_list[i + 1]
                if not (Tier1Measurements._is_sequential(p, c) and Tier1Measurements._is_sequential(c, n)):
                    continue
                if "C" not in p["_atoms"] or "N" not in c["_atoms"] or "CA" not in c["_atoms"] or "C" not in c["_atoms"] or "N" not in n["_atoms"]:
                    continue
                phi = _dihedral(p["_atoms"]["C"], c["_atoms"]["N"], c["_atoms"]["CA"], c["_atoms"]["C"])
                psi = _dihedral(c["_atoms"]["N"], c["_atoms"]["CA"], c["_atoms"]["C"], n["_atoms"]["N"])
                if phi is None or psi is None:
                    continue
                tot += 1
                if Tier1Measurements._is_rama_outlier(c["_name"], phi, psi):
                    out += 1
            return out, tot

        c_out, c_tot = _count_rama(core)
        f_out, f_tot = _count_rama(fringe)
        rama_pct = (c_out / c_tot * 100) if c_tot > 0 else 0
        if f_tot > 0 and f_out > (f_tot * 0.2):
            fatal_fringe = True
        results["LAW-125"] = (
            "PASS" if rama_pct <= 5.0 else "VETO",
            f"{round(rama_pct, 1)}% core outliers ({c_out}/{c_tot})",
            "THRESH:5%", "deterministic",
        )

        # ════════════════════════════════════════════════════
        # LAW-130: Steric Clash Detection (KD-tree)
        # ════════════════════════════════════════════════════
        clash_count = 0
        if len(all_coords) > 1:
            try:
                from scipy.spatial import cKDTree
                tree = cKDTree(all_coords)
                pairs = tree.query_pairs(r=1.5)
                # Filter: exclude bonded atoms (same residue or sequential)
                for i, j in pairs:
                    ai, aj = structure.atoms[i], structure.atoms[j]
                    # Same residue — likely bonded
                    if ai.res_seq == aj.res_seq and ai.chain_id == aj.chain_id:
                        continue
                    # Sequential residues, backbone bond (C-N)
                    if abs(ai.res_seq - aj.res_seq) == 1 and ai.chain_id == aj.chain_id:
                        if {ai.atom_name, aj.atom_name} == {"C", "N"}:
                            continue
                    d = _dist(ai.pos, aj.pos)
                    if d < 1.5 and d > 0.01:
                        clash_count += 1
            except ImportError:
                logger.warning("scipy not available for clash detection")
        clash_pct = (clash_count / max(len(all_coords), 1)) * 100
        results["LAW-130"] = (
            "PASS" if clash_pct < 0.4 else "VETO",
            f"{clash_count} clashes ({round(clash_pct, 2)}%)",
            "THRESH:0.4%", "deterministic",
        )

        # ════════════════════════════════════════════════════
        # LAW-135: Omega (Peptide Bond) Planarity
        # ════════════════════════════════════════════════════
        omega_outliers, omega_total = 0, 0
        for i in range(len(core) - 1):
            r1, r2 = core[i], core[i + 1]
            if not Tier1Measurements._is_sequential(r1, r2):
                continue
            if "CA" in r1["_atoms"] and "C" in r1["_atoms"] and "N" in r2["_atoms"] and "CA" in r2["_atoms"]:
                omega = _dihedral(r1["_atoms"]["CA"], r1["_atoms"]["C"], r2["_atoms"]["N"], r2["_atoms"]["CA"])
                if omega is not None:
                    omega_total += 1
                    # Omega should be ~180 (trans) or ~0 (cis, rare)
                    deviation = min(abs(abs(omega) - 180), abs(omega))
                    if deviation > 30:
                        omega_outliers += 1
        omega_pct = (omega_outliers / omega_total * 100) if omega_total > 0 else 0
        results["LAW-135"] = (
            "PASS" if omega_pct < 3.0 else "VETO",
            f"{round(omega_pct, 1)}% non-planar ({omega_outliers}/{omega_total})",
            "THRESH:3%", "deterministic",
        )

        # ════════════════════════════════════════════════════
        # LAW-145: Chirality (D-amino acid detection)
        # ════════════════════════════════════════════════════
        d_amino = 0
        for r in core:
            if r["_name"] == "GLY":
                continue
            if "N" in r["_atoms"] and "CA" in r["_atoms"] and "C" in r["_atoms"] and "CB" in r["_atoms"]:
                # Improper dihedral N-CA-C-CB: positive for L, negative for D
                imp = _dihedral(r["_atoms"]["N"], r["_atoms"]["CA"], r["_atoms"]["C"], r["_atoms"]["CB"])
                if imp is not None and imp < 0:
                    d_amino += 1
        if d_amino > 0 and fringe:
            fatal_fringe = True
        results["LAW-145"] = (
            "PASS" if d_amino == 0 else "VETO",
            f"{d_amino} D-amino acids detected",
            "THRESH:0", "deterministic",
        )

        # ════════════════════════════════════════════════════
        # LAW-150: Rotamer Audit (Chi1 outliers)
        # ════════════════════════════════════════════════════
        chi1_outliers, chi1_total = 0, 0
        for r in core:
            if r["_name"] in ("GLY", "ALA"):
                continue
            if "N" in r["_atoms"] and "CA" in r["_atoms"] and "CB" in r["_atoms"]:
                # Need CG or OG or SG
                cg = None
                for g_name in ["CG", "CG1", "OG", "OG1", "SG"]:
                    if g_name in r["_atoms"]:
                        cg = r["_atoms"][g_name]
                        break
                if cg is not None:
                    chi1 = _dihedral(r["_atoms"]["N"], r["_atoms"]["CA"], r["_atoms"]["CB"], cg)
                    if chi1 is not None:
                        chi1_total += 1
                        # Common rotamers: gauche+ (~60), gauche- (~-60), trans (~180)
                        closest = min(abs(chi1 - 60), abs(chi1 + 60), abs(abs(chi1) - 180))
                        if closest > 40:
                            chi1_outliers += 1
        chi1_pct = (chi1_outliers / chi1_total * 100) if chi1_total > 0 else 0
        results["LAW-150"] = (
            "PASS" if chi1_pct < 20.0 else "VETO",
            f"{round(chi1_pct, 1)}% outliers ({chi1_outliers}/{chi1_total})",
            "THRESH:20%", "deterministic",
        )

        # ════════════════════════════════════════════════════
        # LAW-155: Voxel Occupancy (connectivity proxy) [HEURISTIC]
        # ════════════════════════════════════════════════════
        if len(all_coords) > 10:
            try:
                from scipy.spatial import cKDTree
                tree = cKDTree(all_coords)
                neighbors = tree.query_ball_point(all_coords, r=3.0)
                mean_neighbors = float(np.mean([len(n) - 1 for n in neighbors]))
                voxel_ok = mean_neighbors >= 2.0
            except ImportError:
                mean_neighbors = 0.0
                voxel_ok = True
        else:
            mean_neighbors = 0.0
            voxel_ok = True
        results["LAW-155"] = (
            "PASS" if voxel_ok else "FAIL",
            f"Mean neighbors: {round(mean_neighbors, 1)}",
            "THRESH:2.0", "heuristic",
        )

        # ════════════════════════════════════════════════════
        # LAW-160: Chain Integrity (CA trace continuity)
        # ════════════════════════════════════════════════════
        chain_breaks = 0
        ca_residues = [r for r in residues if "CA" in r["_atoms"]]
        for i in range(len(ca_residues) - 1):
            r1, r2 = ca_residues[i], ca_residues[i + 1]
            if r1["_chain"] != r2["_chain"]:
                continue
            if not Tier1Measurements._is_sequential(r1, r2):
                continue
            d = _dist(r1["_atoms"]["CA"], r2["_atoms"]["CA"])
            if d > 4.2:
                chain_breaks += 1
        results["LAW-160"] = (
            "PASS" if chain_breaks == 0 else "VETO",
            f"{chain_breaks} breaks > 4.2A",
            "THRESH:4.2A", "deterministic",
        )

        # ════════════════════════════════════════════════════
        # LAW-170: Residue Identity (unknown residues)
        # ════════════════════════════════════════════════════
        unknowns = [r["_name"] for r in residues if r["_name"] not in STANDARD_RESIDUES]
        # Filter out common HETATMs
        real_unknowns = [u for u in unknowns if u not in ("HOH", "WAT", "NA", "CL", "MG", "ZN", "CA", "FE", "SO4", "PO4", "GOL", "EDO")]
        results["LAW-170"] = (
            "PASS" if len(real_unknowns) == 0 else "FAIL",
            f"{len(real_unknowns)} non-standard residues",
            "STANDARD:20", "deterministic",
        )

        # ════════════════════════════════════════════════════
        # LAW-182: Hydrophobic Burial Ratio [HEURISTIC]
        # ════════════════════════════════════════════════════
        hydrophobic = [r for r in core if r["_name"] in HYDROPHOBIC_RESIDUES and "CA" in r["_atoms"]]
        if len(hydrophobic) >= 5 and len(all_coords) > 20:
            # Compute center of mass
            com = np.mean(all_coords, axis=0)
            all_ca = [r["_atoms"]["CA"] for r in residues if "CA" in r["_atoms"]]
            if all_ca:
                max_dist = max(_dist(ca, com) for ca in all_ca)
                if max_dist > 0:
                    buried = sum(1 for r in hydrophobic if _dist(r["_atoms"]["CA"], com) < max_dist * 0.6)
                    ratio = buried / len(hydrophobic) if hydrophobic else 0
                else:
                    ratio = 0
            else:
                ratio = 0
        else:
            ratio = 0.5  # Small structures: assume OK
        results["LAW-182"] = (
            "PASS" if ratio > 0.3 else "FAIL",
            f"Burial ratio: {round(ratio, 2)}",
            "THRESH:0.3", "heuristic",
        )

        # ════════════════════════════════════════════════════
        # LAW-195: Disulfide Geometry
        # ════════════════════════════════════════════════════
        ss_violations = 0
        cys_residues = [r for r in residues if r["_name"] == "CYS" and "SG" in r["_atoms"]]
        for i in range(len(cys_residues)):
            for j in range(i + 1, len(cys_residues)):
                d = _dist(cys_residues[i]["_atoms"]["SG"], cys_residues[j]["_atoms"]["SG"])
                if d < 3.0:  # Close enough to be a disulfide
                    if abs(d - 2.033) > 0.20:
                        ss_violations += 1
        results["LAW-195"] = (
            "PASS" if ss_violations == 0 else "VETO",
            f"{ss_violations} S-S geometry violations",
            "THRESH:2.04A+/-0.20", "deterministic",
        )

        # ════════════════════════════════════════════════════
        # LAW-200: Internal Cavity Detection [HEURISTIC]
        # ════════════════════════════════════════════════════
        void_volume = 0.0
        if len(all_coords) >= 50:
            try:
                # Grid-based void estimation
                mins = all_coords.min(axis=0) - 2.0
                maxs = all_coords.max(axis=0) + 2.0
                vol_total = float(np.prod(maxs - mins))
                # Estimate occupied volume (spheres of ~1.7A radius)
                n_atoms = len(all_coords)
                vol_atoms = n_atoms * (4 / 3) * np.pi * (1.7 ** 3)
                void_volume = max(vol_total - vol_atoms, 0)
            except Exception:
                void_volume = 0.0
        results["LAW-200"] = (
            "PASS" if void_volume < 1000 else "FAIL",
            f"Estimated void: {round(void_volume, 1)} A^3",
            "THRESH:1000A3", "heuristic",
        )

        return results, coverage, fatal_fringe

    # ── Characterization ────────────────────────────────────
    @staticmethod
    def compute_structural_characterization(structure):
        """Compute secondary structure from phi/psi angles."""
        residues = Tier1Measurements._extract(structure)
        total_atoms = len(structure.atoms)
        total_residues = len(residues)

        helix_count, sheet_count, assessed = 0, 0, 0
        for i in range(1, len(residues) - 1):
            p, c, n = residues[i - 1], residues[i], residues[i + 1]
            if not (Tier1Measurements._is_sequential(p, c) and Tier1Measurements._is_sequential(c, n)):
                continue
            if "C" not in p["_atoms"] or "N" not in c["_atoms"] or "CA" not in c["_atoms"] or "C" not in c["_atoms"] or "N" not in n["_atoms"]:
                continue
            phi = _dihedral(p["_atoms"]["C"], c["_atoms"]["N"], c["_atoms"]["CA"], c["_atoms"]["C"])
            psi = _dihedral(c["_atoms"]["N"], c["_atoms"]["CA"], c["_atoms"]["C"], n["_atoms"]["N"])
            if phi is None or psi is None:
                continue
            assessed += 1
            if -120 < phi < -20 and -80 < psi < -10:
                helix_count += 1
            elif -180 < phi < -40 and 70 < psi < 180:
                sheet_count += 1

        if assessed > 0:
            helix_pct = round(helix_count / assessed * 100, 1)
            sheet_pct = round(sheet_count / assessed * 100, 1)
            loop_pct = round(100 - helix_pct - sheet_pct, 1)
        else:
            helix_pct, sheet_pct, loop_pct = 0.0, 0.0, 100.0

        return {
            "helix": helix_pct, "sheet": sheet_pct, "loop": loop_pct,
            "total_atoms": total_atoms, "total_residues": total_residues,
            "method": "phi_psi_rect",
        }

    # ── Confidence ──────────────────────────────────────────
    @staticmethod
    def detect_confidence_source(structure):
        """Detect pLDDT vs B-factor from column 61-66 values."""
        m = structure.confidence.mean_plddt
        if m <= 0:
            return ("none", 0.0, "missing")
        b_values = [a.b_iso for a in structure.atoms if a.b_iso is not None]
        if not b_values:
            return ("none", 0.0, "missing")
        max_b = max(b_values)
        if max_b > 100:
            return ("B-factor", m, "heuristic_bfactor")
        if m > 50 and max_b <= 100:
            return ("pLDDT", m, "heuristic_plddt")
        return ("B-factor", m, "ambiguous")

    @staticmethod
    def decompose_confidence(structure):
        v = [a.b_iso for a in structure.atoms]
        if v:
            return {"mean": round(float(np.mean(v)), 1), "data_available": True}
        return {"mean": 0, "data_available": False}
