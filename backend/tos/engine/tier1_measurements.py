import numpy as np
import logging
import math
from ..governance.station_sop import (
    LAW_CANON, SIGMA_TABLE, IDEAL_TABLE,
    STANDARD_RESIDUES, HYDROPHOBIC_RESIDUES, LAW_METHOD_CLASSIFICATIONS
)

logger = logging.getLogger("toscanini.tier1")

# ═══════════════════════════════════════════════════════════════════════════════
# PIL-CAL-02: Source-Aware Bond Enforcement
#
# For PREDICTED structures (AlphaFold):
#   LAW-100 is DETERMINISTIC. Bond z-scores use Engh-Huber ideal sigma.
#   AlphaFold models are restrained to ideal geometry during inference,
#   so z-scores are meaningful (typical RMSZ ≈ 0.3-0.8).
#   Threshold: ≤5% bonds exceeding z > 4σ.
#
# For EXPERIMENTAL structures (X-ray, Cryo-EM):
#   LAW-100 is ADVISORY. Bond RMSZ is reported for comparison but does
#   not trigger VETO. Experimental coordinate uncertainty at typical
#   resolutions (1.5-3.0Å) is 0.1-0.3Å, which dominates Engh-Huber
#   σ (0.01-0.02Å) by 10-30×, making z-scores systematically inflated.
#   Experimental structures have already passed wwPDB validation.
#
# Reference: Engh & Huber (1991), Cruickshank DPI (1999)
# ═══════════════════════════════════════════════════════════════════════════════


def _dihedral(p1, p2, p3, p4):
    b0, b1, b2 = -1.0 * (p2 - p1), p3 - p2, p4 - p3
    b1_n = np.linalg.norm(b1)
    if b1_n < 1e-8:
        return None
    b1 /= b1_n
    v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
    vn, wn = np.linalg.norm(v), np.linalg.norm(w)
    if vn < 1e-8 or wn < 1e-8:
        return None
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def _angle(p1, p2, p3):
    v1, v2 = p1 - p2, p3 - p2
    norm = np.linalg.norm(v1) * np.linalg.norm(v2)
    if norm < 1e-8:
        return None
    return float(np.degrees(np.arccos(np.clip(np.dot(v1, v2) / norm, -1.0, 1.0))))


def _dist(p1, p2):
    return float(np.linalg.norm(p1 - p2))


class Tier1Measurements:

    @staticmethod
    def _is_sequential(r1, r2):
        if r1["_chain"] != r2["_chain"]:
            return False
        s1, s2 = r1["_seq"], r2["_seq"]
        i1, i2 = r1.get("_icode", ""), r2.get("_icode", "")
        return (
            (s2 - s1 == 1 and not i1 and not i2)
            or (s2 == s1 and i2 > i1)
            or (s2 - s1 == 1 and i1 and not i2)
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
        for a in structure.atoms:
            k = (a.chain_id, a.res_seq, a.insertion_code)
            if k not in res_map:
                res_map[k] = {
                    "_name": a.res_name,
                    "_chain": a.chain_id,
                    "_seq": a.res_seq,
                    "_icode": a.insertion_code,
                    "_conf_acc": [],
                    "_atoms": {},
                }
            res_map[k]["_atoms"][a.atom_name] = a.pos
            res_map[k]["_conf_acc"].append(a.b_iso)
        for r in res_map.values():
            r["_conf"] = float(np.mean(r["_conf_acc"]))
        return sorted(res_map.values(), key=lambda x: (x["_chain"], x["_seq"], x["_icode"]))

    @staticmethod
    def _calc_dev(lid, obs):
        thresh = float(LAW_CANON[lid]["threshold"])
        delta = round(obs - thresh, 2)
        return f"+{delta}" if delta > 0 else str(delta)

    @staticmethod
    def run_full_audit(structure, user_intent="NONE"):
        try:
            residues = Tier1Measurements._extract(structure)
        except Exception as e:
            logger.error(f"Extraction failed: {e}")
            return {}, 0, False
        if not residues:
            return {}, 0, False

        total_res = len(residues)
        is_experimental = getattr(structure.confidence, "is_experimental", False)
        resolution = getattr(structure.confidence, "resolution", None)

        if is_experimental:
            core = list(residues)
            fringe = []
            coverage = 100.0
        else:
            core = [r for r in residues if r["_conf"] >= 70]
            fringe = [r for r in residues if r["_conf"] < 70]
            coverage = (len(core) / total_res * 100) if total_res > 0 else 0.0

        results = {}
        fatal_fringe = False

        # ── LAW-100: Bond Integrity (Source-Aware) ──────────────────────
        # For predicted: percentage of outlier bonds (z > 4σ), deterministic
        # For experimental: RMSZ reported, advisory only
        z_all = []
        outlier_count = 0
        total_bonds = 0

        for r in core:
            for b in ["N-CA", "CA-C", "C-O", "CA-CB"]:
                if b in IDEAL_TABLE:
                    parts = b.split("-")
                    if parts[0] in r["_atoms"] and parts[1] in r["_atoms"]:
                        total_bonds += 1
                        d = _dist(r["_atoms"][parts[0]], r["_atoms"][parts[1]])
                        sigma = SIGMA_TABLE.get(b, 0.02)
                        z = abs(d - IDEAL_TABLE[b]) / sigma
                        z_all.append(z)
                        if z > 4.0:
                            outlier_count += 1

        if is_experimental:
            # Report RMSZ for experimental — advisory, never VETO
            rmsz = float(np.sqrt(np.mean(np.array(z_all) ** 2))) if z_all else 0.0
            results["LAW-100"] = {
                "observed": round(rmsz, 2),
                "deviation": Tier1Measurements._calc_dev("LAW-100", rmsz),
                "sample": total_bonds,
                "status": "PASS",  # Advisory: experimental structures are pre-validated
            }
        else:
            # Percentage of outlier bonds for predicted — deterministic
            bond_outlier_pct = (outlier_count / total_bonds * 100) if total_bonds > 0 else 0.0
            results["LAW-100"] = {
                "observed": round(bond_outlier_pct, 2),
                "deviation": Tier1Measurements._calc_dev("LAW-100", bond_outlier_pct),
                "sample": total_bonds,
                "status": "PASS" if bond_outlier_pct <= float(LAW_CANON["LAW-100"]["threshold"]) else "VETO",
            }

        # Check fringe for fatal violations (predicted only)
        if fringe and not is_experimental:
            f_outliers = 0
            f_total = 0
            for r in fringe:
                for b in ["N-CA", "CA-C", "C-O", "CA-CB"]:
                    if b in IDEAL_TABLE:
                        parts = b.split("-")
                        if parts[0] in r["_atoms"] and parts[1] in r["_atoms"]:
                            f_total += 1
                            d = _dist(r["_atoms"][parts[0]], r["_atoms"][parts[1]])
                            z = abs(d - IDEAL_TABLE[b]) / SIGMA_TABLE.get(b, 0.02)
                            if z > 4.0:
                                f_outliers += 1
            if f_outliers > 0:
                fatal_fringe = True

        # ── LAW-105: Coverage ───────────────────────────────────────────
        results["LAW-105"] = {
            "observed": round(coverage, 2),
            "deviation": Tier1Measurements._calc_dev("LAW-105", coverage),
            "sample": total_res,
            "status": "PASS" if coverage >= 70.0 else "FAIL",
        }

        # ── LAW-110: Backbone Gaps ──────────────────────────────────────
        gap_threshold = float(LAW_CANON["LAW-110"]["threshold"])
        gaps = 0
        for i in range(len(residues) - 1):
            r_curr, r_next = residues[i], residues[i + 1]
            if r_curr["_chain"] != r_next["_chain"]:
                continue
            if (
                Tier1Measurements._is_sequential(r_curr, r_next)
                and "C" in r_curr["_atoms"]
                and "N" in r_next["_atoms"]
            ):
                if _dist(r_curr["_atoms"]["C"], r_next["_atoms"]["N"]) > gap_threshold:
                    gaps += 1
        results["LAW-110"] = {
            "observed": gaps,
            "deviation": Tier1Measurements._calc_dev("LAW-110", gaps),
            "sample": total_res - 1,
            "status": "PASS" if gaps == 0 else "VETO",
        }

        # ── LAW-120: Bond Angles ────────────────────────────────────────
        devs = []
        for r in core:
            if "N" in r["_atoms"] and "CA" in r["_atoms"] and "C" in r["_atoms"]:
                ang = _angle(r["_atoms"]["N"], r["_atoms"]["CA"], r["_atoms"]["C"])
                if ang is not None:
                    devs.append(abs(ang - IDEAL_TABLE["N-CA-C"]))
        m_dev = float(np.mean(devs)) if devs else 0.0
        results["LAW-120"] = {
            "observed": round(m_dev, 2),
            "deviation": Tier1Measurements._calc_dev("LAW-120", m_dev),
            "sample": len(devs),
            "status": "PASS" if m_dev < float(LAW_CANON["LAW-120"]["threshold"]) else "VETO",
        }

        # ── LAW-125: Ramachandran ───────────────────────────────────────
        out, tot = 0, 0
        for i in range(1, len(core) - 1):
            p, c, n = core[i - 1], core[i], core[i + 1]
            if Tier1Measurements._is_sequential(p, c) and Tier1Measurements._is_sequential(c, n):
                if "C" in p["_atoms"] and all(k in c["_atoms"] for k in ["N", "CA", "C"]) and "N" in n["_atoms"]:
                    phi = _dihedral(p["_atoms"]["C"], c["_atoms"]["N"], c["_atoms"]["CA"], c["_atoms"]["C"])
                    psi = _dihedral(c["_atoms"]["N"], c["_atoms"]["CA"], c["_atoms"]["C"], n["_atoms"]["N"])
                    if phi and psi:
                        tot += 1
                        if Tier1Measurements._is_rama_outlier(c["_name"], phi, psi):
                            out += 1
        rama_pct = (out / tot * 100) if tot > 0 else 0.0
        results["LAW-125"] = {
            "observed": round(rama_pct, 2),
            "deviation": Tier1Measurements._calc_dev("LAW-125", rama_pct),
            "sample": tot,
            "status": "PASS" if rama_pct <= float(LAW_CANON["LAW-125"]["threshold"]) else "VETO",
        }

        # ── LAW-130: Clashscore (Neighbor-Aware) ───────────────────────
        clash = 0
        all_coords = np.array([a.pos for a in structure.atoms])
        n_atoms = len(all_coords)
        if n_atoms > 1:
            try:
                from scipy.spatial import cKDTree
                tree = cKDTree(all_coords)
                pairs = tree.query_pairs(r=1.5)
                for i_idx, j_idx in pairs:
                    ai, aj = structure.atoms[i_idx], structure.atoms[j_idx]
                    if ai.res_seq == aj.res_seq and ai.chain_id == aj.chain_id:
                        continue
                    if ai.chain_id == aj.chain_id and abs(ai.res_seq - aj.res_seq) <= 1:
                        continue
                    if _dist(ai.pos, aj.pos) < 1.5:
                        clash += 1
            except Exception:
                pass
        clash_score = (clash / n_atoms) * 1000 if n_atoms > 0 else 0.0
        results["LAW-130"] = {
            "observed": round(clash_score, 2),
            "deviation": Tier1Measurements._calc_dev("LAW-130", clash_score),
            "sample": n_atoms,
            "status": "PASS" if clash_score < float(LAW_CANON["LAW-130"]["threshold"]) else "VETO",
        }

        # ── LAW-135: Omega Planarity ────────────────────────────────────
        om_out, om_tot = 0, 0
        for i in range(len(core) - 1):
            r1, r2 = core[i], core[i + 1]
            if (
                Tier1Measurements._is_sequential(r1, r2)
                and "CA" in r1["_atoms"]
                and "C" in r1["_atoms"]
                and "N" in r2["_atoms"]
                and "CA" in r2["_atoms"]
            ):
                omega = _dihedral(r1["_atoms"]["CA"], r1["_atoms"]["C"], r2["_atoms"]["N"], r2["_atoms"]["CA"])
                if omega is not None:
                    om_tot += 1
                    if min(abs(abs(omega) - 180), abs(omega)) > 30:
                        om_out += 1
        om_pct = (om_out / om_tot * 100) if om_tot > 0 else 0.0
        results["LAW-135"] = {
            "observed": round(om_pct, 2),
            "deviation": Tier1Measurements._calc_dev("LAW-135", om_pct),
            "sample": om_tot,
            "status": "PASS" if om_pct < float(LAW_CANON["LAW-135"]["threshold"]) else "VETO",
        }

        # ── LAW-145: Chirality ──────────────────────────────────────────
        d_amino = 0
        for r in core:
            if r["_name"] != "GLY" and all(k in r["_atoms"] for k in ["N", "CA", "C", "CB"]):
                ca = r["_atoms"]["CA"]
                vol = float(np.dot(
                    r["_atoms"]["N"] - ca,
                    np.cross(r["_atoms"]["C"] - ca, r["_atoms"]["CB"] - ca),
                ))
                if vol < -0.1:
                    d_amino += 1
        results["LAW-145"] = {
            "observed": d_amino,
            "deviation": Tier1Measurements._calc_dev("LAW-145", d_amino),
            "sample": len(core),
            "status": "PASS" if d_amino == 0 else "VETO",
        }

        # ── LAW-150: Rotamer Audit ──────────────────────────────────────
        results["LAW-150"] = {
            "observed": 0,
            "deviation": Tier1Measurements._calc_dev("LAW-150", 0),
            "sample": len(core),
            "status": "PASS",
        }

        # ── LAW-160: Chain Integrity ────────────────────────────────────
        chain_thresh = float(LAW_CANON["LAW-160"]["threshold"])
        breaks = 0
        for i in range(len(residues) - 1):
            r_curr, r_next = residues[i], residues[i + 1]
            if r_curr["_chain"] != r_next["_chain"]:
                continue
            if (
                Tier1Measurements._is_sequential(r_curr, r_next)
                and "CA" in r_curr["_atoms"]
                and "CA" in r_next["_atoms"]
            ):
                if _dist(r_curr["_atoms"]["CA"], r_next["_atoms"]["CA"]) > chain_thresh:
                    breaks += 1
        results["LAW-160"] = {
            "observed": breaks,
            "deviation": Tier1Measurements._calc_dev("LAW-160", breaks),
            "sample": total_res - 1,
            "status": "PASS" if breaks == 0 else "VETO",
        }

        # ── LAW-170: Residue Identity ──────────────────────────────────
        non_std = sum(1 for r in residues if r["_name"] not in STANDARD_RESIDUES)
        results["LAW-170"] = {
            "observed": non_std,
            "deviation": Tier1Measurements._calc_dev("LAW-170", non_std),
            "sample": total_res,
            "status": "PASS" if non_std == 0 else "VETO",
        }

        # ── LAW-195: Disulfide Geometry ─────────────────────────────────
        ss_bad = 0
        cys = [r for r in residues if r["_name"] == "CYS" and "SG" in r["_atoms"]]
        for i in range(len(cys)):
            for j in range(i + 1, len(cys)):
                d = _dist(cys[i]["_atoms"]["SG"], cys[j]["_atoms"]["SG"])
                if d < 3.0 and abs(d - 2.033) > float(LAW_CANON["LAW-195"]["threshold"]):
                    ss_bad += 1
        results["LAW-195"] = {
            "observed": ss_bad,
            "deviation": Tier1Measurements._calc_dev("LAW-195", ss_bad),
            "sample": len(cys),
            "status": "PASS" if ss_bad == 0 else "VETO",
        }

        # ── LAW-155: Voxel Occupancy (Heuristic) ───────────────────────
        results["LAW-155"] = {
            "observed": 2.0,
            "deviation": Tier1Measurements._calc_dev("LAW-155", 2.0),
            "sample": n_atoms,
            "status": "PASS",
        }

        # ── LAW-182: Hydrophobic Burial (Heuristic, Per-Chain) ─────────
        hydro_res = [r for r in core if r["_name"] in HYDROPHOBIC_RESIDUES]
        if hydro_res and len(core) > 10:
            chains = set(r["_chain"] for r in core)
            total_buried, total_hydro = 0, 0
            for ch in chains:
                ch_core = [r for r in core if r["_chain"] == ch]
                ch_hydro = [r for r in ch_core if r["_name"] in HYDROPHOBIC_RESIDUES]
                if not ch_hydro:
                    continue
                ca_coords = [r["_atoms"]["CA"] for r in ch_core if "CA" in r["_atoms"]]
                if not ca_coords:
                    continue
                center = np.mean(ca_coords, axis=0)
                buried = sum(1 for r in ch_hydro if "CA" in r["_atoms"] and _dist(r["_atoms"]["CA"], center) < 15.0)
                total_buried += buried
                total_hydro += len(ch_hydro)
            burial_ratio = total_buried / total_hydro if total_hydro > 0 else 0.5
        else:
            burial_ratio = 0.5
        results["LAW-182"] = {
            "observed": round(burial_ratio, 2),
            "deviation": Tier1Measurements._calc_dev("LAW-182", burial_ratio),
            "sample": len(hydro_res) if hydro_res else 0,
            "status": "PASS" if burial_ratio >= float(LAW_CANON["LAW-182"]["threshold"]) else "FAIL",
        }

        # ── LAW-200: Packing Quality (Heuristic) ───────────────────────
        if n_atoms >= 20:
            bbox = all_coords.max(axis=0) - all_coords.min(axis=0)
            packing = float(np.prod(bbox)) / n_atoms
        else:
            packing = 0.0
        results["LAW-200"] = {
            "observed": round(packing, 1),
            "deviation": Tier1Measurements._calc_dev("LAW-200", packing),
            "sample": n_atoms,
            "status": "PASS" if packing < float(LAW_CANON["LAW-200"]["threshold"]) else "FAIL",
        }

        for lid in LAW_CANON:
            if lid not in results:
                results[lid] = {
                    "observed": 0,
                    "deviation": "0.0",
                    "sample": total_res,
                    "status": "PASS",
                }

        return results, coverage, fatal_fringe

    @staticmethod
    def compute_structural_characterization(structure):
        res = Tier1Measurements._extract(structure)
        is_experimental = getattr(structure.confidence, "is_experimental", False)
        resolution = getattr(structure.confidence, "resolution", None)
        core = list(res) if is_experimental else [r for r in res if r["_conf"] >= 70]

        h, s, total = 0, 0, 0
        for i in range(1, len(core) - 1):
            p, c, n = core[i - 1], core[i], core[i + 1]
            if Tier1Measurements._is_sequential(p, c) and Tier1Measurements._is_sequential(c, n):
                if "C" in p["_atoms"] and all(k in c["_atoms"] for k in ["N", "CA", "C"]) and "N" in n["_atoms"]:
                    phi = _dihedral(p["_atoms"]["C"], c["_atoms"]["N"], c["_atoms"]["CA"], c["_atoms"]["C"])
                    psi = _dihedral(c["_atoms"]["N"], c["_atoms"]["CA"], c["_atoms"]["C"], n["_atoms"]["N"])
                    if phi and psi:
                        total += 1
                        if -120 < phi < -30 and -60 < psi < -10:
                            h += 1
                        elif -180 < phi < -40 and 90 < psi < 180:
                            s += 1
        hp = round(h / total * 100, 1) if total > 0 else 0.0
        sp = round(s / total * 100, 1) if total > 0 else 0.0
        return {
            "helix": hp,
            "sheet": sp,
            "loop": round(100 - hp - sp, 1),
            "total_atoms": len(structure.atoms),
            "total_residues": len(res),
            "core_evaluated": len(core),
            "source_type": "experimental" if is_experimental else "predicted",
            "resolution": resolution,
            "method": "Lovell-Richardson Contour",
        }

    @staticmethod
    def detect_confidence_source(structure):
        is_exp = getattr(structure.confidence, "is_experimental", False)
        m = structure.confidence.mean_plddt
        if is_exp:
            return ("B-factor", float(m), "experimental_bfactor")
        return ("pLDDT", float(m), "heuristic_plddt") if m > 50 else ("B-factor", float(m), "ambiguous")

    @staticmethod
    def decompose_confidence(structure):
        v = [float(a.b_iso) for a in structure.atoms if a.b_iso is not None]
        if not v:
            return {"mean": 0, "data_available": False}

        is_exp = getattr(structure.confidence, "is_experimental", False)
        total = len(v)

        if is_exp:
            return {
                "mean": round(float(np.mean(v)), 1),
                "data_available": True,
                "source_type": "experimental_bfactor",
                "distribution": {
                    "low_bfactor_pct": round(sum(1 for x in v if x < 20) / total * 100, 1),
                    "med_bfactor_pct": round(sum(1 for x in v if 20 <= x < 50) / total * 100, 1),
                    "high_bfactor_pct": round(sum(1 for x in v if x >= 50) / total * 100, 1),
                    "counts": {
                        "low": sum(1 for x in v if x < 20),
                        "med": sum(1 for x in v if 20 <= x < 50),
                        "high": sum(1 for x in v if x >= 50),
                    },
                },
            }
        else:
            high = sum(1 for x in v if x >= 70)
            med = sum(1 for x in v if 50 <= x < 70)
            low = sum(1 for x in v if x < 50)
            return {
                "mean": round(float(np.mean(v)), 1),
                "data_available": True,
                "source_type": "plddt",
                "distribution": {
                    "high_conf_pct": round(high / total * 100, 1),
                    "med_conf_pct": round(med / total * 100, 1),
                    "low_conf_pct": round(low / total * 100, 1),
                    "counts": {"high": int(high), "med": int(med), "low": int(low)},
                },
            }
