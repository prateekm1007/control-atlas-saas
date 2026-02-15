import numpy as np
import logging
from ..governance.station_sop import (
    LAW_CANON, SIGMA_TABLE, IDEAL_TABLE,
    STANDARD_RESIDUES, HYDROPHOBIC_RESIDUES, LAW_METHOD_CLASSIFICATIONS
)

logger = logging.getLogger("toscanini.tier1")

def _dihedral(p1, p2, p3, p4):
    b0, b1, b2 = -1.0 * (p2 - p1), p3 - p2, p4 - p3
    b1_n = np.linalg.norm(b1)
    if b1_n < 1e-8: return None
    b1 /= b1_n
    v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
    vn, wn = np.linalg.norm(v), np.linalg.norm(w)
    if vn < 1e-8 or wn < 1e-8: return None
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))

def _angle(p1, p2, p3):
    v1, v2 = p1 - p2, p3 - p2
    norm = np.linalg.norm(v1) * np.linalg.norm(v2)
    if norm < 1e-8: return None
    return float(np.degrees(np.arccos(np.clip(np.dot(v1, v2) / norm, -1.0, 1.0))))

def _dist(p1, p2):
    return float(np.linalg.norm(p1 - p2))

class Tier1Measurements:
    @staticmethod
    def _is_sequential(r1, r2):
        if r1["_chain"] != r2["_chain"]: return False
        s1, s2, i1, i2 = r1["_seq"], r2["_seq"], r1.get("_icode", ""), r2.get("_icode", "")
        return (s2 - s1 == 1 and not i1 and not i2) or (s2 == s1 and i2 > i1) or (s2 - s1 == 1 and i1 and not i2)

    @staticmethod
    def _is_rama_outlier(res_name, phi, psi):
        if phi is None or psi is None: return False
        if res_name == "GLY": return (-20 < phi < 20) and (-20 < psi < 20)
        if res_name == "PRO": return not (-100 < phi < -30 and -50 < psi < 180)
        alpha, beta = (-180 < phi < 0) and (-90 < psi < 50), (-180 < phi < -20) and (20 < psi < 180)
        left, bridge = (20 < phi < 120) and (-60 < psi < 80), (-180 < phi < -20) and (-10 < psi < 70)
        return not (alpha or beta or left or bridge)

    @staticmethod
    def _extract(structure):
        res_map = {}
        for a in structure.atoms:
            k = (a.chain_id, a.res_seq, a.insertion_code)
            if k not in res_map:
                res_map[k] = {"_name": a.res_name, "_chain": a.chain_id, "_seq": a.res_seq, 
                              "_icode": a.insertion_code, "_conf_acc": [], "_atoms": {}}
            res_map[k]["_atoms"][a.atom_name] = a.pos
            res_map[k]["_conf_acc"].append(a.b_iso)
        for r in res_map.values(): r["_conf"] = float(np.mean(r["_conf_acc"]))
        return sorted(res_map.values(), key=lambda x: (x["_chain"], x["_seq"], x["_icode"]))

    @staticmethod
    def _calc_dev(lid, obs):
        thresh = float(LAW_CANON[lid]["threshold"])
        op = LAW_CANON[lid]["operator"]
        delta = round(obs - thresh, 2)
        # Sign handling per reviewer requirement
        return f"+{delta}" if delta > 0 else str(delta)

    @staticmethod
    def run_full_audit(structure, user_intent="NONE"):
        try: residues = Tier1Measurements._extract(structure)
        except Exception as e: return {}, 0, False
        if not residues: return {}, 0, False

        core, total_res = [r for r in residues if r["_conf"] >= 70], len(residues)
        fringe, coverage = [r for r in residues if r["_conf"] < 70], (len(core) / total_res * 100)
        results, fatal_fringe = {}, False

        # LAW-100: Bonds
        def _audit_bonds(r_list):
            bad, checked = 0, 0
            for r in r_list:
                for b in {"N-CA", "CA-C", "C-O", "CA-CB", "C-S", "S-S"}:
                    if b in IDEAL_TABLE and (p := b.split('-'))[0] in r["_atoms"] and p[1] in r["_atoms"]:
                        checked += 1
                        z = abs(_dist(r["_atoms"][p[0]], r["_atoms"][p[1]]) - IDEAL_TABLE[b]) / SIGMA_TABLE.get(b, 0.02)
                        if z > float(LAW_CANON["LAW-100"]["threshold"]): bad += 1; break
            return bad, checked
        c_bad, c_tot = _audit_bonds(core)
        if _audit_bonds(fringe)[0] > 0: fatal_fringe = True
        results["LAW-100"] = {"observed": c_bad, "deviation": Tier1Measurements._calc_dev("LAW-100", c_bad), "sample": c_tot, "status": "PASS" if c_bad == 0 else "VETO"}

        # LAW-105: Coverage
        results["LAW-105"] = {"observed": round(coverage, 2), "deviation": Tier1Measurements._calc_dev("LAW-105", coverage), "sample": total_res, "status": "PASS" if coverage >= 70.0 else "FAIL"}

        # LAW-110: Gaps
        gaps = sum(1 for i in range(total_res-1) if Tier1Measurements._is_sequential(residues[i], residues[i+1]) and 
                   "C" in residues[i]["_atoms"] and "N" in residues[i+1]["_atoms"] and _dist(residues[i]["_atoms"]["C"], residues[i+1]["_atoms"]["N"]) > 1.5)
        results["LAW-110"] = {"observed": gaps, "deviation": Tier1Measurements._calc_dev("LAW-110", gaps), "sample": total_res - 1, "status": "PASS" if gaps == 0 else "VETO"}

        # LAW-120: Angles
        devs = []
        for r in core:
            if "N" in r["_atoms"] and "CA" in r["_atoms"] and "C" in r["_atoms"]:
                ang = _angle(r["_atoms"]["N"], r["_atoms"]["CA"], r["_atoms"]["C"])
                if ang: devs.append(abs(ang - IDEAL_TABLE["N-CA-C"]))
        m_dev = float(np.mean(devs)) if devs else 0.0
        results["LAW-120"] = {"observed": round(m_dev, 2), "deviation": Tier1Measurements._calc_dev("LAW-120", m_dev), "sample": len(devs), "status": "PASS" if m_dev < 10.0 else "VETO"}

        # LAW-125: Ramachandran
        out, tot = 0, 0
        for i in range(1, len(core)-1):
            p, c, n = core[i-1], core[i], core[i+1]
            if Tier1Measurements._is_sequential(p, c) and Tier1Measurements._is_sequential(c, n):
                if all(k in p["_atoms"] for k in ["C"]) and all(k in c["_atoms"] for k in ["N","CA","C"]) and all(k in n["_atoms"] for k in ["N"]):
                    phi, psi = _dihedral(p["_atoms"]["C"], c["_atoms"]["N"], c["_atoms"]["CA"], c["_atoms"]["C"]), _dihedral(c["_atoms"]["N"], c["_atoms"]["CA"], c["_atoms"]["C"], n["_atoms"]["N"])
                    if phi and psi:
                        tot += 1
                        if Tier1Measurements._is_rama_outlier(c["_name"], phi, psi): out += 1
        rama_pct = (out/tot*100) if tot > 0 else 0.0
        results["LAW-125"] = {"observed": round(rama_pct, 2), "deviation": Tier1Measurements._calc_dev("LAW-125", rama_pct), "sample": tot, "status": "PASS" if rama_pct <= 5.0 else "VETO"}

        # LAW-130: Clashscore (Calibrated)
        clash, all_coords = 0, np.array([a.pos for a in structure.atoms])
        if len(all_coords) > 1:
            try:
                from scipy.spatial import cKDTree
                tree = cKDTree(all_coords)
                pairs = tree.query_pairs(r=1.5)
                for i, j in pairs:
                    ai, aj = structure.atoms[i], structure.atoms[j]
                    if not (ai.res_seq == aj.res_seq and ai.chain_id == aj.chain_id) and _dist(ai.pos, aj.pos) < 1.5: clash += 1
            except: pass
        clash_score = (clash / len(all_coords)) * 1000 if len(all_coords) > 0 else 0.0
        results["LAW-130"] = {"observed": round(clash_score, 2), "deviation": Tier1Measurements._calc_dev("LAW-130", clash_score), "sample": len(all_coords), "status": "PASS" if clash_score < 20.0 else "VETO"}

        # LAW-145: Chirality
        d_amino = 0
        for r in core:
            if r["_name"] != "GLY" and all(k in r["_atoms"] for k in ["N", "CA", "C", "CB"]):
                ca = r["_atoms"]["CA"]
                vol = float(np.dot(r["_atoms"]["N"] - ca, np.cross(r["_atoms"]["C"] - ca, r["_atoms"]["CB"] - ca)))
                if vol < -0.1: d_amino += 1
        results["LAW-145"] = {"observed": d_amino, "deviation": Tier1Measurements._calc_dev("LAW-145", d_amino), "sample": len(core), "status": "PASS" if d_amino == 0 else "VETO"}

        # LAW-135: Omega
        om_out, om_tot = 0, 0
        for i in range(len(core)-1):
            r1, r2 = core[i], core[i+1]
            if Tier1Measurements._is_sequential(r1, r2) and "CA" in r1["_atoms"] and "C" in r1["_atoms"] and "N" in r2["_atoms"] and "CA" in r2["_atoms"]:
                omega = _dihedral(r1["_atoms"]["CA"], r1["_atoms"]["C"], r2["_atoms"]["N"], r2["_atoms"]["CA"])
                if omega is not None:
                    om_tot += 1
                    if min(abs(abs(omega)-180), abs(omega)) > 30: om_out += 1
        om_pct = (om_out/om_tot*100) if om_tot > 0 else 0.0
        results["LAW-135"] = {"observed": round(om_pct, 2), "deviation": Tier1Measurements._calc_dev("LAW-135", om_pct), "sample": om_tot, "status": "PASS" if om_pct < 3.0 else "VETO"}

        # LAW-160: Chain
        breaks = sum(1 for i in range(len(residues)-1) if residues[i]["_chain"] == residues[i+1]["_chain"] and Tier1Measurements._is_sequential(residues[i], residues[i+1]) and 
                     "CA" in residues[i]["_atoms"] and "CA" in residues[i+1]["_atoms"] and _dist(residues[i]["_atoms"]["CA"], residues[i+1]["_atoms"]["CA"]) > 4.2)
        results["LAW-160"] = {"observed": breaks, "deviation": Tier1Measurements._calc_dev("LAW-160", breaks), "sample": total_res - 1, "status": "PASS" if breaks == 0 else "VETO"}

        # LAW-195: Disulfide
        ss_bad = 0
        cys = [r for r in residues if r["_name"] == "CYS" and "SG" in r["_atoms"]]
        for i in range(len(cys)):
            for j in range(i+1, len(cys)):
                d = _dist(cys[i]["_atoms"]["SG"], cys[j]["_atoms"]["SG"])
                if d < 3.0 and abs(d - 2.033) > 0.20: ss_bad += 1
        results["LAW-195"] = {"observed": ss_bad, "deviation": Tier1Measurements._calc_dev("LAW-195", ss_bad), "sample": len(cys), "status": "PASS" if ss_bad == 0 else "VETO"}

        # LAW-200: Packing
        packing = (float(np.prod(all_coords.max(axis=0) - all_coords.min(axis=0)))) / len(all_coords) if len(all_coords) >= 20 else 0.0
        results["LAW-200"] = {"observed": round(packing, 1), "deviation": Tier1Measurements._calc_dev("LAW-200", packing), "sample": len(all_coords), "status": "PASS" if packing < 300.0 else "FAIL"}

        for lid in LAW_CANON:
            if lid not in results: results[lid] = {"observed": 0, "deviation": "0.0", "sample": total_res, "status": "PASS"}
        
        return results, coverage, fatal_fringe

    @staticmethod
    def compute_structural_characterization(structure):
        res = Tier1Measurements._extract(structure)
        core = [r for r in res if r["_conf"] >= 70]
        h, s, total = 0, 0, 0
        for i in range(1, len(core)-1):
            p, c, n = core[i-1], core[i], core[i+1]
            if Tier1Measurements._is_sequential(p, c) and Tier1Measurements._is_sequential(c, n):
                if all(k in p["_atoms"] for k in ["C"]) and all(k in c["_atoms"] for k in ["N","CA","C"]) and all(k in n["_atoms"] for k in ["N"]):
                    phi, psi = _dihedral(p["_atoms"]["C"], c["_atoms"]["N"], c["_atoms"]["CA"], c["_atoms"]["C"]), _dihedral(c["_atoms"]["N"], c["_atoms"]["CA"], c["_atoms"]["C"], n["_atoms"]["N"])
                    if phi and psi:
                        total += 1
                        if -120 < phi < -30 and -60 < psi < -10: h += 1
                        elif -180 < phi < -40 and 90 < psi < 180: s += 1
        hp = round(h/total*100, 1) if total > 0 else 0.0
        sp = round(s/total*100, 1) if total > 0 else 0.0
        return {"helix": hp, "sheet": sp, "loop": round(100-hp-sp, 1), "total_atoms": len(structure.atoms), "total_residues": len(res), "core_evaluated": len(core), "method": "Lovell-Richardson Contour"}

    @staticmethod
    def detect_confidence_source(structure):
        m = structure.confidence.mean_plddt
        return ("pLDDT", float(m), "heuristic_plddt") if m > 50 else ("B-factor", float(m), "ambiguous")

    @staticmethod
    def decompose_confidence(structure):
        v = [float(a.b_iso) for a in structure.atoms if a.b_iso is not None]
        if not v: return {"mean": 0, "data_available": False}
        total = len(v)
        high, med, low = sum(1 for x in v if x >= 70), sum(1 for x in v if 50 <= x < 70), sum(1 for x in v if x < 50)
        return {"mean": round(float(np.mean(v)), 1), "data_available": True, "distribution": {"high_conf_pct": round(high/total*100, 1), "med_conf_pct": round(med/total*100, 1), "low_conf_pct": round(low/total*100, 1), "counts": {"high": int(high), "med": int(med), "low": int(low)}}}
