import numpy as np
import re as _re
import logging
from ..governance.station_sop import LAW_CANON, STANDARD_RESIDUES, HYDROPHOBIC_RESIDUES, TOTAL_LAWS, LAW_METHOD_CLASSIFICATIONS

logger = logging.getLogger("toscanini.tier1")

def _sanitize_measurement(text):
    s = str(text); s = _re.sub(r'np\.\w+\(([^)]+)\)', r'\1', s); return s

def _make_result(status, measurement, detail, law_id):
    return (status, _sanitize_measurement(measurement), detail, LAW_METHOD_CLASSIFICATIONS.get(law_id, "unknown"))

def _dihedral(p1, p2, p3, p4):
    b0, b1, b2 = -1.0*(p2-p1), p3-p2, p4-p3
    b1 /= (np.linalg.norm(b1) + 1e-12)
    v, w = b0 - np.dot(b0, b1)*b1, b2 - np.dot(b2, b1)*b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))

def _angle(p1, p2, p3):
    v1, v2 = p1-p2, p3-p2
    return float(np.degrees(np.arccos(np.clip(np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)+1e-12), -1.0, 1.0))))

def _are_sequential(r1, r2):
    if r1["_chain"] != r2["_chain"]: return False
    return (r2["_seq"] - r1["_seq"]) == 1 or (r2["_seq"] == r1["_seq"] and r2.get("_icode","") > r1.get("_icode",""))

def _extract_residues(structure):
    res_map = {}
    for atom in structure.atoms:
        key = (atom.chain_id, atom.res_seq, atom.insertion_code)
        if key not in res_map:
            res_map[key] = {"_name": atom.res_name, "_chain": atom.chain_id, "_seq": atom.res_seq, "_icode": atom.insertion_code, "_conf_accum": []}
        res_map[key][atom.atom_name] = atom.pos
        res_map[key]["_conf_accum"].append(atom.b_iso)
    sorted_res = sorted(res_map.values(), key=lambda r: (r["_chain"], r["_seq"], r["_icode"]))
    for r in sorted_res: r["_conf"] = np.mean(r["_conf_accum"])
    return sorted_res

def _build_atom_residue_index(structure):
    atom_data = []
    for atom in structure.atoms:
        atom_data.append((atom.chain_id, atom.res_seq, atom.insertion_code, atom.atom_name, list(atom.pos)))
    atom_data.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
    coords = np.array([a[4] for a in atom_data], dtype=float)
    res_keys = [(a[0], a[1], a[2]) for a in atom_data]
    return coords, res_keys

class Tier1Measurements:
    @staticmethod
    def run_full_audit(structure, user_intent="NONE"):
        """UPDATED SIGNATURE: Accepts user_intent for context-aware laws."""
        results = {}
        try: residues = _extract_residues(structure)
        except Exception as e: return {lid: ("FAIL", f"Error: {e}", "SYSTEM", "error") for lid in LAW_CANON}
        
        # 1. Epistemic Gates
        results["LAW-105"] = Tier1Measurements._check_reliability_coverage(residues, user_intent)
        
        # 2. Geometric Core
        results["LAW-100"] = Tier1Measurements._check_bond_lengths(residues)
        results["LAW-110"] = Tier1Measurements._check_backbone_gaps(residues)
        results["LAW-120"] = Tier1Measurements._check_bond_angles(residues)
        results["LAW-125"] = Tier1Measurements._check_ramachandran(residues)
        results["LAW-130"] = Tier1Measurements._check_steric_clashes(structure)
        results["LAW-135"] = Tier1Measurements._check_omega_planarity(residues)
        
        # 3. Chemical Integrity
        results["LAW-145"] = Tier1Measurements._check_chirality(residues)
        results["LAW-150"] = Tier1Measurements._check_rotamers(residues)
        results["LAW-160"] = Tier1Measurements._check_chain_integrity(residues)
        results["LAW-170"] = Tier1Measurements._check_residue_identity(residues)
        results["LAW-182"] = Tier1Measurements._check_hydrophobic_burial(residues)
        results["LAW-195"] = Tier1Measurements._check_disulfide_geometry(residues)
        results["LAW-200"] = Tier1Measurements._check_internal_cavities(structure)
        
        # 4. Advisory
        results["LAW-155"] = Tier1Measurements._check_voxel_occupancy(structure)
        return results

    @staticmethod
    def _check_reliability_coverage(residues, intent):
        lid = "LAW-105"; high_conf = [r for r in residues if r["_conf"] >= 70.0]
        pct = (len(high_conf)/len(residues))*100 if residues else 0
        # MINI-PROTEIN ADAPTATION
        count_floor = 20 if intent in ["LINEAR", "LINKERS"] else 40
        status = "PASS" if (pct >= 30.0 and len(high_conf) >= count_floor) else "FAIL"
        return _make_result(status, f"Coverage={round(pct,1)}%.", f"GATE:{count_floor}res", lid)

    @staticmethod
    def _check_bond_lengths(residues):
        lid = "LAW-100"; devs = []
        for r in residues:
            for b, i in {"N-CA": 1.47, "CA-C": 1.52, "C-O": 1.23}.items():
                if b.split("-")[0] in r and b.split("-")[1] in r:
                    devs.append(abs(np.linalg.norm(r[b.split("-")[0]]-r[b.split("-")[1]]) - i))
        for i in range(len(residues)-1):
            if _are_sequential(residues[i], residues[i+1]) and "C" in residues[i] and "N" in residues[i+1]:
                d = np.linalg.norm(residues[i]["C"]-residues[i+1]["N"])
                if d < 2.0: devs.append(abs(d-1.33))
        m = round(float(np.mean(devs)), 4) if devs else 0
        return _make_result("PASS" if m < 0.05 else "FAIL", f"Mean={m}A.", "THRESH:0.05A", lid)

    @staticmethod
    def _check_backbone_gaps(residues):
        lid = "LAW-110"; gaps = [r["_seq"] for i, r in enumerate(residues[:-1]) if _are_sequential(r, residues[i+1]) and np.linalg.norm(r["C"]-residues[i+1]["N"]) > 1.5]
        return _make_result("PASS" if not gaps else "FAIL", f"{len(gaps)} gaps.", "GAPS:0", lid)

    @staticmethod
    def _check_bond_angles(residues):
        lid = "LAW-120"; devs = []
        for r in residues:
            if all(a in r for a in ["N","CA","C"]): devs.append(abs(_angle(r["N"], r["CA"], r["C"]) - 111.0))
        for i in range(len(residues)-1):
            if _are_sequential(residues[i], residues[i+1]) and "CA" in residues[i] and "C" in residues[i] and "N" in residues[i+1]:
                devs.append(abs(_angle(residues[i]["CA"], residues[i]["C"], residues[i+1]["N"]) - 116.6))
        m = round(float(np.mean(devs)), 2) if devs else 0
        return _make_result("PASS" if m < 4.6 else "FAIL", f"Mean={m}deg.", "THRESH:4.6deg", lid)

    @staticmethod
    def _check_ramachandran(residues):
        lid = "LAW-125"; outliers, total = [], 0
        for i in range(1, len(residues)-1):
            p, c, n = residues[i-1], residues[i], residues[i+1]
            if not (_are_sequential(p, c) and _are_sequential(c, n)): continue
            if c["_name"] == "GLY" or "C" not in p or not all(a in c for a in ["N","CA","C"]) or "N" not in n: continue
            if c["_conf"] < 70: continue
            phi, psi = _dihedral(p["C"], c["N"], c["CA"], c["C"]), _dihedral(c["N"], c["CA"], c["C"], n["N"])
            alpha, beta = (-180 < phi < 0) and (-90 < psi < 50), (-180 < phi < -20) and (20 < psi < 180)
            left, bridge = (20 < phi < 120) and (-60 < psi < 80), (-180 < phi < -20) and (-10 < psi < 70)
            total += 1
            if not (alpha or beta or left or bridge): outliers.append(c["_seq"])
        if total == 0: return _make_result("FAIL", "Zero high-conf core.", "NO_DATA", lid)
        pct = round((len(outliers)/total)*100, 1)
        return _make_result("PASS" if pct <= 5.0 else "FAIL", f"{pct}% outliers.", "THRESH:5%", lid)

    @staticmethod
    def _check_steric_clashes(structure):
        lid = "LAW-130"; coords, res_keys = _build_atom_residue_index(structure); n = len(coords)
        from scipy.spatial import cKDTree
        tree = cKDTree(coords); pairs = tree.query_pairs(r=2.6); clashes = []
        unique = sorted(list(set(res_keys))); adj = set()
        for i in range(len(unique)-1):
            if unique[i][0] == unique[i+1][0]: adj.add((unique[i], unique[i+1]))
        for i, j in pairs:
            ki, kj = res_keys[i], res_keys[j]
            if ki == kj or (ki, kj) in adj or (kj, ki) in adj: continue
            clashes.append((i,j))
        pct = round((len(clashes)/n)*100, 2)
        return _make_result("PASS" if pct < 0.4 else "FAIL", f"{len(clashes)} clashes ({pct}%).", "THRESH:0.4%", lid)

    @staticmethod
    def _check_omega_planarity(residues):
        lid = "LAW-135"; outliers, total = 0, 0
        for i in range(len(residues)-1):
            r1, r2 = residues[i], residues[i+1]
            if not _are_sequential(r1, r2): continue
            if r2["_conf"] < 70: continue
            if all(a in r1 for a in ["CA","C"]) and all(a in r2 for a in ["N","CA"]):
                omega = _dihedral(r1["CA"], r1["C"], r2["N"], r2["CA"])
                dev = min(abs(omega-180), abs(omega+180)) if r2["_name"] != "PRO" else min(abs(omega-180), abs(omega+180), abs(omega))
                total += 1
                if dev > 30.0: outliers += 1
        if total == 0: return _make_result("FAIL", "Zero high-conf core.", "NO_DATA", lid)
        pct = round((outliers/total)*100, 1)
        return _make_result("PASS" if pct < 3.0 else "FAIL", f"{pct}% outliers.", "THRESH:3%", lid)

    @staticmethod
    def _check_chirality(residues):
        lid = "LAW-145"; v = [r["_seq"] for r in residues if r["_name"]!="GLY" and "CB" in r and np.dot(r["N"]-r["CA"], np.cross(r["C"]-r["CA"], r["CB"]-r["CA"])) < 0]
        return _make_result("PASS" if not v else "FAIL", f"{len(v)} D-amino acids.", "THRESH:L-only", lid)

    @staticmethod
    def _check_rotamers(residues):
        lid = "LAW-150"; o, t = 0, 0
        for r in residues:
            if r["_name"] in ("GLY","ALA") or not all(a in r for a in ["N","CA","CB"]): continue
            g = next((x for x in ["CG","CG1","OG","OG1","SG"] if x in r), None)
            if g:
                t += 1
                if min(abs(_dihedral(r["N"], r["CA"], r["CB"], r[g])-i) for i in [-60, 60, 180]) > 40.0: o += 1
        p = round((o/t)*100, 1) if t else 0
        return _make_result("PASS" if p < 25 else "FAIL", f"{p}% outliers.", "THRESH:25%", lid)

    @staticmethod
    def _check_voxel_occupancy(structure):
        lid = "LAW-155"; coords = np.array([list(a.pos) for a in structure.atoms])
        r = round(len(set(map(tuple, np.floor(coords/2.0).astype(int))))/len(coords), 3)
        return _make_result("PASS" if r < 0.95 else "FAIL", f"Ratio={r}", "CONNECTIVITY", lid)

    @staticmethod
    def _check_chain_integrity(residues):
        lid = "LAW-160"; b = 0
        for i in range(len(residues)-1):
            if _are_sequential(residues[i], residues[i+1]) and "CA" in residues[i] and "CA" in residues[i+1]:
                if np.linalg.norm(residues[i]["CA"]-residues[i+1]["CA"]) > 4.2: b += 1
        return _make_result("PASS" if b == 0 else "FAIL", f"{b} breaks", "THRESH:4.2A", lid)

    @staticmethod
    def _check_residue_identity(residues):
        lid = "LAW-170"; u = [r["_seq"] for r in residues if r["_name"] not in STANDARD_RESIDUES]
        return _make_result("PASS" if not u else "FAIL", f"{len(u)} unknowns", "STANDARD:20", lid)

    @staticmethod
    def _check_hydrophobic_burial(residues):
        lid = "LAW-182"; all_ca = [r["CA"] for r in residues if "CA" in r]
        if len(all_ca) < 10: return _make_result("PASS", "Small", "N/A", lid)
        med = np.median([np.linalg.norm(c - np.mean(all_ca, axis=0)) for c in all_ca])
        b, t = 0, 0
        for r in residues:
            if r["_name"] in HYDROPHOBIC_RESIDUES and "CA" in r:
                t += 1
                if np.linalg.norm(r["CA"]-np.mean(all_ca, axis=0)) < med: b += 1
        ratio = round(b/t, 2) if t else 0
        return _make_result("PASS" if ratio > 0.3 else "FAIL", f"Ratio={ratio}", "THRESH:0.3", lid)

    @staticmethod
    def _check_disulfide_geometry(residues):
        lid = "LAW-195"; cys = [r["SG"] for r in residues if r["_name"]=="CYS" and "SG" in r]; v = 0
        for i in range(len(cys)):
            for j in range(i+1, len(cys)):
                if np.linalg.norm(cys[i]-cys[j]) < 3.0 and abs(np.linalg.norm(cys[i]-cys[j])-2.04) > 0.20: v += 1
        return _make_result("PASS" if v == 0 else "FAIL", f"{v} violations", "THRESH:2.04A", lid)

    @staticmethod
    def _check_internal_cavities(structure):
        lid = "LAW-200"; coords = np.array([list(a.pos) for a in structure.atoms])
        if len(coords) < 50: return _make_result("PASS", "Small", "N/A", lid)
        mins, maxs = coords.min(0)-3, coords.max(0)+3; res=1.8; dims = np.ceil((maxs-mins)/res).astype(int)
        grid = np.zeros(dims, bool)
        for c in coords:
            idx = np.floor((c-mins)/res).astype(int)
            for dx, dy, dz in np.ndindex((3,3,3)):
                ni = (idx[0]+dx-1, idx[1]+dy-1, idx[2]+dz-1)
                if all(0<=ni[k]<dims[k] for k in range(3)): grid[ni] = True
        bx, ax = np.cumsum(grid,0)>0, np.cumsum(grid[::-1],0)[::-1]>0
        by, ay = np.cumsum(grid,1)>0, np.cumsum(grid[:,::-1],1)[:,::-1]>0
        bz, az = np.cumsum(grid,2)>0, np.cumsum(grid[:,:,::-1],2)[:,:,::-1]>0
        void = int(np.sum(bx & ax & by & ay & bz & az & ~grid)) * (res**3)
        return _make_result("PASS" if void < 1000 else "FAIL", f"Void={round(void,1)}A3", "THRESH:1000A3", lid)

    @staticmethod
    def compute_structural_characterization(structure):
        try:
            residues = _extract_residues(structure)
            helix, sheet, total = 0, 0, 0
            for i in range(1, len(residues)-1):
                p, c, n = residues[i-1], residues[i], residues[i+1]
                if not (_are_sequential(p, c) and _are_sequential(c, n)): continue
                if "C" not in p or not all(a in c for a in ["N","CA","C"]) or "N" not in n: continue
                phi, psi = _dihedral(p["C"], c["N"], c["CA"], c["C"]), _dihedral(c["N"], c["CA"], c["C"], n["N"])
                total += 1
                if -120 < phi < -30 and -70 < psi < 20: helix += 1
                elif -180 < phi < -40 and 60 < psi < 180: sheet += 1
            h_pct = round((helix/total)*100,1) if total else 0
            s_pct = round((sheet/total)*100,1) if total else 0
            return {"helix": h_pct, "sheet": s_pct, "loop": round(100-h_pct-s_pct,1), "total_atoms": len(structure.atoms), "total_residues": len(residues), "method": "phi_psi_rect"}
        except: return {"helix":0, "sheet":0, "loop":100, "total_atoms":0, "total_residues":0, "method":"failed"}

    @staticmethod
    def detect_confidence_source(structure):
        m = float(getattr(structure.confidence, 'mean_plddt', 0.0))
        return ("pLDDT", m, "explicit") if m > 0 else ("B-factor", 0.0, "missing")

    @staticmethod
    def decompose_confidence(structure):
        try:
            v = [float(a.b_iso) for a in structure.atoms]
            return {"mean": round(np.mean(v),1), "data_available": True} if v else {"mean":0, "data_available":False}
        except: return {"mean":0, "data_available":False}
