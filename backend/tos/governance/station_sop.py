import hashlib as _hashlib

STATION_METADATA = {
    "version": "22.5.3",
    "institution": "TOSCANINI Structural Governance",
    "posture": "CONSERVATIVE",
    "min_coverage_pass": 70.0,
    "ml_sigma_multiplier": 1.25,
    "schema_version": "1.0.0"
}

# ═══════════════════════════════════════════════════════════════════════════════
# PIL-CHM-03: Chemical Canon (Unit Standardization)
#
# CALIBRATION NOTES (v22.5.2):
#
# LAW-100 (Bond Integrity):
#   Previous: 4.0σ — too strict for macromolecular crystallography.
#   At 1.7Å resolution, coordinate uncertainty is ~0.1-0.2Å, which translates
#   to z-scores of 5-15σ against Engh-Huber ideals (σ = 0.01-0.02Å).
#   MolProbity flags bonds at 4σ as *noteworthy* but does not reject structures
#   unless >5% of bonds exceed 4σ. We now measure the PERCENTAGE of outlier
#   bonds (z > 4σ) against total bonds, with threshold at 5%.
#   Reference: Engh & Huber (1991), Acta Cryst. A47, 392-400.
#
# LAW-110 (Backbone Gap):
#   Previous: 1.5Å — catches normal peptide bonds at medium resolution.
#   Inter-residue C(i)-N(i+1) peptide bond: ideal 1.329Å, but coordinate
#   uncertainty at 1.5-2.5Å resolution means observed values of 1.3-1.7Å
#   are routine. True backbone gaps (missing residues) show C-N > 3.0Å.
#   Recalibrated to 2.0Å: catches real gaps, ignores coordinate imprecision.
#   Reference: Word et al. (1999), J. Mol. Biol. 285, 1735-1747.
#
# LAW-160 (Chain Integrity):
#   Previous: 4.2Å — too tight for normal Cα spacing variation.
#   Trans-peptide Cα(i)-Cα(i+1): ideal 3.8Å, observed range 3.6-4.3Å
#   in well-refined structures. Cis-proline: ~2.9Å. True chain breaks
#   show CA-CA > 5.0Å. Recalibrated to 4.5Å.
#   Reference: Ramachandran & Sasisekharan (1968), Adv. Protein Chem. 23.
#
# LAW-182 (Hydrophobic Burial):
#   Previous: 0.3 ratio with 12Å radius — too dependent on protein size.
#   Multi-chain complexes (e.g., 4HHB tetramer) spread the centroid,
#   making legitimate hydrophobic cores appear unburied. Adjusted to
#   per-chain evaluation with 15Å radius.
#   NOTE: This remains heuristic (LAW-182 is classified heuristic).
# ═══════════════════════════════════════════════════════════════════════════════

LAW_CANON = {
    "LAW-100": {"title": "Bond Integrity", "threshold": 5.0, "unit": "%", "operator": "<=", "scope": "core residues", "type": "Percentage"},
    "LAW-105": {"title": "Reliability Coverage", "threshold": 70.0, "unit": "%", "operator": ">=", "scope": "total residues", "type": "Percentage"},
    "LAW-110": {"title": "Backbone Gap", "threshold": 2.0, "unit": "Å", "operator": "<=", "scope": "sequential pairs", "type": "Scalar"},
    "LAW-120": {"title": "Bond Angle", "threshold": 18.0, "unit": "°", "operator": "<=", "scope": "core residues", "type": "RMSD"},
    "LAW-125": {"title": "Ramachandran", "threshold": 5.0, "unit": "%", "operator": "<=", "scope": "core residues", "type": "Percentage"},
    "LAW-130": {"title": "Clashscore", "threshold": 20.0, "unit": "clashes/1000 atoms", "operator": "<=", "scope": "total atoms", "type": "Rate"},
    "LAW-135": {"title": "Omega Planarity", "threshold": 3.0, "unit": "%", "operator": "<=", "scope": "core residues", "type": "Percentage"},
    "LAW-145": {"title": "Chirality", "threshold": 0.0, "unit": "count", "operator": "=", "scope": "core residues", "type": "Count"},
    "LAW-150": {"title": "Rotamer Audit", "threshold": 20.0, "unit": "%", "operator": "<=", "scope": "core residues", "type": "Percentage"},
    "LAW-160": {"title": "Chain Integrity", "threshold": 4.5, "unit": "Å", "operator": "<=", "scope": "sequential C-alphas", "type": "Scalar"},
    "LAW-170": {"title": "Residue Identity", "threshold": 0.0, "unit": "count", "operator": "=", "scope": "total residues", "type": "Count"},
    "LAW-182": {"title": "Hydrophobic Burial", "threshold": 0.3, "unit": "ratio", "operator": ">=", "scope": "hydrophobic core", "type": "Rate"},
    "LAW-195": {"title": "Disulfide Geometry", "threshold": 0.20, "unit": "Å", "operator": "<=", "scope": "cys pairs", "type": "Scalar"},
    "LAW-155": {"title": "Voxel Occupancy", "threshold": 2.0, "unit": "V", "operator": ">=", "scope": "atomic grid", "type": "Scalar"},
    "LAW-200": {"title": "Packing Quality", "threshold": 300.0, "unit": "Å³/atom", "operator": "<=", "scope": "bounding box", "type": "Rate"}
}

LAW_METHOD_CLASSIFICATIONS = {lid: "deterministic" if lid not in ["LAW-155", "LAW-182", "LAW-200"] else "heuristic" for lid in LAW_CANON}
DETERMINISTIC_LAWS = sorted([l for l, m in LAW_METHOD_CLASSIFICATIONS.items() if m == "deterministic"])
HEURISTIC_LAWS = sorted([l for l, m in LAW_METHOD_CLASSIFICATIONS.items() if m == "heuristic"])
TOTAL_LAWS, DETERMINISTIC_COUNT, HEURISTIC_COUNT = len(LAW_CANON), len(DETERMINISTIC_LAWS), len(HEURISTIC_LAWS)

LAW_CATEGORIES = {lid: "geometric" if lid not in ["LAW-145", "LAW-150", "LAW-170", "LAW-182", "LAW-195"] else "chemical" for lid in LAW_CANON}
CATEGORY_DISPLAY = {"geometric": "GEOMETRIC INVARIANTS", "chemical": "CHEMICAL INVARIANTS"}

def get_laws_by_category():
    g = {}
    for l, c in LAW_CATEGORIES.items(): g.setdefault(c, []).append(l)
    return {k: sorted(v) for k, v in g.items()}

SCORE_DEFINITIONS = {
    "DETERMINISTIC_SCORE": {"title": "Deterministic Integrity", "explanation": "Physical compliance relative to crystallographic ideals."},
    "ADVISORY_SCORE": {"title": "Heuristic Advisory", "explanation": "Statistical proxies for structural plausibility."},
    "ML_CONFIDENCE": {"title": "ML Confidence", "explanation": "Model-reported mean pLDDT."},
    "STRATEGIC_SCORE": {"title": "EPI Priority Index", "explanation": "Non-probabilistic prioritization metric for lead selection."}
}

BAYESIAN_FORMULA = "P = clamp(S6, 0, 1) x W_arch x (1 - M_S8)"

ARCHITECTURE_WEIGHTS = {
    "NONE": 1.0, "LINEAR": 0.95, "MULTIVALENT": 0.85, "ENGAGEMENT": 0.80,
    "MOTIFS": 0.90, "LINKERS": 0.88, "CONFORMATIONAL": 0.75, "METAL": 0.70
}

SIGMA_TABLE = {"N-CA": 0.010, "CA-C": 0.021, "C-O": 0.020, "C-N": 0.014, "CA-CB": 0.020, "C-S": 0.025, "S-S": 0.016, "N-CA-C": 2.5, "CA-C-N": 1.5, "O-C-N": 1.6}
IDEAL_TABLE = {"N-CA": 1.470, "CA-C": 1.525, "C-O": 1.231, "C-N": 1.329, "CA-CB": 1.530, "C-S": 1.803, "S-S": 2.033, "N-CA-C": 111.0, "CA-C-N": 116.6, "O-C-N": 123.0}

STANDARD_RESIDUES = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
HYDROPHOBIC_RESIDUES = {"ALA","VAL","ILE","LEU","MET","PHE","TRP","PRO"}

def compute_canon_hash():
    s = "".join([f"{k}|{LAW_METHOD_CLASSIFICATIONS[k]}" for k in sorted(LAW_CANON.keys())])
    return _hashlib.sha256(s.encode()).hexdigest()[:16]

LAW_CANON_HASH = compute_canon_hash()
