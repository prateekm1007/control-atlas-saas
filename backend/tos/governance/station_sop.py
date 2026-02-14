import hashlib as _hashlib

STATION_METADATA = {
    "version": "22.4.2",
    "institution": "TOSCANINI Structural Governance",
    "posture": "CONSERVATIVE",
    "scope": "GENERATIVE_AGNOSTIC",
    "format_fidelity": "High-Fidelity PDF v3",
    "min_coverage_pass": 70.0,
    "fringe_escalation_threshold": 20.0
}

# 1. THE RAW CANON (15 Laws)
LAW_CANON = {
    "LAW-100": {"title": "Bond Integrity", "principle": "Mean deviation < 0.05A.", "threshold": 0.05, "unit": "A", "severity": "FATAL"},
    "LAW-105": {"title": "Reliability Coverage", "principle": "Core >= 70% floor required.", "threshold": 70.0, "unit": "%", "severity": "FATAL"},
    "LAW-110": {"title": "Backbone Gap", "principle": "Connectivity < 1.5A.", "threshold": 1.5, "unit": "A", "severity": "FATAL"},
    "LAW-120": {"title": "Bond Angle", "principle": "3-point angles within 4.6 deg.", "threshold": 4.6, "unit": "deg", "severity": "FATAL"},
    "LAW-125": {"title": "Ramachandran", "principle": "Phi/Psi outliers < 5% in core.", "threshold": 5.0, "unit": "%", "severity": "FATAL"},
    "LAW-130": {"title": "Steric Clash", "principle": "Rejects overlaps > 5.0 score.", "threshold": 5.0, "unit": "score", "severity": "FATAL"},
    "LAW-135": {"title": "Omega Planarity", "principle": "Peptide planarity < 3% in core.", "threshold": 3.0, "unit": "%", "severity": "FATAL"},
    "LAW-145": {"title": "Chirality", "principle": "L-amino acid consistency.", "threshold": 0.0, "unit": "count", "severity": "FATAL"},
    "LAW-150": {"title": "Rotamer Audit", "principle": "Sidechain chi1 consistency.", "threshold": 20.0, "unit": "%", "severity": "ADVISORY"},
    "LAW-160": {"title": "Chain Integrity", "principle": "C-alpha trace continuity.", "threshold": 4.2, "unit": "A", "severity": "FATAL"},
    "LAW-170": {"title": "Residue Identity", "principle": "Standard AA mapping.", "threshold": 0, "unit": "unk", "severity": "ADVISORY"},
    "LAW-182": {"title": "Hydrophobic Burial", "principle": "Burial ratio > 0.3.", "threshold": 0.3, "unit": "ratio", "severity": "ADVISORY"},
    "LAW-195": {"title": "Disulfide Geometry", "principle": "S-S bond length (2.04A).", "threshold": 0.20, "unit": "A", "severity": "FATAL"},
    "LAW-155": {"title": "Voxel Occupancy", "principle": "Atomic connectivity.", "threshold": 2.0, "unit": "V", "severity": "ADVISORY"},
    "LAW-200": {"title": "Internal Cavity", "principle": "Voids < 1000 A^3.", "threshold": 1000, "unit": "A3", "severity": "ADVISORY"}
}

# 2. METHOD CLASSIFICATIONS
LAW_METHOD_CLASSIFICATIONS = {lid: "deterministic" if lid not in ["LAW-155", "LAW-182", "LAW-200"] else "heuristic" for lid in LAW_CANON}

# 3. COMPUTED LISTS & COUNTS (Requirement for main.py)
DETERMINISTIC_LAWS = sorted([l for l, m in LAW_METHOD_CLASSIFICATIONS.items() if m == "deterministic"])
HEURISTIC_LAWS = sorted([l for l, m in LAW_METHOD_CLASSIFICATIONS.items() if m == "heuristic"])
TOTAL_LAWS = len(LAW_CANON)
DETERMINISTIC_COUNT = len(DETERMINISTIC_LAWS)
HEURISTIC_COUNT = len(HEURISTIC_LAWS)

# 4. CATEGORIES (Requirement for pdf_generator.py)
LAW_CATEGORIES = {lid: "geometric" if lid not in ["LAW-145", "LAW-150", "LAW-170", "LAW-182", "LAW-195"] else "chemical" for lid in LAW_CANON}
CATEGORY_DISPLAY = {"geometric": "GEOMETRIC INVARIANTS", "chemical": "CHEMICAL INVARIANTS"}

def get_laws_by_category():
    g = {}
    for l, c in LAW_CATEGORIES.items(): g.setdefault(c, []).append(l)
    return g

# 5. SCORE DEFINITIONS (Requirement for dashboard and pdf)
SCORE_DEFINITIONS = {
    "DETERMINISTIC_SCORE": {"title": "Deterministic Integrity", "explanation": "VETO gate. Requires 100% core integrity."},
    "ADVISORY_SCORE": {"title": "Heuristic Advisory", "explanation": "Informational metrics."},
    "ML_CONFIDENCE": {"title": "ML Confidence", "explanation": "Mean pLDDT."},
    "STRATEGIC_SCORE": {"title": "EPI Priority Index", "explanation": "Prioritization metric."}
}

# 6. PIPELINE CONSTANTS
BAYESIAN_FORMULA = "P = clamp(S6, 0, 1) x W_arch x (1 - M_S8)"
ARCHITECTURE_WEIGHTS = {"NONE": 1.0, "LINEAR": 0.95, "MULTIVALENT": 0.85, "ENGAGEMENT": 0.80, "MOTIFS": 0.90, "LINKERS": 0.88, "CONFORMATIONAL": 0.75, "METAL": 0.70}
STANDARD_RESIDUES = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
HYDROPHOBIC_RESIDUES = {"ALA","VAL","ILE","LEU","MET","PHE","TRP","PRO"}
CONFIDENCE_EXCLUSION_POLICY = {"trigger": "mean_confidence <= 0", "action": "exclude"}

# 7. PHYSICS CONSTANTS (Requirement for tier1_measurements.py)
SIGMA_TABLE = {"N-CA": 0.010, "CA-C": 0.021, "C-O": 0.020, "C-N": 0.014, "CA-CB": 0.020, "C-S": 0.025, "S-S": 0.016, "N-CA-C": 2.5, "CA-C-N": 1.5, "O-C-N": 1.6}
IDEAL_TABLE = {"N-CA": 1.470, "CA-C": 1.525, "C-O": 1.231, "C-N": 1.329, "CA-CB": 1.530, "C-S": 1.803, "S-S": 2.033, "N-CA-C": 111.0, "CA-C-N": 116.6, "O-C-N": 123.0}

def compute_canon_hash():
    s = "".join([f"{k}|{LAW_METHOD_CLASSIFICATIONS[k]}" for k in sorted(LAW_CANON.keys())])
    return _hashlib.sha256(s.encode()).hexdigest()[:16]
LAW_CANON_HASH = compute_canon_hash()

