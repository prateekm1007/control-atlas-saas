"""
TOSCANINI STATION SOP // INVISIBLE GOVERNANCE FILE
Strict adherence to these constants is mandatory for all modules.
"""

STATION_METADATA = {
    "version": "21.36.0",
    "institution": "TOSCANINI (TM) Structural Governance System",
    "confidentiality": "Confidential – Institutional Use Only",
    "format_fidelity": "LaTeX-Standard-Dossier-v1",
    "upload_limit_mb": 25,
    "allowed_extensions": [".pdb", ".cif", ".mmcif"]
}

LAW_CANON = {
    "LAW-100": {"title": "Bond Length Conformity", "principle": "Check bond lengths deviate < 0.05A from Engh-Huber standards."},
    "LAW-110": {"title": "Backbone Continuity", "principle": "Verification of peptide bond connectivity. Rejects gaps > 1.5A."},
    "LAW-120": {"title": "Bond Angle Conformity", "principle": "Validation of 3-point atomic angles within 4.6 degrees of ideal."},
    "LAW-130": {"title": "Steric Clash Exclusion", "principle": "Voxel-based proximity check. Rejects overlaps > 0.4A."},
    "LAW-140": {"title": "Residue Packing Density", "principle": "Evaluates core compactness via hydrophobic collapse norms."},
    "LAW-155": {"title": "Voxel Occupancy Continuity", "principle": "Grid-based verification of atomic cluster connectivity."},
    "LAW-160": {"title": "Chain Integrity", "principle": "Topological validation of C-alpha trace to ensure single fold."},
    "LAW-170": {"title": "Residue Identity Validation", "principle": "Verifies every residue maps to a valid amino acid stereoisomer."},
    "LAW-182": {"title": "Hydrophobic Burial Ratio", "principle": "Calculates ratio of non-polar surface area burial."},
    "LAW-195": {"title": "Disulfide Geometry Validation", "principle": "Checks S-S bond lengths (2.04A) and orientations."},
    "LAW-200": {"title": "Internal Cavity Assessment", "principle": "Spatial audit for significant internal voids (>200 A^3)."}
}

SCORE_DEFINITIONS = {
    "PHYSICAL_SCORE": {
        "title": "Physical Integrity",
        "explanation": "Adherence to Tier-1 physical invariants governing geometry and feasibility.",
        "impact": "If <100%, the structure is automatically vetoed."
    },
    "ML_CONFIDENCE": {
        "title": "ML Confidence (pLDDT)",
        "explanation": "Predicted Local Distance Difference Test score (0-100).",
        "impact": "Confidence cannot override physical law violations."
    },
    "STRATEGIC_SCORE": {
        "title": "EPI Priority Index",
        "explanation": "Bayesian success probability math used for wet-lab allocation.",
        "impact": "Ranks candidates for resource allocation."
    }
}

BAYESIAN_FORMULA = "P_final = (P_base + P_phys + P_ml - Tax) x W_deriv x (1 - M_S8)"
