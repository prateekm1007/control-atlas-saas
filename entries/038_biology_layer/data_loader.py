"""
Entry 038 — Biology Layer Data Loader
Simulates access to DepMap (Essentiality) and GTEx (Expression).
"""

# Mock Knowledge Graph for Oncology Targets
BIOLOGY_KNOWLEDGE = {
    "KRAS": {
        "essentiality_score": 1.0, # Essential in KRAS-mutant lines
        "expression_lung": "HIGH",
        "cancer_driver": True
    },
    "TP53": {
        "essentiality_score": 0.2, # Loss of function driver (restoration needed)
        "expression_lung": "MEDIUM",
        "cancer_driver": True
    },
    "EGFR": {
        "essentiality_score": 0.9, # Essential in EGFR-mutant lines
        "expression_lung": "HIGH",
        "cancer_driver": True
    },
    "GAPDH": {
        "essentiality_score": 0.0, # Housekeeping/Non-driver context
        "expression_lung": "HIGH",
        "cancer_driver": False
    }
}

def get_target_data(gene_name: str) -> dict:
    """Retrieve biological context for a gene."""
    return BIOLOGY_KNOWLEDGE.get(gene_name.upper(), {})
