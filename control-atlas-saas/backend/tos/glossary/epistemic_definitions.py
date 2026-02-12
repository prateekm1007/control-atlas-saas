DEFINITIONS = {
    "PHYSICAL_SCORE": {
        "title": "Physical Integrity Score (Invariant)",
        "explanation": "Derived from deterministic Euclidean geometry. It measures if the atoms 'fit' in 3D space without violating the Pauli Exclusion Principle (steric clashes) or covalent bond geometry. A 100% score means the structure is chemically permissible."
    },
    "CONFIDENCE_SCORE": {
        "title": "ML Confidence (pLDDT)",
        "explanation": "Calculated via the 'Predicted Local Distance Difference Test' (pLDDT). The model (AlphaFold/ESMFold) uses a dedicated head to predict the per-residue lDDT-Cα score. A score of 54.8% indicates the model predicts the backbone is likely correct, but side-chain orientations and local tertiary 'packing' are speculative. For Drug Discovery, scores <70 in binding pockets render docking simulations unreliable."
    },
    "ROUTING_BANNER": {
        "title": "Forensic Routing",
        "explanation": "A 'Physically Valid but Low Confidence' structure is a high-risk synthesis candidate. It obeys physics, but lacks the evidential support required for high-stakes Pharma spend."
    }
}
