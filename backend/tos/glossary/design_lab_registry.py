ARCHITECTURES = {
    "LINEAR": {
        "success_criteria": "Requires high helicity propensity and defined pocket depth.",
        "failure_modes": ["Aggregation via exposed hydrophobic patches", "Proteolytic instability"],
        "trade_off": "Fast kinetics but lower affinity compared to Multivalent."
    },
    "MULTIVALENT": {
        "success_criteria": "Requires symmetric epitopes and precise linker length matching.",
        "failure_modes": ["Entropy penalty from excessive linker length", "Self-association"],
        "trade_off": "1000x Avidity gain vs significant synthesis complexity."
    },
    "ENGAGEMENT": {
        "success_criteria": "Requires high shape complementarity (Sc) and surface match.",
        "failure_modes": ["Induced-fit energy penalty", "Lack of specific polar registry"],
        "trade_off": "High specificity but sensitive to conformational noise."
    },
    "MOTIFS": {
        "success_criteria": "Requires validated pharmacophore alignment.",
        "failure_modes": ["Side-chain rotamer clashing", "Electronic misalignment"],
        "trade_off": "Potent binding but narrow search space."
    },
    "LINKERS": {
        "success_criteria": "Adherence to LAW-150 (<3x EAAAK).",
        "failure_modes": ["Lever-paradox collapse", "Torsional entropy bleed"],
        "trade_off": "High modularity vs risk of folding interference."
    },
    "CONFORMATIONAL": {
        "success_criteria": "N/C-terminal capping and salt-bridge registry.",
        "failure_modes": ["Global fold destabilization", "Disulfide mismatch"],
        "trade_off": "Extreme stability but difficult to engineer dynamically."
    },
    "METAL": {
        "success_criteria": "Tetrahedral or Octahedral coordination registry.",
        "failure_modes": ["Ion dissociation", "Misfolding in low-metal environments"],
        "trade_off": "Absolute topological rigidity vs high environmental sensitivity."
    }
}
