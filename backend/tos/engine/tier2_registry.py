TIER2_CANON = {
    "T2-301": {
        "name": "Packing Efficiency",
        "why_tier2": "Structural quality metric. Poor packing doesn't violate physics, but predicts low thermodynamic stability.",
        "rationale": "Measures internal atomic voids. High packing is a hallmark of well-designed cores.",
        "interpretation": "High optimization potential if >= 90%."
    },
    "T2-302": {
        "name": "Solvent Accessibility",
        "why_tier2": "Manufacturing risk metric. Exposed hydrophobic patches are physically possible but cause aggregation.",
        "rationale": "Evaluates the burial of non-polar residues to prevent solution-phase collapse.",
        "interpretation": "Usable if >= 70%."
    }
}
