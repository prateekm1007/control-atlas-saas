CANON = {
    "LAW-120": {"title": "Peptide Bond Sanity", "principle": "Covalent bond length (1.33A).", "rationale": "Peptide bonds are resonance hybrids. Violations indicate quantum impossibility."},
    "LAW-155": {"title": "Steric Clash Prohibition", "principle": "Atoms cannot overlap.", "rationale": "The Pauli exclusion principle forbids electron cloud overlap. Non-bonded atoms < 2.0A create physically impossible geometry."},
    "LAW-160": {"title": "Backbone Continuity", "principle": "Consistent C-alpha spacing (~3.8A).", "rationale": "Deviations > 4.5A indicate a torn or hallucinated chain."},
    "LAW-182": {"title": "Hydrophobic Burial", "principle": "Non-polar residues must be buried.", "rationale": "Surface hydrophobic patches lead to immediate aggregation."},
    "LAW-190": {"title": "Ring Planarity", "principle": "Aromatic rings must remain flat.", "rationale": "Delocalized pi-bonding enforces planarity."},
    "LAW-195": {"title": "Disulfide Geometry", "principle": "S-S distance approx 2.05A.", "rationale": "Incorrect geometry traps designs in non-functional states."},
    "LAW-200": {"title": "Cavity Collapse", "principle": "Functional pockets must exist.", "rationale": "Pocket collapse renders the design dead."},
    "LAW-210": {"title": "Helix Crossing", "principle": "Canonical packing angles.", "rationale": "Ridges-into-grooves packing geometry constraints."},
    "LAW-220": {"title": "Beta Strand Registry", "principle": "Precise H-bond alignment.", "rationale": "Beta-sheet stability requires spatial registry."},
    "LAW-230": {"title": "Atomic Valence", "principle": "Quantum valid coordination.", "rationale": "Atoms have fixed coordination numbers."}
}
def list_all_law_ids(): return list(CANON.keys())
def get_law_explanation(lid): return CANON.get(lid, {"title": lid, "principle": "N/A", "rationale": "N/A"})
