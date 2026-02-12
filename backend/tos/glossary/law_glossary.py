def list_all_law_ids():
    return ["LAW-100", "LAW-110", "LAW-120", "LAW-130", "LAW-140", "LAW-155", "LAW-160", "LAW-170", "LAW-182", "LAW-195", "LAW-200"]

def get_law_explanation(law_id):
    explanations = {
        "LAW-100": {"title": "Bond Length Conformity", "principle": "Deterministic check ensuring C-C, C-N, and C-O bond lengths deviate less than 0.05A from Engh-Huber standard values."},
        "LAW-110": {"title": "Backbone Continuity", "principle": "Verification of peptide bond connectivity. Rejects structures with unphysical C-N gaps (>1.5A) or generative 'tears'."},
        "LAW-120": {"title": "Bond Angle Conformity", "principle": "Validation of 3-point atomic angles. Ensures geometry is within 4.6 degrees of ideal hybridization states."},
        "LAW-130": {"title": "Steric Clash Exclusion", "principle": "Voxel-based proximity check. Rejects any non-bonding atom pairs that overlap by more than 0.4A of their van der Waals radii."},
        "LAW-140": {"title": "Residue Packing Density", "principle": "Evaluates the internal core compactness. Rejects models with excessive 'hollow' space that violates hydrophobic collapse norms."},
        "LAW-155": {"title": "Voxel Occupancy Continuity", "principle": "Grid-based density verification. Rejects structures with disconnected atomic clusters or 'ghost' atoms."},
        "LAW-160": {"title": "Chain Integrity", "principle": "Topological validation of C-alpha trace. Ensures the program is analyzing a single, valid fold without fragments."},
        "LAW-170": {"title": "Residue Identity Validation", "principle": "Verifies every residue maps to a valid 20-amino acid stereoisomer. Rejects unknown or malformed residue types."},
        "LAW-182": {"title": "Hydrophobic Burial Ratio", "principle": "Calculates the ratio of non-polar surface area. Rejects structures where hydrophobic sidechains are over-exposed to solvent."},
        "LAW-195": {"title": "Disulfide Geometry Validation", "principle": "Checks S-S bond lengths (2.04A) and dihedral angles. Rejects non-physical disulfide orientations."},
        "LAW-200": {"title": "Internal Cavity Assessment", "principle": "Spatial audit for significant internal voids (>200 cubic Angstroms) that would render the structure unstable."}
    }
    return explanations.get(law_id, {"title": "Physical Invariant", "principle": "Standard Chemical Canon Law."})
