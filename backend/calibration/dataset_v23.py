"""
TOSCANINI Phase 2 Calibration Dataset (v23)

50 structures selected for statistical coverage across:
- Experimental methods: X-ray, Cryo-EM, NMR
- Resolution ranges: <1.5A, 1.5-2.5A, 2.5-3.5A, >3.5A
- Structural diversity: single chain, multi-chain, membrane, enzyme
- Known edge cases: disulfides, cis-prolines, insertion codes

Selection criteria:
- All structures deposited in RCSB PDB or AlphaFold DB
- All have been validated by wwPDB validation pipeline
- Resolution documented where applicable
- No structures from the 7 benchmark set (avoid circular calibration)
"""

CALIBRATION_DATASET = {
    "xray_ultrahigh": [
        {"pdb_id": "2VB1", "resolution": 0.48, "desc": "Crambin ultra-high res"},
        {"pdb_id": "1EJG", "resolution": 0.54, "desc": "Cholesterol oxidase"},
        {"pdb_id": "2OV0", "resolution": 0.75, "desc": "Lysozyme ultra-high"},
        {"pdb_id": "1I0T", "resolution": 0.89, "desc": "Ferredoxin"},
        {"pdb_id": "3NIR", "resolution": 1.00, "desc": "Concanavalin A"},
        {"pdb_id": "1GCI", "resolution": 1.10, "desc": "Glucose isomerase"},
        {"pdb_id": "2ACE", "resolution": 1.10, "desc": "Acetylcholinesterase"},
        {"pdb_id": "1BPI", "resolution": 1.10, "desc": "BPTI"},
    ],
    "xray_medium": [
        {"pdb_id": "1TEN", "resolution": 1.60, "desc": "Tenascin"},
        {"pdb_id": "3LZM", "resolution": 1.70, "desc": "T4 lysozyme"},
        {"pdb_id": "1AKE", "resolution": 2.00, "desc": "Adenylate kinase"},
        {"pdb_id": "2CBA", "resolution": 2.00, "desc": "Carbonic anhydrase"},
        {"pdb_id": "1HHO", "resolution": 2.10, "desc": "Deoxyhemoglobin"},
        {"pdb_id": "3PGK", "resolution": 2.10, "desc": "Phosphoglycerate kinase"},
        {"pdb_id": "1GFL", "resolution": 1.90, "desc": "GFP"},
        {"pdb_id": "1MBO", "resolution": 2.00, "desc": "Myoglobin"},
        {"pdb_id": "2LYZ", "resolution": 2.00, "desc": "Hen lysozyme"},
        {"pdb_id": "1PPE", "resolution": 1.80, "desc": "Porcine elastase"},
    ],
    "xray_low": [
        {"pdb_id": "1AON", "resolution": 3.00, "desc": "GroEL-GroES complex"},
        {"pdb_id": "1FIN", "resolution": 2.80, "desc": "CDK2-cyclin A"},
        {"pdb_id": "1DDL", "resolution": 3.50, "desc": "Dihydroorotase"},
        {"pdb_id": "1NEK", "resolution": 2.80, "desc": "Nek2 kinase"},
        {"pdb_id": "1BRS", "resolution": 2.80, "desc": "Barnase-barstar"},
        {"pdb_id": "2HBS", "resolution": 2.50, "desc": "Sickle cell hemoglobin"},
    ],
    "cryo_em": [
        {"pdb_id": "5A63", "resolution": 2.80, "desc": "TRPV1 channel"},
        {"pdb_id": "6J6J", "resolution": 2.20, "desc": "Human 26S proteasome"},
        {"pdb_id": "6WLC", "resolution": 2.90, "desc": "SARS-CoV-2 nucleocapsid"},
        {"pdb_id": "7K3G", "resolution": 2.60, "desc": "SARS-CoV-2 spike closed"},
        {"pdb_id": "3J9I", "resolution": 3.60, "desc": "80S ribosome"},
        {"pdb_id": "5T4O", "resolution": 3.40, "desc": "Spliceosome"},
        {"pdb_id": "6NB6", "resolution": 3.00, "desc": "ABC transporter"},
        {"pdb_id": "7JTL", "resolution": 2.30, "desc": "SARS-CoV-2 RBD-ACE2"},
    ],
    "nmr": [
        {"pdb_id": "1D3Z", "resolution": None, "desc": "Ubiquitin NMR"},
        {"pdb_id": "2KOX", "resolution": None, "desc": "GB1 domain NMR"},
        {"pdb_id": "1GB1", "resolution": None, "desc": "Protein G B1 NMR"},
        {"pdb_id": "2JOF", "resolution": None, "desc": "Trp-cage miniprotein NMR"},
        {"pdb_id": "1L2Y", "resolution": None, "desc": "Villin headpiece NMR"},
        {"pdb_id": "2N0A", "resolution": None, "desc": "Alpha-synuclein NMR"},
        {"pdb_id": "1RVS", "resolution": None, "desc": "BPTI NMR"},
        {"pdb_id": "2MJB", "resolution": None, "desc": "Amyloid beta NMR"},
    ],
    "predicted": [
        {"pdb_id": "AF-P04637-F1", "resolution": None, "desc": "p53 tumor suppressor"},
        {"pdb_id": "AF-P68871-F1", "resolution": None, "desc": "Hemoglobin beta"},
        {"pdb_id": "AF-P0DTC2-F1", "resolution": None, "desc": "SARS-CoV-2 spike"},
        {"pdb_id": "AF-Q9Y6K9-F1", "resolution": None, "desc": "NF-kB essential modulator"},
        {"pdb_id": "AF-P21802-F1", "resolution": None, "desc": "FGFR2"},
        {"pdb_id": "AF-O15111-F1", "resolution": None, "desc": "IKK-alpha"},
        {"pdb_id": "AF-P62988-F1", "resolution": None, "desc": "Ubiquitin"},
        {"pdb_id": "AF-P02768-F1", "resolution": None, "desc": "Serum albumin"},
        {"pdb_id": "AF-P10636-F1", "resolution": None, "desc": "Tau protein"},
        {"pdb_id": "AF-P00533-F1", "resolution": None, "desc": "EGFR"},
    ],
}


def get_all_targets():
    """Return flat list of (pdb_id, category, resolution, description)."""
    targets = []
    for category, entries in CALIBRATION_DATASET.items():
        for entry in entries:
            targets.append((
                entry["pdb_id"],
                category,
                entry["resolution"],
                entry["desc"]
            ))
    return targets


if __name__ == "__main__":
    targets = get_all_targets()
    print(f"Total calibration targets: {len(targets)}")
    for cat, entries in CALIBRATION_DATASET.items():
        print(f"  {cat}: {len(entries)}")
