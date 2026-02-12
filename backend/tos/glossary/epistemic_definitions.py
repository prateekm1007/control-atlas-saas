DEFINITIONS = {
    "PHYSICAL_SCORE": {
        "title": "Physical Integrity",
        "explanation": "Deterministic adherence to Tier-1 physical invariants governing bond geometry, sterics, and topology. 100% indicates zero structural impossibilities detected.",
        "impact": "Structures <100% are physically impossible and are automatically vetoed."
    },
    "ML_CONFIDENCE": {
        "title": "ML Confidence (pLDDT)",
        "explanation": "Predicted Local Distance Difference Test (0-100). Measures the model's internal certainty of the local fold topology.",
        "impact": "Confidence cannot override physical law violations."
    },
    "STRATEGIC_SCORE": {
        "title": "EPI Priority Index",
        "explanation": "Expected Prioritization Index. Bayesian success probability: P = (P_base + P_phys + P_ml - Tax) x (1 - M_S8).",
        "impact": "Used for ranking candidates for wet-lab resource allocation."
    },
    "PILLAR_MANIFEST": [
        ("PIL-VIS-01: Dual-Engine Visualizer", "Auto-switch styles (Cartoon/PASS vs Stick/VETO)."),
        ("PIL-LED-02: Diagnostic Ledger", "Law-by-law evidence expansion with measurements."),
        ("PIL-CHM-03: Chemical Canon", "11 Tier-1 Invariants governed by a context-aware registry."),
        ("PIL-DSC-04: Discovery Console", "Exhaustive resolution of UniProt and AFDB identifiers."),
        ("PIL-ACQ-05: Industrial Acquisition", "Multi-format fetch loop with mmCIF and PDB fallback."),
        ("PIL-ESM-06: ESM Atlas Integration", "Integration with zero-latency transformer-based folding."),
        ("PIL-EPI-07: Epistemic Duality", "Explicit separation of Physical Fact from AI Estimation."),
        ("PIL-VOX-08: Voxelized Physics", "O(N) spatial-grid reasoning for multimer clash audits."),
        ("PIL-BBK-09: Backbone Continuity", "LAW-160 monitor for generative backbone tears."),
        ("PIL-INS-10: Forensic Inspection", "Zoom-to-coordinate failure marker autopsy."),
        ("PIL-NAR-11: PhD Narrative Engine", "16-tier Gemini hierarchy failure recovery and expert synthesis."),
        ("PIL-MTH-12: Strategic Success Math", "Bayesian success probability math (S6 x S8)."),
        ("PIL-ARC-13: Architecture Lab", "9 geometry-derived classification categories (Warheads)."),
        ("PIL-PAR-14: Unified Parity", "UI and PDF scores must be bit-equivalent."),
        ("PIL-AUT-15: Zero-Button Automation", "One-click 'Fetch -> Audit -> Seal' workflow."),
        ("PIL-NKG-16: Negative Knowledge Graph", "Persistent JSONL memory of physical failure modes."),
        ("PIL-SYM-17: Visual Symmetry", "Symmetrical, centered dual-view topological projections."),
        ("PIL-NOT-18: Notary Seal", "SHA-256 Hashing of the coordinate artifact."),
        ("PIL-LCY-19: Program Lifecycle", "Forward-only state machine (S0–S7)."),
        ("PIL-IMM-20: Immutable Ingestion", "Data frozen in frozen=True structural dataclasses."),
        ("PIL-HYG-21: Institutional Hygiene", "Type-safe guards and ASCII-only sanitization.")
    ]
}
