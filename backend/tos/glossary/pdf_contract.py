PDF_CONTRACT = {
    "version": "17.2.0",
    "pages": [
        {"index": 1, "name": "Notary Front Page", "required": ["verdict.binary", "verdict.physical_score", "provenance.audit_id"]},
        {"index": 2, "name": "Topological Projection", "required": ["visualization.ca_coords"]},
        {"index": 3, "name": "Physics Ledger", "required": ["tier1.laws"]},
        {"index": 4, "name": "Engineering Lenses", "required": ["tier2.laws"], "condition": "governance.engineering_enabled"},
        {"index": 5, "name": "Strategic Success Probability", "required": ["tier3.probability", "tier3.formula"]},
        {"index": 6, "name": "Architecture Decision & Rejects", "required": ["architecture.category", "architecture.deep_dive"]},
        {"index": 7, "name": "Empirical Provenance", "required": ["tier3.priors_used"]}
    ]
}
