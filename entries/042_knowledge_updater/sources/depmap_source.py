"""
Checks for DepMap (Essentiality) updates.
Mocked for v1 - enables architecture without huge download.
"""

# Simulating an API check
MOCK_UPDATES = {
    "KRAS": {"score": 1.0, "version": "23Q4"},
    "TP53": {"score": 0.2, "version": "23Q4"},
    # Simulate a new update for EGFR
    "EGFR": {"score": 0.95, "version": "24Q1"} 
}

def check_update(gene_name, current_version):
    """Check if a newer DepMap version exists."""
    data = MOCK_UPDATES.get(gene_name.upper())
    if not data:
        return None
        
    if data["version"] > current_version:
        return data
    return None
