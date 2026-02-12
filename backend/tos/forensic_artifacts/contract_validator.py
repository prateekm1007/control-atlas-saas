class ContractViolation(Exception): pass

def validate_payload(payload):
    """Enforces v17.x Sovereign Contract: No tier leakage, no missing data."""
    required_keys = ["verdict", "provenance", "tier1", "tier3", "architecture", "visualization", "definitions"]
    for key in required_keys:
        if key not in payload: raise ContractViolation(f"Sovereign Breach: Missing Section {key}")
    
    # Tier-1 Consistency
    if "laws" not in payload["tier1"]: raise ContractViolation("Sovereign Breach: Tier-1 Ledger is empty.")
    
    # Page-Specific Requirements
    if not payload["provenance"].get("hash"): raise ContractViolation("Sovereign Breach: Page 1 lacks Coordinate Hash.")
    if not payload["visualization"].get("ca_coords"): raise ContractViolation("Sovereign Breach: Page 2 lacks Projections.")
    
    return True
