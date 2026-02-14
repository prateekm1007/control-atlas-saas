import sys, os
sys.path.insert(0, './backend')
from tos.governance.station_sop import (
    LAW_CANON, DETERMINISTIC_LAWS, HEURISTIC_LAWS, DETERMINISTIC_COUNT
)

def run_contract_audit():
    print(f"Executing Sovereign Contract Audit...")
    
    # 1. Consistency Check: No Magic Numbers
    total_derived = len(DETERMINISTIC_LAWS) + len(HEURISTIC_LAWS)
    assert len(LAW_CANON) == total_derived, f"Law mismatch: Canon({len(LAW_CANON)}) != Tiers({total_derived})"
    
    # 2. AUTHORITATIVE Check
    assert DETERMINISTIC_COUNT == len(DETERMINISTIC_LAWS), "Denominator mismatch for veto gate"
    
    print(f"✅ Audit Success: {len(LAW_CANON)} Laws validated ({DETERMINISTIC_COUNT} Authority, {len(HEURISTIC_LAWS)} Advisory).")

if __name__ == "__main__":
    try:
        run_contract_audit()
    except Exception as e:
        print(f"❌ CONTRACT VIOLATION: {e}")
        sys.exit(1)
