"""
Test: Remediation ZIP Determinism

Verifies that the same audit result produces byte-identical ZIP output.
This is a critical governance requirement — artifacts must be deterministic.
"""

import hashlib
from remediation_generator import generate_remediation_zip


def test_determinism():
    """Generate same audit twice, confirm byte-identical ZIPs."""
    
    # Fixed audit result (no randomness)
    audit = {
        'verdict': {'binary': 'VETO', 'coverage_pct': 65.0, 'deterministic_score': 72,
                    'det_total': 12, 'det_passed': 7},
        'governance': {'audit_id': 'DETERMINISM_TEST', 'timestamp_utc': '2025-01-01T00:00:00Z'},
        'provenance': {'source': 'test_structure'},
        'characterization': {'total_residues': 150},
        'tier1': {'laws': [
            {'law_id': 'LAW-125', 'title': 'Ramachandran', 'method': 'deterministic',
             'observed': 8.2, 'threshold': 5.0, 'operator': '<=', 'status': 'VETO'},
            {'law_id': 'LAW-150', 'title': 'Rotamer Audit', 'method': 'deterministic',
             'observed': 24.1, 'threshold': 20.0, 'operator': '<=', 'status': 'VETO'},
            {'law_id': 'LAW-130', 'title': 'Clashscore', 'method': 'deterministic',
             'observed': 28.5, 'threshold': 20.0, 'operator': '<=', 'status': 'VETO'},
            {'law_id': 'LAW-110', 'title': 'Backbone Gap', 'method': 'deterministic',
             'observed': 1, 'threshold': 0, 'operator': '=', 'status': 'VETO'},
        ]},
    }

    # Generate twice
    zip1 = generate_remediation_zip(audit)
    zip2 = generate_remediation_zip(audit)

    # Compute SHA-256 hashes
    hash1 = hashlib.sha256(zip1).hexdigest()
    hash2 = hashlib.sha256(zip2).hexdigest()

    # Verify
    assert hash1 == hash2, f"NON-DETERMINISTIC: {hash1} != {hash2}"
    assert len(zip1) == len(zip2), f"SIZE MISMATCH: {len(zip1)} != {len(zip2)}"
    assert zip1 == zip2, "BYTE-LEVEL MISMATCH"

    print(f"✅ DETERMINISM TEST PASSED")
    print(f"   SHA-256: {hash1}")
    print(f"   Size: {len(zip1)} bytes")
    print(f"   Files: remediation_report.json, rosetta_relax.xml, rosetta.flags,")
    print(f"          loop_modeling.xml, openmm_equilibrate.py, README.txt")


if __name__ == "__main__":
    test_determinism()
