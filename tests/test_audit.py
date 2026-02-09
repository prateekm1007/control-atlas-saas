import pytest
from tools.native_audit import SovereignJudge

def test_contact_density_calculation():
    judge = SovereignJudge()
    result = judge.calculate_contact_density(
        "structures/champ005/structure.cif", "A", "B"
    )
    assert result['rho'] == 199
    assert result['clashes'] == 7
    assert result['min_distance_A'] == 1.63

def test_hash_integrity():
    import hashlib
    with open("structures/champ005/structure.cif", "rb") as f:
        sha = hashlib.sha256(f.read()).hexdigest()
    assert sha == "8c0e283004eda5a18f711fc0c3e5e59ea26d49571af1a8516e4136386002e131"
