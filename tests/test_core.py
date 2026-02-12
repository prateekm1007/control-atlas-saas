import json
from pathlib import Path

def test_structures():
    assert Path("structures/champ005/structure.cif").exists()
    assert Path("structures/twinrod_v2/structure.cif").exists()

def test_metrics():
    with open("structures/validation_summary.json") as f:
        data = json.load(f)
    assert data["leads"]["twinrod_v2"]["post_minimization"]["clashes"] == 0
