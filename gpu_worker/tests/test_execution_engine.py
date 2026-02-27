"""
Toscanini Phase B2 — Execution Engine Test Suite (Week 8)
Run: python3 tests/test_execution_engine.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

passed = 0
failed = 0

def run_test(name, fn):
    global passed, failed
    try:
        fn()
        print(f"  ✓ {name}")
        passed += 1
    except Exception as e:
        print(f"  ✗ {name}: {e}")
        failed += 1

print("\n=== EXECUTION ENGINE TEST SUITE (Week 8) ===\n")

# Protocol selector tests
def t1():
    from worker.protocol_selector import select_protocol, PROTOCOL_ROSETTA
    result = select_protocol(["LAW-125", "LAW-150"])
    assert result["protocol"] == PROTOCOL_ROSETTA
run_test("Rama+Rotamer → Rosetta", t1)

def t2():
    from worker.protocol_selector import select_protocol, PROTOCOL_OPENMM
    result = select_protocol(["LAW-182", "LAW-200"])
    assert result["protocol"] == PROTOCOL_OPENMM
run_test("Burial+Packing → OpenMM", t2)

def t3():
    from worker.protocol_selector import select_protocol, PROTOCOL_BOTH
    result = select_protocol(["LAW-125", "LAW-182"])
    assert result["protocol"] == PROTOCOL_BOTH
run_test("Rama+Burial → Both", t3)

def t4():
    from worker.protocol_selector import select_protocol, PROTOCOL_LOOP
    result = select_protocol(["LAW-110"])
    assert result["protocol"] == PROTOCOL_LOOP
run_test("Backbone gap → Loop", t4)

def t5():
    from worker.protocol_selector import select_protocol, PROTOCOL_OPENMM
    result = select_protocol([])
    assert result["protocol"] == PROTOCOL_OPENMM
run_test("No violations → OpenMM default", t5)

def t6():
    from worker.protocol_selector import get_openmm_config
    config = get_openmm_config(["LAW-182"])
    assert config["sim_steps"] == 2500000  # 5ns for burial
run_test("OpenMM 5ns for burial violations", t6)

def t7():
    from worker.protocol_selector import get_openmm_config
    config = get_openmm_config(["LAW-125"])
    assert config["sim_steps"] == 1000000  # 2ns default
run_test("OpenMM 2ns for geometry violations", t7)

def t8():
    from worker.protocol_selector import get_rosetta_scoreterms
    weights = get_rosetta_scoreterms(["LAW-125", "LAW-150", "LAW-130"])
    assert weights["rama_weight"] == "2.0"
    assert weights["dun_weight"]  == "2.0"
    assert weights["rep_weight"]  == "1.5"
run_test("Rosetta weights upregulated for violations", t8)

def t9():
    from worker.protocol_selector import get_rosetta_scoreterms
    weights = get_rosetta_scoreterms([])
    assert weights["rama_weight"] == "1.0"
    assert weights["dun_weight"]  == "1.0"
run_test("Rosetta default weights when no violations", t9)

def t10():
    from worker.execution_engine import generate_job_id
    ids = {generate_job_id() for _ in range(100)}
    assert len(ids) == 100  # All unique
run_test("Job ID uniqueness (100 samples)", t10)

def t11():
    from worker.execution_engine import _generate_rosetta_xml
    xml = _generate_rosetta_xml(
        "TESTAUDIT",
        ["LAW-125", "LAW-150"],
        {"rama_weight": "2.0", "dun_weight": "2.0",
         "rep_weight": "0.55", "omega_weight": "1.0",
         "geom_weight": "1.0", "dslf_weight": "1.0"}
    )
    assert "<ROSETTASCRIPTS>" in xml
    assert "TESTAUDIT" in xml
    assert 'weight="2.0"' in xml
run_test("Rosetta XML generation", t11)

def t12():
    from worker.protocol_selector import select_protocol
    result = select_protocol(["LAW-125", "LAW-150", "LAW-130",
                               "LAW-182", "LAW-200"])
    assert "estimated_minutes" in result
    assert result["estimated_minutes"] > 0
run_test("Protocol selector returns time estimate", t12)

print(f"\n=== RESULTS: {passed}/{passed+failed} passed ===\n")
if failed > 0:
    sys.exit(1)
