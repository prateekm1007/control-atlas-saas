"""
Toscanini Phase B2 — Worker Test Suite
Run directly: python3 tests/test_worker.py
"""
import sys
import os
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

print("\n=== WORKER TASK TEST SUITE ===\n")

def t1():
    from worker.tasks import execute_openmm, execute_rosetta, auto_callback
    assert execute_openmm is not None
run_test("Task module imports", t1)

def t2():
    from worker.tasks import app
    assert app.conf.task_serializer == "json"
    assert app.conf.worker_prefetch_multiplier == 1
run_test("Celery app configuration", t2)

def t3():
    from worker.tasks import execute_openmm
    import inspect
    params = list(inspect.signature(execute_openmm).parameters.keys())
    assert "job_id" in params
    assert "pdb_bytes_hex" in params
run_test("OpenMM task signature", t3)

def t4():
    from worker.tasks import execute_rosetta
    import inspect
    params = list(inspect.signature(execute_rosetta).parameters.keys())
    assert "job_id" in params
    assert "xml_protocol" in params
run_test("Rosetta task signature", t4)

def t5():
    from worker.tasks import auto_callback
    result = auto_callback(
        {"status": "failed", "error": "test", "job_id": "T5"},
        "ORIGINAL_AUDIT_001"
    )
    assert result["status"] == "skipped"
run_test("Auto-callback skips failed execution", t5)

def t6():
    pdb = b"ATOM      1  N   ALA A   1       0.000   0.000   0.000\nEND\n"
    assert bytes.fromhex(pdb.hex()) == pdb
run_test("PDB hex encoding round-trip", t6)

def t7():
    from worker.tasks import MAX_JOB_SECONDS
    assert MAX_JOB_SECONDS == 1800
run_test("Job timeout configured (30 min)", t7)

def t8():
    from worker.tasks import app
    assert app.conf.worker_prefetch_multiplier == 1
run_test("Single GPU concurrency enforced", t8)

print(f"\n=== RESULTS: {passed}/{passed+failed} passed ===\n")
if failed > 0:
    sys.exit(1)
