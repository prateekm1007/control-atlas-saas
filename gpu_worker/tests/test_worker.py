"""
Toscanini Phase B2 — Worker Test Suite (Week 7)
Run: python3 tests/test_worker.py
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

print("\n=== WORKER TEST SUITE (Week 7) ===\n")

def t1():
    from worker.tasks import execute_openmm, execute_rosetta, auto_callback
    assert all(x is not None for x in [execute_openmm, execute_rosetta, auto_callback])
run_test("Task imports", t1)

def t2():
    from worker.tasks import app
    assert app.conf.task_serializer == "json"
    assert app.conf.worker_prefetch_multiplier == 1
run_test("Celery config", t2)

def t3():
    from worker.tasks import execute_openmm
    import inspect
    params = list(inspect.signature(execute_openmm).parameters.keys())
    assert "job_id" in params and "pdb_bytes_hex" in params
run_test("OpenMM task signature", t3)

def t4():
    from worker.tasks import execute_rosetta
    import inspect
    params = list(inspect.signature(execute_rosetta).parameters.keys())
    assert "job_id" in params and "xml_protocol" in params
run_test("Rosetta task signature", t4)

def t5():
    from worker.tasks import auto_callback
    result = auto_callback({"status": "failed", "job_id": "T5"}, "AUD001")
    assert result["status"] == "skipped"
run_test("Auto-callback skips failed", t5)

def t6():
    pdb = b"ATOM      1  N   ALA A   1       0.000   0.000   0.000\nEND\n"
    assert bytes.fromhex(pdb.hex()) == pdb
run_test("PDB hex round-trip", t6)

def t7():
    from worker.tasks import MAX_JOB_SECONDS
    assert MAX_JOB_SECONDS == 1800
run_test("30 min timeout", t7)

def t8():
    from worker.job_state import (STATE_QUEUED, STATE_RUNNING,
                                   STATE_SUCCESS, STATE_FAILED, STATE_TIMEOUT)
    assert STATE_QUEUED == "queued"
    assert STATE_SUCCESS == "success"
    assert STATE_FAILED == "failed"
run_test("State constants", t8)

def t9():
    from worker.job_state import create_job, get_job, update_job
    import inspect
    for fn in [create_job, get_job, update_job]:
        assert "job_id" in inspect.signature(fn).parameters
run_test("State manager API", t9)

def t10():
    from worker.health import full_health_check
    result = full_health_check()
    assert all(k in result for k in ["redis", "gpu", "openmm", "worker_ready"])
    assert result["worker_ready"] is True
run_test("Health monitor structure", t10)

def t11():
    from worker.health import check_openmm
    result = check_openmm()
    assert "status" in result
run_test("OpenMM health check", t11)

def t12():
    from worker.job_state import get_job_duration
    import inspect
    assert "job" in inspect.signature(get_job_duration).parameters
run_test("Job duration calculator", t12)

print(f"\n=== RESULTS: {passed}/{passed+failed} passed ===\n")
if failed > 0:
    sys.exit(1)
