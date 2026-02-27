"""
Toscanini Phase B2 — Credits + Cost Test Suite
Run: python3 tests/test_credits.py
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

print("\n=== CREDITS + COST TEST SUITE (Week 11) ===\n")

def t1():
    from worker.credits import get_credits, BETA_CREDITS_EMAIL
    import secrets
    uid = f"test_{secrets.token_hex(4)}@test.com"
    rec = get_credits(uid)
    assert rec["credits_total"]     == BETA_CREDITS_EMAIL
    assert rec["credits_remaining"] == BETA_CREDITS_EMAIL
    assert rec["tier"]              == "beta"
run_test("New email user gets beta credits", t1)

def t2():
    from worker.credits import get_credits, BETA_CREDITS_ANON
    import secrets
    ip = f"192.168.{secrets.token_hex(1)}.{secrets.token_hex(1)}"
    rec = get_credits(ip)
    assert rec["credits_total"]     == BETA_CREDITS_ANON
    assert rec["credits_remaining"] == BETA_CREDITS_ANON
run_test("Anonymous user gets anon credits", t2)

def t3():
    from worker.credits import check_credits
    import secrets
    uid = f"test_{secrets.token_hex(4)}@test.com"
    result = check_credits(uid, "openmm")
    assert result["allowed"]  is True
    assert result["cost"]     == 5
run_test("Credit check allows new user", t3)

def t4():
    from worker.credits import get_credits, deduct_credits, BETA_CREDITS_EMAIL
    import secrets
    uid = f"test_{secrets.token_hex(4)}@test.com"
    get_credits(uid)
    deduct_credits(uid, "openmm", "JOB001")
    rec = get_credits(uid)
    assert rec["credits_remaining"] == BETA_CREDITS_EMAIL - 1
    assert rec["jobs_submitted"]    == 1
run_test("Credit deduction works", t4)

def t5():
    from worker.credits import get_credits, deduct_credits, check_credits
    import secrets
    uid = f"exhausted_{secrets.token_hex(4)}@test.com"
    get_credits(uid)
    # Exhaust all credits
    for i in range(10):
        deduct_credits(uid, "openmm", f"JOB{i:03d}")
    result = check_credits(uid, "openmm")
    assert result["allowed"]           is False
    assert result["credits_remaining"] == 0
run_test("Exhausted credits blocked", t5)

def t6():
    from worker.credits import get_usage_stats
    import secrets
    uid = f"stats_{secrets.token_hex(4)}@test.com"
    stats = get_usage_stats(uid)
    assert "credits_remaining" in stats
    assert "jobs_submitted"    in stats
    assert "gpu_minutes_used"  in stats
run_test("Usage stats structure", t6)

def t7():
    from worker.cost_tracker import estimate_cost
    est = estimate_cost("openmm")
    assert est["estimated_minutes"] == 5.0
    assert est["estimated_usd"]     > 0
    assert "gpu_instance"           in est
run_test("Cost estimation", t7)

def t8():
    from worker.cost_tracker import estimate_cost
    costs = {p: estimate_cost(p)["estimated_usd"] for p in ["openmm","rosetta","both","loop"]}
    assert costs["both"] > costs["openmm"]
    assert costs["loop"] > costs["both"]
run_test("Cost ordering (loop > both > openmm)", t8)

def t9():
    from worker.cost_tracker import record_job_cost, get_total_costs
    import time, secrets
    jid = "COST_" + secrets.token_hex(4).upper()
    t0  = time.time() - 300  # 5 minutes ago
    t1  = time.time()
    rec = record_job_cost(jid, "openmm", t0, t1, success=True)
    assert rec["job_id"]           == jid
    assert rec["duration_minutes"] > 0
    assert rec["actual_usd"]       > 0
run_test("Job cost recording", t9)

def t10():
    from worker.cost_tracker import get_total_costs
    summary = get_total_costs(days=30)
    assert "total_jobs"        in summary
    assert "total_usd"         in summary
    assert "total_gpu_minutes" in summary
    assert "by_protocol"       in summary
run_test("Cost aggregation structure", t10)

print(f"\n=== RESULTS: {passed}/{passed+failed} passed ===\n")
if failed > 0:
    sys.exit(1)
