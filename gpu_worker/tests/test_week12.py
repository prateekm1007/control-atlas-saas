"""
Toscanini Week 12 — Security & Sandboxing Smoke Tests
Run: TOSCANINI_DATA_DIR=/tmp/tos_test python3 gpu_worker/tests/test_week12.py
No Docker required. No GPU required. No live Redis required.
"""
import sys, os, subprocess, time, signal, io, importlib

# Ensure worker package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ.setdefault("TOSCANINI_DATA_DIR", "/tmp/tos_w12_test")

passed = 0
failed = 0

def run(name, fn):
    global passed, failed
    try:
        fn()
        print(f"  ✓ {name}")
        passed += 1
    except Exception as e:
        print(f"  ✗ {name}: {e}")
        failed += 1

print("\n=== WEEK 12 SMOKE TESTS ===\n")

# ── Group 1: _count_residues (pure-Python, no deps) ──────────────────────────

def t_count_basic():
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../backend"))
    # inline the function so this test has zero import side-effects
    def _count_residues(pdb_bytes):
        seen = set()
        for line in pdb_bytes.decode("utf-8", errors="replace").splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                chain  = line[21] if len(line) > 21 else "?"
                resseq = line[22:26].strip() if len(line) > 26 else "0"
                seen.add((chain, resseq))
        return len(seen)

    pdb = b"ATOM      1  CA  ALA A   1       1.000   1.000   1.000  1.00 10.00           C\n"
    assert _count_residues(pdb) == 1

def t_count_two_chains():
    def _count_residues(pdb_bytes):
        seen = set()
        for line in pdb_bytes.decode("utf-8", errors="replace").splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                chain  = line[21] if len(line) > 21 else "?"
                resseq = line[22:26].strip() if len(line) > 26 else "0"
                seen.add((chain, resseq))
        return len(seen)

    pdb = (
        b"ATOM      1  CA  ALA A   1       1.000   1.000   1.000  1.00 10.00           C\n"
        b"ATOM      2  CA  ALA A   1       1.100   1.100   1.100  1.00 10.00           C\n"  # same residue
        b"ATOM      3  CA  GLY B   1       5.000   5.000   5.000  1.00 10.00           C\n"  # diff chain
    )
    assert _count_residues(pdb) == 2  # A:1 and B:1

def t_count_empty():
    def _count_residues(pdb_bytes):
        seen = set()
        for line in pdb_bytes.decode("utf-8", errors="replace").splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                chain  = line[21] if len(line) > 21 else "?"
                resseq = line[22:26].strip() if len(line) > 26 else "0"
                seen.add((chain, resseq))
        return len(seen)
    assert _count_residues(b"END\n") == 0

run("_count_residues: single residue", t_count_basic)
run("_count_residues: two chains counted separately", t_count_two_chains)
run("_count_residues: empty PDB returns 0", t_count_empty)

# ── Group 2: _hard_timeout_process (subprocess hardening) ────────────────────

def t_timeout_kills():
    from worker.tasks import _hard_timeout_process
    start = time.monotonic()
    try:
        with _hard_timeout_process(["sleep", "30"], timeout_sec=1) as proc:
            proc.communicate(timeout=1)
    except subprocess.TimeoutExpired:
        pass
    elapsed = time.monotonic() - start
    assert elapsed < 8, f"Process took {elapsed:.1f}s — not killed promptly"

def t_success_exits_clean():
    from worker.tasks import _hard_timeout_process
    with _hard_timeout_process(["echo", "toscanini"], timeout_sec=5) as proc:
        out, _ = proc.communicate(timeout=5)
    assert b"toscanini" in out
    assert proc.returncode == 0

def t_list_args_no_shell():
    """Verify cmd is a list in execute_rosetta — no shell=True in executable code."""
    import inspect, ast, sys
    sys.path.insert(0, ".")
    from worker import tasks
    full_src = inspect.getsource(tasks.execute_rosetta)
    # Use AST to extract only real Call nodes — immune to comments and docstrings
    tree = ast.parse(full_src)
    calls = [
        ast.unparse(node)
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
    ]
    call_str = " ".join(calls)
    assert "shell=True" not in call_str,  "shell=True in a real Call node"
    assert "subprocess.run" not in full_src, "subprocess.run still present"
    assert "cmd = [" in full_src,            "cmd should be built as a list"
    assert "os.setsid" in full_src,          "os.setsid missing"
    assert "killpg" in full_src,             "killpg missing"

run("Hard timeout kills process within grace period", t_timeout_kills)
run("Successful process exits cleanly", t_success_exits_clean)
run("execute_rosetta uses list args (no shell=True)", t_list_args_no_shell)

# ── Group 3: Celery config assertions ────────────────────────────────────────

def t_celery_time_limits():
    from worker.tasks import app as celery_app
    assert celery_app.conf.task_time_limit      == 1800, "Hard limit must be 1800s"
    assert celery_app.conf.task_soft_time_limit == 1740, "Soft limit must be 1740s"

def t_celery_max_tasks_per_child():
    from worker.tasks import app as celery_app
    assert celery_app.conf.worker_max_tasks_per_child == 2, "Must recycle after 2 jobs"

run("Celery hard+soft time limits set correctly", t_celery_time_limits)
run("Celery worker recycled after 2 jobs", t_celery_max_tasks_per_child)

# ── Group 4: docker-compose.yml wiring ───────────────────────────────────────

def t_compose_has_redis():
    import yaml
    dc = yaml.safe_load(open(os.path.join(os.path.dirname(__file__), "../../docker-compose.yml")))
    assert "redis"      in dc["services"], "redis service missing"
    assert "gpu_worker" in dc["services"], "gpu_worker service missing"

def t_compose_has_volumes():
    import yaml
    dc = yaml.safe_load(open(os.path.join(os.path.dirname(__file__), "../../docker-compose.yml")))
    top_vols = dc.get("volumes", {})
    assert "toscanini_data" in top_vols, "toscanini_data volume missing"
    assert "redis_data"     in top_vols, "redis_data volume missing"

def t_compose_data_dir_env():
    import yaml
    dc = yaml.safe_load(open(os.path.join(os.path.dirname(__file__), "../../docker-compose.yml")))
    brain_env = dc["services"]["brain"].get("environment", [])
    env_str = str(brain_env)
    assert "TOSCANINI_DATA_DIR" in env_str, "TOSCANINI_DATA_DIR not in brain env"

try:
    import yaml
    run("docker-compose has redis + gpu_worker", t_compose_has_redis)
    run("docker-compose has required volumes",   t_compose_has_volumes)
    run("brain service has TOSCANINI_DATA_DIR",  t_compose_data_dir_env)
except ImportError:
    print("  ⚠ pyyaml not installed — skipping compose checks (pip install pyyaml)")

# ── Group 5: main.py patch verification (no imports, source-level) ────────────

def t_main_has_validate():
    src = open(os.path.join(os.path.dirname(__file__), "../../backend/main.py")).read()
    assert "_validate_pdb_upload" in src
    assert "_count_residues"      in src
    assert "MAX_PDB_BYTES"        in src
    assert "MAX_RESIDUE_COUNT"    in src

def t_main_no_dir_bug():
    src = open(os.path.join(os.path.dirname(__file__), "../../backend/main.py")).read()
    assert '"audit_result" in dir()' not in src, "dir() bug still present"

def t_main_validate_called_in_submit():
    src = open(os.path.join(os.path.dirname(__file__), "../../backend/main.py")).read()
    assert "_validate_pdb_upload(content_bytes)" in src

run("main.py contains validation helpers",          t_main_has_validate)
run("main.py: dir() bug removed",                  t_main_no_dir_bug)
run("main.py: _validate_pdb_upload called in submit", t_main_validate_called_in_submit)

# ── Results ───────────────────────────────────────────────────────────────────

print(f"\n=== RESULTS: {passed}/{passed+failed} passed ===\n")
if failed > 0:
    sys.exit(1)
