#!/usr/bin/env python3
"""Acceptance tests AT-01..AT-16 for the EC-16..EC-26 enforcement scaffolding.

EC-20 dogfood: these are real pytest node IDs; the package names them in
ACCEPTANCE_TESTS.json and `verify_claimed_shas.py --check-tests` executes them on the
exact candidate. Plain stdlib + PyYAML; every checker runs as a subprocess with
argument arrays (no shell). Exit status carries the verdict; counts are printed.
"""
import json
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
CI = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(CI))  # repo root (branch checkout)
PY = sys.executable

PASS, FAIL = [], []


def run(script, *args, cwd=REPO):
    return subprocess.run([PY, os.path.join(CI, script), *args],
                          capture_output=True, text=True, cwd=cwd, shell=False)


def check(at_id, name, cond, detail=""):
    (PASS if cond else FAIL).append(f"{at_id} {name}" + (f" -- {detail}" if not cond else ""))


def tmpfile(content: str, suffix=".txt") -> str:
    f = tempfile.NamedTemporaryFile("w", suffix=suffix, delete=False)
    f.write(content)
    f.close()
    return f.name


def rev(args, cwd=REPO):
    return subprocess.run(["git", *args], capture_output=True, text=True, cwd=cwd).stdout.strip()


# ---------------------------------------------------------------- AT-01 (EC-16/19)

def test_at01_claimed_shas():
    head = rev(["rev-parse", "HEAD"])
    good = tmpfile(f"Claimed-SHAs:\n  {head}\nDirective-ID: TEST\nCloses-items: TEST\nStatus-summary: PUSHED\n")
    bad = tmpfile("Claimed-SHAs:\n  " + "0" * 40 + "\n")
    r = run("verify_claimed_shas.py", "--body", good, "--local")
    check("AT-01", "real HEAD passes offline verification", r.returncode == 0, r.stdout)
    r = run("verify_claimed_shas.py", "--body", bad, "--local")
    check("AT-01", "absent SHA fails", r.returncode == 1 and "not on remote" in r.stdout, r.stdout)
    none = tmpfile("no claimed shas here\n")
    r = run("verify_claimed_shas.py", "--body", none, "--local")
    check("AT-01", "missing field fails closed", r.returncode == 1 and "field missing" in r.stdout, r.stdout)
    explicit = tmpfile("Claimed-SHAs: none\nDirective-ID: TEST\nCloses-items: T\nStatus-summary: NOT_RUN\n")
    r = run("verify_claimed_shas.py", "--body", explicit, "--local")
    check("AT-01", "explicit none declaration passes", r.returncode == 0, r.stdout)
    mixed = tmpfile("Claimed-SHAs: none\n  " + head + "\n")
    r = run("verify_claimed_shas.py", "--body", mixed, "--local")
    check("AT-01", "none mixed with SHAs fails closed", r.returncode == 1 and "mixes" in r.stdout, r.stdout)
    r = run("verify_claimed_shas.py", "--body", bad, "--require-directive-id-only")
    check("AT-01", "directive-only mode ignores SHA absence, flags missing Directive-ID",
          r.returncode == 1 and "Directive-ID" in r.stdout and "not on remote" not in r.stdout, r.stdout)


# ---------------------------------------------------------------- AT-02 (EC-21 local)

def test_at02_declared_base():
    head = rev(["rev-parse", "HEAD"])
    base = rev(["rev-parse", "HEAD~1"])
    ok = tmpfile('{"status": 200, "body": {"sha": "%s"}}' % base)
    refs = tmpfile(json.dumps({"main-line": "HEAD~1", "science-line": "HEAD", "policy-line": "HEAD"}))
    r = run("check_declared_base.py", "--base", base, "--head", head,
            "--stack-parent", "main-line", "--repo-slug", "x/y",
            "--api-json", ok, "--stack-refs", refs)
    check("AT-02", "ancestor base passes with matching stack", r.returncode == 0, r.stdout)
    r = run("check_declared_base.py", "--base", "0" * 40, "--head", head,
            "--stack-parent", "main-line", "--repo-slug", "x/y",
            "--api-json", '{"status": 404}', "--stack-refs", refs)
    check("AT-02", "absent base fails closed", r.returncode == 1 and "LOCAL_OBJECT_ABSENT" in r.stdout, r.stdout)


# ---------------------------------------------------------------- AT-03 (EC-25 lint)

def at03_lint_dir():
    bad = tmpfile("""
name: bad
on:
  push:
    branches: [main]
jobs:
  deploy:
    runs-on: ubuntu-latest
    steps: [{run: echo}]
""")
    ok = tmpfile("""
name: good
on:
  pull_request:
jobs:
  deploy:
    runs-on: ubuntu-latest
    steps: [{run: echo}]
""")
    ok2 = tmpfile("""
name: good2
on:
  workflow_run:
    workflows: [x]
jobs:
  unit-tests:
    runs-on: ubuntu-latest
    steps: [{run: echo}]
""")
    d = tempfile.mkdtemp()
    for i, src in enumerate([bad, ok, ok2]):
        with open(src) as fh:
            content = fh.read()
        with open(os.path.join(d, f"w{i}.yml"), "w") as fh:
            fh.write(content)
    return d


def test_at03_deploy_lint():
    d = at03_lint_dir()
    r = run("workflow_deploy_trigger_lint.py", d)
    check("AT-03", "deploy-on-push flagged", r.returncode == 1 and "w0.yml" in r.stdout, r.stdout)
    check("AT-03", "deploy-on-PR and non-deploy jobs not flagged", "w1.yml" not in r.stdout
          and "w2.yml" not in r.stdout, r.stdout)
    r = run("workflow_deploy_trigger_lint.py", d, "--report-only")
    check("AT-03", "--report-only exits 0 with findings listed",
          r.returncode == 0 and "w0.yml" in r.stdout, r.stdout)


# ---------------------------------------------------------------- AT-04 (EC-22 basic)

def test_at04_status_vocabulary():
    bad = tmpfile("STATUS: RUNNING\nSTATUS: DONE\nSTATUS: RUN_PASS\n"
                  "Gitleaks came back clean scan\nFailures are pre-existing\nProgress: 92%\n")
    r = run("status_vocabulary_lint.py", bad, "--runs-manifest", '{"runs": {}}')
    check("AT-04", "unknown status tokens flagged", "RUNNING" in r.stdout and "DONE" in r.stdout, r.stdout)
    check("AT-04", "RUN_PASS without id/args flagged", "without run_id" in r.stdout, r.stdout)
    check("AT-04", "clean scan without citation flagged", "'clean scan'" in r.stdout, r.stdout)
    check("AT-04", "pre-existing without citation flagged", "'pre-existing'" in r.stdout, r.stdout)
    check("AT-04", "percentage without run records flagged", "not supported by run records" in r.stdout, r.stdout)
    good = tmpfile("STATUS: NOT_RUN\nSTATUS: RUN_FAIL\nSTATUS: PUSHED\nSTATUS: MERGED\n")
    r = run("status_vocabulary_lint.py", good, "--runs-manifest", '{"runs": {}}')
    check("AT-04", "non-RUN_PASS vocabulary passes", r.returncode == 0, r.stdout)


# ---------------------------------------------------------------- AT-05 (EC-18 sequence)

def test_at05_sequence():
    bad = tmpfile("## 2026-10-20T10:00:00Z | E2 | action | second\n"
                  "## 2026-10-20T09:00:00Z | E1 | action | first\n")
    backdated = tmpfile("## 2026-10-19T00:00:00Z | E1 | read | marker\n")
    recon = tmpfile("## 2026-10-19T00:00:00Z | E1 | read | marker | RECONSTRUCTED\n")
    r = run("evidence_sequence_check.py", bad)
    check("AT-05", "backwards timestamps flagged", r.returncode == 1 and "backwards" in r.stdout, r.stdout)
    r = run("evidence_sequence_check.py", backdated, "--commit-time", "2026-10-20T12:00:00Z")
    check("AT-05", "backdated entry without label flagged",
          r.returncode == 1 and "RECONSTRUCTED" in r.stdout, r.stdout)
    r = run("evidence_sequence_check.py", recon, "--commit-time", "2026-10-20T12:00:00Z")
    check("AT-05", "RECONSTRUCTED backfill accepted", r.returncode == 0, r.stdout)


# ---------------------------------------------------------------- AT-06 (EC-26 positive)

def _write_mini_index(d, rules):
    reg = {"tests": {f"AT-{i:02d}": {"name": f"t{i}", "node": f"x.py::test_t{i}"} for i in range(1, 9)}}
    open(os.path.join(d, "registry.json"), "w").write(json.dumps(reg))
    os.makedirs(os.path.join(d, "docs"), exist_ok=True)
    rows = "\n".join(f"| {rid} | {st} | procedural:AGENTS.md | AT-01 |" for rid, st in rules)
    idx = f"""# index
<!-- BEGIN RULE INVENTORY -->
reading-note-path: docs/note.md
pinned-blobs: b6bb35e84956d84a8c9d9e65734c212baf3ed075
| Rule ID | Status | Enforcement point(s) | Acceptance test(s) |
|---|---|---|---|
{rows}
<!-- END RULE INVENTORY -->
"""
    os.makedirs(os.path.join(d, "docs", "directives"), exist_ok=True)
    open(os.path.join(d, "docs", "directives", "INDEX.md"), "w").write(idx)
    open(os.path.join(d, "AGENTS.md"), "w").write("# agents\n")
    return idx


def test_at06_reading_note():
    d = tempfile.mkdtemp()
    _write_mini_index(d, [("EC-16", "PROPOSED"), ("EC-17", "PROPOSED")])
    os.makedirs(os.path.join(d, "docs"), exist_ok=True)
    note = ("blob b6bb35e84956d84a8c9d9e65734c212baf3ed075 read\n"
            "| EC-16 | -> check A |\n| EC-17 | -> check B |\n")
    open(os.path.join(d, "docs", "note.md"), "w").write(note)
    r = run("reading_note_map_check.py", "--index", "docs/directives/INDEX.md",
            "--repo-root", d, "--registry", "registry.json", "--include-proposed")
    check("AT-06", "complete note passes with proposed included", r.returncode == 0, r.stdout)
    r2 = run("reading_note_map_check.py", "--index", "docs/directives/INDEX.md",
             "--repo-root", d, "--registry", "registry.json")
    check("AT-06", "without --include-proposed the ACTIVE-only set is empty -> fail closed",
          r2.returncode == 1 and "empty" in r2.stdout.lower(), r2.stdout)


# ---------------------------------------------------------------- AT-07 (PR template)

def test_at07_pr_template():
    tpl = open(os.path.join(REPO, "proposals", "ci", "PULL_REQUEST_TEMPLATE.md"),
               encoding="utf-8").read()
    for field in ("Directive-ID", "Claimed-SHAs", "Base-SHA", "Stack-parent", "Acceptance-Tests", "Status-summary"):
        check("AT-07", f"PR template carries {field}", field + ":" in tpl)


# ---------------------------------------------------------------- AT-08 (EC-26 fail-closed set)

def test_at08_reading_note_failclosed():
    d = tempfile.mkdtemp()
    idx = _write_mini_index(d, [("EC-16", "PROPOSED")])
    os.makedirs(os.path.join(d, "docs"), exist_ok=True)
    open(os.path.join(d, "AGENTS.md"), "w").write("# agents\n")

    def run_idx(**kw):
        base = ["reading_note_map_check.py", "--index", "docs/directives/INDEX.md",
                "--repo-root", d, "--registry", "registry.json", "--include-proposed"]
        return run(*base, **kw)

    # 1. missing note
    r = run_idx()
    check("AT-08", "missing note fails closed", r.returncode == 1 and "missing at canonical path" in r.stdout, r.stdout)
    # 2. note missing blob citation
    open(os.path.join(d, "docs", "note.md"), "w").write("| EC-16 | -> c |\n")
    r = run_idx()
    check("AT-08", "note missing blob citation fails", r.returncode == 1 and "pinned blob" in r.stdout, r.stdout)
    # 3. note missing rule mapping
    open(os.path.join(d, "docs", "note.md"), "w").write("blob b6bb35e84956d84a8c9d9e65734c212baf3ed075 read\n")
    r = run_idx()
    check("AT-08", "note missing rule mapping fails", r.returncode == 1 and "does not map required rule EC-16" in r.stdout, r.stdout)
    # 4. empty rule set
    d2 = tempfile.mkdtemp()
    _write_mini_index(d2, [])
    os.makedirs(os.path.join(d2, "docs"), exist_ok=True)
    open(os.path.join(d2, "AGENTS.md"), "w").write("# a\n")
    open(os.path.join(d2, "docs", "note.md"), "w").write(
        "blob b6bb35e84956d84a8c9d9e65734c212baf3ed075\n")
    r = run("reading_note_map_check.py", "--index", "docs/directives/INDEX.md",
            "--repo-root", d2, "--registry", "registry.json", "--include-proposed")
    check("AT-08", "empty inventory fails closed", r.returncode == 1 and "empty" in r.stdout.lower(), r.stdout)
    # 5. duplicate IDs
    d3 = tempfile.mkdtemp()
    _write_mini_index(d3, [("EC-16", "PROPOSED")])
    os.makedirs(os.path.join(d3, "docs"), exist_ok=True)
    p = os.path.join(d3, "docs", "directives", "INDEX.md")
    t = open(p).read().replace("<!-- END RULE INVENTORY",
                               "| EC-16 | PROPOSED | procedural:AGENTS.md | AT-01 |\n<!-- END RULE INVENTORY")
    open(p, "w").write(t)
    open(os.path.join(d3, "AGENTS.md"), "w").write("# a\n")
    open(os.path.join(d3, "docs", "note.md"), "w").write(
        "blob b6bb35e84956d84a8c9d9e65734c212baf3ed075\n| EC-16 | x |\n")
    r = run("reading_note_map_check.py", "--index", "docs/directives/INDEX.md",
            "--repo-root", d3, "--registry", "registry.json", "--include-proposed")
    check("AT-08", "duplicate rule IDs fail closed", r.returncode == 1 and "duplicate" in r.stdout.lower(), r.stdout)
    # 6. malformed row
    d4 = tempfile.mkdtemp()
    _write_mini_index(d4, [("EC-16", "PROPOSED")])
    os.makedirs(os.path.join(d4, "docs"), exist_ok=True)
    p = os.path.join(d4, "docs", "directives", "INDEX.md")
    t = open(p).read().replace("| EC-16 | PROPOSED | procedural:AGENTS.md | AT-01 |",
                               "| EC-16 | PROPOSED | only-three-columns |")
    open(p, "w").write(t)
    r = run("reading_note_map_check.py", "--index", "docs/directives/INDEX.md",
            "--repo-root", d4, "--registry", "registry.json", "--include-proposed")
    check("AT-08", "malformed row fails closed", r.returncode == 1 and "malformed" in r.stdout.lower(), r.stdout)
    # 7. unknown status token
    d5 = tempfile.mkdtemp()
    _write_mini_index(d5, [("EC-16", "SOMEDAY")])
    os.makedirs(os.path.join(d5, "docs"), exist_ok=True)
    open(os.path.join(d5, "AGENTS.md"), "w").write("# a\n")
    r = run("reading_note_map_check.py", "--index", "docs/directives/INDEX.md",
            "--repo-root", d5, "--registry", "registry.json", "--include-proposed")
    check("AT-08", "unknown status token fails closed", r.returncode == 1 and "unknown status" in r.stdout.lower(), r.stdout)
    # 8. named enforcement file does not exist
    d6 = tempfile.mkdtemp()
    _write_mini_index(d6, [("EC-16", "PROPOSED")])
    os.makedirs(os.path.join(d6, "docs"), exist_ok=True)
    p = os.path.join(d6, "docs", "directives", "INDEX.md")
    t = open(p).read().replace("procedural:AGENTS.md", "procedural:NO_SUCH_FILE.md")
    open(p, "w").write(t)
    open(os.path.join(d6, "docs", "note.md"), "w").write(
        "blob b6bb35e84956d84a8c9d9e65734c212baf3ed075\n| EC-16 | x |\n")
    r = run("reading_note_map_check.py", "--index", "docs/directives/INDEX.md",
            "--repo-root", d6, "--registry", "registry.json", "--include-proposed")
    check("AT-08", "missing enforcement file fails closed",
          r.returncode == 1 and "does not exist" in r.stdout, r.stdout)
    # 9. named AT not in registry
    d7 = tempfile.mkdtemp()
    _write_mini_index(d7, [("EC-16", "PROPOSED")])
    os.makedirs(os.path.join(d7, "docs"), exist_ok=True)
    p = os.path.join(d7, "docs", "directives", "INDEX.md")
    t = open(p).read().replace("| AT-01 |", "| AT-99 |")
    open(p, "w").write(t)
    open(os.path.join(d7, "AGENTS.md"), "w").write("# a\n")
    open(os.path.join(d7, "docs", "note.md"), "w").write(
        "blob b6bb35e84956d84a8c9d9e65734c212baf3ed075\n| EC-16 | x |\n")
    r = run("reading_note_map_check.py", "--index", "docs/directives/INDEX.md",
            "--repo-root", d7, "--registry", "registry.json", "--include-proposed")
    check("AT-08", "AT id missing from registry fails closed",
          r.returncode == 1 and "AT-99" in r.stdout, r.stdout)


# ---------------------------------------------------------------- AT-09 (EC-20 execution)

def test_at09_acceptance_tests_executed():
    proj = tempfile.mkdtemp()
    subprocess.run(["git", "init", "-q"], cwd=proj)
    subprocess.run(["git", "config", "user.email", "t@t"], cwd=proj)
    subprocess.run(["git", "config", "user.name", "t"], cwd=proj)
    open(os.path.join(proj, "test_ok.py"), "w").write(
        "def test_first_passes():\n    assert True\n\ndef test_second_passes():\n    assert 1 == 1\n")
    open(os.path.join(proj, "test_bad.py"), "w").write("def test_fails():\n    assert 0\n")
    subprocess.run(["git", "add", "."], cwd=proj, capture_output=True)
    subprocess.run(["git", "commit", "-qm", "c"], cwd=proj, capture_output=True)
    ev = os.path.join(proj, "evidence.json")
    # 1. missing identifier fails closed and is NOT executed
    body = tmpfile("Acceptance-Tests:\n  test_ok.py::test_first_passes\n  test_ok.py::DOES_NOT_EXIST\n")
    r = run("verify_claimed_shas.py", "--body", body, "--check-tests", "--evidence-out", ev, cwd=proj)
    check("AT-09", "missing named test fails (never silently skipped)",
          r.returncode == 1 and "not collectable" in r.stdout, r.stdout)
    # 2. failing named test fails
    body = tmpfile("Acceptance-Tests:\n  test_bad.py::test_fails\n")
    r = run("verify_claimed_shas.py", "--body", body, "--check-tests", "--evidence-out", ev, cwd=proj)
    check("AT-09", "failing named test fails the check", r.returncode == 1 and "FAILED on the candidate" in r.stdout, r.stdout)
    # 3. passing named tests produce a complete evidence record
    body = tmpfile("Acceptance-Tests:\n  test_ok.py::test_first_passes\n  test_ok.py::test_second_passes\n")
    r = run("verify_claimed_shas.py", "--body", body, "--check-tests", "--evidence-out", ev, cwd=proj)
    check("AT-09", "passing named tests pass", r.returncode == 0, r.stdout)
    e = json.load(open(ev))
    cand = rev(["rev-parse", "HEAD"], cwd=proj)
    check("AT-09", "evidence records candidate SHA", e.get("candidate_sha") == cand, json.dumps(e)[:200])
    check("AT-09", "evidence records command as argv LIST (no shell interpolation)",
          isinstance(e.get("command_argv"), list) and e.get("shell_used") is False, json.dumps(e.get("command_argv", ""))[:120])
    check("AT-09", "evidence records counts and exit status",
          e.get("counts", {}).get("tests") == 2 and e.get("exit_status") == 0, json.dumps(e.get("counts")))


# ---------------------------------------------------------------- AT-10 (EC-21 remote/stack)

def _synth_repo():
    proj = tempfile.mkdtemp()
    g = lambda *a: subprocess.run(["git", *a], cwd=proj, capture_output=True, text=True)
    g("init", "-q"); g("config", "user.email", "t@t"); g("config", "user.name", "t")
    open(os.path.join(proj, "f.txt"), "w").write("base\n")
    g("add", "."); g("commit", "-qm", "base")
    base = rev(["rev-parse", "HEAD"], cwd=proj)
    g("branch", "main")   # stack ref for main-line
    open(os.path.join(proj, "f.txt"), "w").write("child\n")
    g("add", "."); g("commit", "-qm", "child")
    head = rev(["rev-parse", "HEAD"], cwd=proj)
    return proj, base, head


def test_at10_declared_base_remote_stack():
    proj, base, head = _synth_repo()
    refs = tmpfile(json.dumps({"main-line": "main", "science-line": "main", "policy-line": "main"}))
    ok200 = tmpfile('{"status": 200, "body": {"sha": "%s"}}' % base)
    # 1. valid base: local + remote 200 exact + ancestral + correct stack
    r = run("check_declared_base.py", "--base", base, "--head", head, "--stack-parent", "main-line",
            "--repo-slug", "x/y", "--api-json", ok200, "--stack-refs", refs, cwd=proj)
    check("AT-10", "valid base passes (local+remote+ancestry+stack)",
          r.returncode == 0, r.stdout)
    # 2. locally present but remote 404 -> NOT_ON_REMOTE (local is not remote proof)
    r = run("check_declared_base.py", "--base", base, "--head", head, "--stack-parent", "main-line",
            "--repo-slug", "x/y", "--api-json", '{"status": 404}', "--stack-refs", refs, cwd=proj)
    check("AT-10", "local object + remote 404 => NOT_ON_REMOTE",
          r.returncode == 1 and "NOT_ON_REMOTE" in r.stdout, r.stdout)
    # 3. remote 403 => VERIFICATION_FAILURE (never misreported as absent)
    r = run("check_declared_base.py", "--base", base, "--head", head, "--stack-parent", "main-line",
            "--repo-slug", "x/y", "--api-json", '{"status": 403}', "--stack-refs", refs, cwd=proj)
    check("AT-10", "403 => REMOTE_VERIFICATION_FAILURE (fail closed)",
          r.returncode == 1 and "REMOTE_VERIFICATION_FAILURE" in r.stdout and "NOT_ON_REMOTE" not in r.stdout, r.stdout)
    # 4. returned SHA mismatch => error
    r = run("check_declared_base.py", "--base", base, "--head", head, "--stack-parent", "main-line",
            "--repo-slug", "x/y", "--api-json", '{"status": 200, "body": {"sha": "%s"}}' % ("f" * 40),
            "--stack-refs", refs, cwd=proj)
    check("AT-10", "returned SHA mismatch fails closed",
          r.returncode == 1 and "does not match requested" in r.stdout, r.stdout)
    # 5. wrong stack-parent (declared science-line; synthetic repo has no science-unique history)
    r = run("check_declared_base.py", "--base", base, "--head", head, "--stack-parent", "science-line",
            "--repo-slug", "x/y", "--api-json", ok200, "--stack-refs", refs, cwd=proj)
    check("AT-10", "wrong stack-parent fails (STACK_MISMATCH)",
          r.returncode == 1 and "STACK_MISMATCH" in r.stdout, r.stdout)
    # 6. malformed SHA
    r = run("check_declared_base.py", "--base", "abc123", "--head", head, "--stack-parent", "main-line",
            "--repo-slug", "x/y", "--api-json", ok200, "--stack-refs", refs, cwd=proj)
    check("AT-10", "malformed SHA fails (MALFORMED_BASE)",
          r.returncode == 1 and "MALFORMED_BASE" in r.stdout, r.stdout)
    # 7. incomplete history (shallow) fails closed
    open(os.path.join(proj, ".git", "shallow"), "w").write("")
    r = run("check_declared_base.py", "--base", base, "--head", head, "--stack-parent", "main-line",
            "--repo-slug", "x/y", "--api-json", ok200, "--stack-refs", refs, cwd=proj)
    check("AT-10", "shallow repo fails (INCOMPLETE_HISTORY)",
          r.returncode == 1 and "INCOMPLETE_HISTORY" in r.stdout, r.stdout)
    # 8. unrelated but locally-present base (fixture lies about remote to isolate ancestry)
    proj2, base2, head2 = _synth_repo()
    proj3, base3, _ = _synth_repo()
    refs2 = tmpfile(json.dumps({"main-line": "main", "science-line": "main", "policy-line": "main"}))
    ok3 = tmpfile('{"status": 200, "body": {"sha": "%s"}}' % base3)
    r = run("check_declared_base.py", "--base", base3, "--head", "HEAD", "--stack-parent", "main-line",
            "--repo-slug", "x/y", "--api-json", ok3, "--stack-refs", refs2, cwd=proj2)
    check("AT-10", "unrelated base fails ancestry",
          r.returncode == 1 and "ANCESTRY" in r.stdout, r.stdout)


# ---------------------------------------------------------------- AT-11 (EC-22 manifest)

def test_at11_status_manifest():
    od = "c" * 40
    h1 = "5" * 40
    manifest = tmpfile(json.dumps({"runs": {
        "39376279686": {"head_sha": od, "conclusion": "success"},
        "39376518935": {"head_sha": h1, "conclusion": "failure"}}}))
    bad = tmpfile(
        f"STATUS: RUN_PASS(run_id=99999999999, sha={od})\n"          # fabricated run id
        f"STATUS: RUN_PASS(run_id=39376279686, sha={h1})\n"          # wrong SHA
        f"STATUS: RUN_PASS(run_id=39376518935, sha={h1})\n"          # non-success conclusion
        "STATUS: WEIRD\n"                                            # unknown token
        "Progress: 92% (run 39376279686)\n"                          # unsupported pct
        "Gitleaks came back clean scan\n"                            # uncited clean scan
        "Failures are pre-existing\n")                               # unsupported pre-existing
    r = run("status_vocabulary_lint.py", bad, "--runs-manifest", manifest)
    check("AT-11", "fabricated run ID flagged", "not found in the runs manifest" in r.stdout, r.stdout)
    check("AT-11", "wrong-SHA run flagged", "not the cited sha" in r.stdout, r.stdout)
    check("AT-11", "non-success run flagged", "not success" in r.stdout, r.stdout)
    check("AT-11", "unknown token flagged", "WEIRD" in r.stdout, r.stdout)
    check("AT-11", "unsupported percentage flagged", "not supported by run records" in r.stdout, r.stdout)
    check("AT-11", "uncited clean scan flagged", "'clean scan'" in r.stdout, r.stdout)
    check("AT-11", "unsupported pre-existing flagged", "'pre-existing'" in r.stdout, r.stdout)
    check("AT-11", "any finding -> exit 1", r.returncode == 1, f"exit={r.returncode}")
    good = tmpfile(
        f"STATUS: NOT_RUN\nSTATUS: RUN_PASS(run_id=39376279686, sha={od})\n"
        f"Scan claim: clean scan (scan=gitleaks run 39376279686, sha={od})\n"
        f"Failure is pre-existing (base={'2' * 40})\n"
        f"STATUS: RUN_PASS(run_id=39376279686, sha={od}) pct=50\n")
    r = run("status_vocabulary_lint.py", good, "--runs-manifest", manifest)
    check("AT-11", "manifest-backed statuses pass with computed pct", r.returncode == 0, r.stdout)
    # without a manifest, RUN_PASS must fail closed
    r = run("status_vocabulary_lint.py", good)
    check("AT-11", "RUN_PASS without manifest fails closed",
          r.returncode == 1 and "no unverified progress" in r.stdout, r.stdout)


# ---------------------------------------------------------------- AT-12 (EC-16 strictness)

def test_at12_strict_sha_lists():
    head = rev(["rev-parse", "HEAD"])
    # 1. one malformed entry rejects the WHOLE list (the valid SHA is not checked)
    mixed = tmpfile(f"Claimed-SHAs:\n  {head}\n  abc123\n  {'a' * 41}\n")
    r = run("verify_claimed_shas.py", "--body", mixed, "--repo", "x/y",
            "--api-json", '{"status": 200, "body": {"sha": "%s"}}' % head)
    check("AT-12", "any malformed entry rejects the whole list",
          r.returncode == 1 and "whole list is rejected" in r.stdout and "abc123" in r.stdout, r.stdout)
    check("AT-12", "no partial verification happened (valid SHA not reported OK)",
          "verified on remote" not in r.stdout, r.stdout)
    # 2. remote 404 = absent (distinct)
    body404 = tmpfile(f"Claimed-SHAs:\n  {head}\n")
    r = run("verify_claimed_shas.py", "--body", body404, "--repo", "x/y", "--api-json", '{"status": 404}')
    check("AT-12", "404 reported as not on remote",
          r.returncode == 1 and "not on remote" in r.stdout, r.stdout)
    # 3. remote 403/429/500 = VERIFICATION_FAILURE (never misreported as absent)
    for code in (401, 403, 429, 500):
        r = run("verify_claimed_shas.py", "--body", body404, "--repo", "x/y",
                "--api-json", '{"status": %d}' % code)
        check("AT-12", f"HTTP {code} => VERIFICATION_FAILURE (fail closed)",
              r.returncode == 1 and "VERIFICATION_FAILURE" in r.stdout and "not on remote" not in r.stdout, r.stdout)
    # 4. returned SHA mismatch = fail closed
    r = run("verify_claimed_shas.py", "--body", body404, "--repo", "x/y",
            "--api-json", '{"status": 200, "body": {"sha": "%s"}}' % ("e" * 40))
    check("AT-12", "returned-SHA mismatch fails closed",
          r.returncode == 1 and "does not match requested" in r.stdout, r.stdout)
    # 5. short and uppercase entries rejected by name
    ugly = tmpfile("Claimed-SHAs:\n  " + head[:12] + "\n  " + head.upper() + "\n")
    r = run("verify_claimed_shas.py", "--body", ugly, "--repo", "x/y",
            "--api-json", '{"status": 200, "body": {"sha": "%s"}}' % head)
    check("AT-12", "short and malformed entries named in the rejection",
          r.returncode == 1 and r.stdout.count("malformed Claimed-SHAs entry") >= 2, r.stdout)


# ---------------------------------------------------------------- AT-13 (EC-18/24 append-only)

def test_at13_append_only_pushlog():
    older_path = os.path.join(REPO, "docs", "evidence", "push-log", "PUSH_LOG.md")
    older = open(older_path, "rb").read()
    # 1. real current log vs 17eb881 version: append-only holds
    r = run("evidence_sequence_check.py", older_path, "--parent-prefix",
            tmpfile(older.decode("utf-8"), ".md"))
    check("AT-13", "real PUSH_LOG preserves the 17eb881 version as prefix", r.returncode == 0, r.stdout)
    # 2. synthetic: pure append passes
    d = tempfile.mkdtemp()
    older_f = os.path.join(d, "old.md"); cur_f = os.path.join(d, "cur.md")
    open(older_f, "w").write("## 2026-10-20T06:00:00Z | E1 | a | b\n")
    open(cur_f, "w").write("## 2026-10-20T06:00:00Z | E1 | a | b\n## 2026-10-20T07:00:00Z | E2 | a | c\n")
    r = run("evidence_sequence_check.py", cur_f, "--parent-prefix", older_f)
    check("AT-13", "pure append passes parent-prefix", r.returncode == 0, r.stdout)
    # 3. edited prior entry fails
    open(cur_f, "w").write("## 2026-10-20T06:00:00Z | E1 | a | EDITED\n## 2026-10-20T07:00:00Z | E2 | a | c\n")
    r = run("evidence_sequence_check.py", cur_f, "--parent-prefix", older_f)
    check("AT-13", "edited prior entry fails", r.returncode == 1 and "MODIFIED" in r.stdout, r.stdout)
    # 4. deleted prior entry fails
    open(cur_f, "w").write("## 2026-10-20T07:00:00Z | E2 | a | c\n")
    r = run("evidence_sequence_check.py", cur_f, "--parent-prefix", older_f)
    check("AT-13", "deleted prior entry fails", r.returncode == 1, r.stdout)
    # 5. pushlog truthfulness: placeholder CONFIRM rejected; real SHA accepted
    ph = tmpfile("## 2026-10-20T06:06:00Z | P-011 | push | branch=x sha=self-referential-tip credential=OWNER-PAT-1\n"
                 "## 2026-10-20T09:00:00Z | P-011-CONFIRM | push-confirm | branch=x remote_sha=17eb88184c9f461c6bc9ee592f916f7ff4c921c4 credential=OWNER-PAT-1 evidence=\"git ls-remote origin refs/heads/x\"\n")
    r = run("evidence_sequence_check.py", ph, "--pushlog-schema")
    check("AT-13", "legacy placeholder push corrected by a well-formed CONFIRM passes (the P-011 case)",
          r.returncode == 0, r.stdout)
    ph2 = tmpfile("## 2026-10-20T09:00:00Z | P-9-CONFIRM | push-confirm | branch=x remote_sha=self-referential-tip credential=OWNER-PAT-1 evidence=\"git ls-remote origin refs/heads/x\"\n")
    r = run("evidence_sequence_check.py", ph2, "--pushlog-schema")
    check("AT-13", "CONFIRM with placeholder remote_sha fails truthfulness",
          r.returncode == 1 and "truthfulness" in r.stdout, r.stdout)
    # 6. the real P-011-CONFIRM entry carries the real 40-hex remote SHA
    log = open(older_path, encoding="utf-8").read()
    m = re.search(r"P-011-CONFIRM.*remote_sha=([0-9a-f]{40})", log)
    check("AT-13", "P-011-CONFIRM carries the real full SHA 17eb881...",
          bool(m) and m.group(1) == "17eb88184c9f461c6bc9ee592f916f7ff4c921c4", "regex miss")
    r = run("evidence_sequence_check.py", older_path, "--pushlog-schema")
    check("AT-13", "real PUSH_LOG passes --pushlog-schema", r.returncode == 0, r.stdout)


# ---------------------------------------------------------------- AT-14 (EC-25 remediation)

def test_at14_deploy_remediation():
    fx = os.path.join(REPO, "proposals", "remediation", "fixtures")
    rem = os.path.join(REPO, "proposals", "remediation")
    d = tempfile.mkdtemp()
    for fn in ("deploy_on_push_VIOLATION.yml", "toscanini_gate_AIN_DEFECT.yml"):
        open(os.path.join(d, fn), "w").write(open(os.path.join(fx, fn)).read())
    r = run("workflow_deploy_trigger_lint.py", d)
    check("AT-14", "main's push-triggered deploy fixture is FLAGGED",
          r.returncode == 1 and "deploy_on_push_VIOLATION.yml" in r.stdout, r.stdout)
    d2 = tempfile.mkdtemp()
    open(os.path.join(d2, "deploy.dispatch-only.yml"), "w").write(
        open(os.path.join(rem, "deploy.dispatch-only.yml")).read())
    r = run("workflow_deploy_trigger_lint.py", d2)
    check("AT-14", "dispatch-only remediated deploy is ACCEPTED by the lint", r.returncode == 0, r.stdout)
    import yaml
    g = yaml.safe_load(open(os.path.join(rem, "toscanini-gate.branch-triggers.yml")))
    on = g.get(True, g.get("on"))
    check("AT-14", "proposed gate push filter matches nested branches",
          on["push"]["branches"] == ["**"], str(on["push"]))
    check("AT-14", "proposed gate pull_request filter repaired to [main]",
          on["pull_request"]["branches"] == ["main"], str(on["pull_request"]))
    body_lines = [l for l in open(os.path.join(rem, "toscanini-gate.branch-triggers.yml"))
                  if not l.lstrip().startswith("#")]
    check("AT-14", "proposed gate body contains no 'ain]' defect",
          not any("ain]" in l for l in body_lines), "found in non-comment lines")
    fxg = yaml.safe_load(open(os.path.join(fx, "toscanini_gate_AIN_DEFECT.yml")))
    fx_on = fxg.get(True, fxg.get("on"))
    check("AT-14", "fixture reproduces main's defects (push '*', pull_request 'ain]')",
          fx_on["push"]["branches"] == ["*"] and fx_on["pull_request"]["branches"] == ["ain]"],
          str(fx_on))
    check("AT-14", "remediation README states the lint-activation precondition",
          "PRECONDITION" in open(os.path.join(REPO, "proposals", "ci", "README.md")).read())


# ---------------------------------------------------------------- AT-15 (EC-23 structure)

def test_at15_gitleaks_pinned_full_coverage():
    y = open(os.path.join(REPO, "proposals", "ci", "gitleaks-all-branches.yml")).read()
    m = re.search(r"gitleaks/gitleaks-action@([0-9a-f]+)", y)
    check("AT-15", "action pinned to a full 40-hex commit (never a tag)",
          bool(m) and len(m.group(1)) == 40, "no 40-hex pin found")
    check("AT-15", "SARIF upload disabled (least privilege, no security-events needed)",
          'GITLEAKS_ENABLE_UPLOAD_ARTIFACT: "false"' in y)
    check("AT-15", "permissions: contents: read declared",
          re.search(r"permissions:\s*\n\s*contents: read", y) is not None)
    check("AT-15", "all-refs job scans --log-opts=--all (evidence branches included)",
          "gitleaks-all-refs" in y and '--log-opts="--all"' in y)
    check("AT-15", "binary download verified against the official checksum",
          "sha256sum -c -" in y and "551f6fc83ea457d62a0d98237cbad105af8d557003051f41f3e7ca7b3f2470eb" in y)
    check("AT-15", "seeded-secret positive + clean negative controls present",
          "seeded-secret-detection-test" in y and "NOT detected" in y and "clean repo exits 0" in y)
    check("AT-15", "shallow checkout fails closed before any scan",
          "is-shallow-repository" in y)
    import yaml as _y
    doc = _y.safe_load(y)
    check("AT-15", "three jobs defined (checkout scan, all-refs scan, detection test)",
          set(doc["jobs"].keys()) == {"gitleaks", "gitleaks-all-refs", "seeded-secret-detection-test"},
          str(set(doc["jobs"].keys())))


# ---------------------------------------------------------------- AT-16 (inventory counts)

def test_at16_inventory_counts():
    ci_dir = os.path.join(REPO, "proposals", "ci")
    ymls = sorted(f for f in os.listdir(ci_dir) if f.endswith((".yml", ".yaml")))
    check("AT-16", "seven workflow YAML drafts exist", len(ymls) == 7, str(ymls))
    check("AT-16", "expected seven names", set(ymls) == {
        "check-declared-base.yml", "evidence-sequence-check.yml", "gitleaks-all-branches.yml",
        "no-deploy-on-push-lint.yml", "reading-note-map-check.yml", "status-vocabulary-lint.yml",
        "verify-claimed-shas.yml"}, str(ymls))
    checkers = sorted(f for f in os.listdir(CI)
                      if f.endswith(".py") and os.path.isfile(os.path.join(CI, f)))
    check("AT-16", "six checker scripts exist", len(checkers) == 6, str(checkers))
    reg = json.load(open(os.path.join(HERE, "ACCEPTANCE_TESTS.json")))["tests"]
    idx = open(os.path.join(REPO, "docs", "directives", "INDEX.md")).read()
    idx_ats = set(re.findall(r"AT-\d+", idx))
    check("AT-16", "every AT named in INDEX exists in the registry",
          idx_ats <= set(reg.keys()), str(idx_ats - set(reg.keys())))
    src = open(os.path.join(HERE, "test_acceptance.py")).read()
    fns = set(re.findall(r"def (test_at\d+_\w+)", src))
    nodes = {v["node"].split("::")[-1] for v in reg.values()}
    check("AT-16", "registry node IDs all exist in the test source",
          nodes <= fns, str(nodes - fns))
    receipt = open(os.path.join(REPO, "docs", "evidence", "ec16-ec26", "RECEIPT.md")).read()
    check("AT-16", "receipt carries the corrected counts (7 workflows / 6 checkers)",
          "7 workflow" in receipt and "6 checker" in receipt, "correction entry missing")


# ------------------------------------------------- AT-18 (EC-20 scoped collection, Round-19)

def test_at18_scoped_collection():
    """Regression for the EC-20 activation failures (runs 39417485854/88957/92355).

    Root cause: repo-wide `pytest --collect-only` imports EVERY test module under
    the root; legacy script-style modules that sys.exit() at import time raise
    SystemExit inside pytest's collector -> INTERNALERROR -> collection aborts ->
    named-but-actually-present tests were reported "not collectable". The fix
    scopes collection to the files that contain the named tests.
    """
    proj = tempfile.mkdtemp()
    subprocess.run(["git", "init", "-q"], cwd=proj)
    subprocess.run(["git", "config", "user.email", "t@t"], cwd=proj)
    subprocess.run(["git", "config", "user.name", "t"], cwd=proj)
    os.makedirs(os.path.join(proj, "legacytests"))
    # a poison module: executes at import and exits nonzero (gpu_worker-style)
    open(os.path.join(proj, "legacytests", "test_poison.py"), "w").write(
        "import sys\nprint('legacy script suite running')\nsys.exit(1)\n")
    # a good pytest module alongside it
    open(os.path.join(proj, "legacytests", "test_good.py"), "w").write(
        "def test_named_passes():\n    assert True\n")
    subprocess.run(["git", "add", "."], cwd=proj, capture_output=True)
    subprocess.run(["git", "commit", "-qm", "c"], cwd=proj, capture_output=True)
    ev = os.path.join(proj, "evidence.json")
    # 1. a poison module elsewhere in the tree must NOT abort collection of the
    #    named test (this is the exact regression: it used to fail the EC-20 job)
    body = tmpfile("Acceptance-Tests:\n  legacytests/test_good.py::test_named_passes\n")
    r = run("verify_claimed_shas.py", "--body", body, "--check-tests", "--evidence-out", ev, cwd=proj)
    check("AT-18", "unrelated import-time-sys.exit module does not abort named-test collection",
          r.returncode == 0, r.stdout)
    e = json.load(open(ev))
    check("AT-18", "scoped collection executed the named test",
          e.get("counts", {}).get("tests") == 1 and e.get("exit_status") == 0, json.dumps(e.get("counts")))
    check("AT-18", "evidence records the per-file collection scope",
          e.get("collection_scope_files") == ["legacytests/test_good.py"], json.dumps(e.get("collection_scope_files")))
    # 2. naming a test in the poison module still fails closed (not silently skipped)
    body = tmpfile("Acceptance-Tests:\n  legacytests/test_poison.py::test_anything\n")
    r = run("verify_claimed_shas.py", "--body", body, "--check-tests", "--evidence-out", ev, cwd=proj)
    check("AT-18", "test inside an unimportable module fails closed as not collectable",
          r.returncode == 1 and "not collectable" in r.stdout, r.stdout)
    # 3. naming a file that does not exist fails closed with its own message
    body = tmpfile("Acceptance-Tests:\n  no/such/file.py::test_x\n")
    r = run("verify_claimed_shas.py", "--body", body, "--check-tests", "--evidence-out", ev, cwd=proj)
    check("AT-18", "named test file that does not exist fails closed",
          r.returncode == 1 and "does not exist" in r.stdout, r.stdout)


# --------------------------------------- AT-17 (EC-26 ACTIVE-liveness, Round-19)

def test_at17_active_workflow_liveness():
    """ACTIVE rules must map to LIVE workflows, not inert proposals drafts.

    Round-19 direction 6 (EC-26): the checker previously verified only that a
    referenced file/job EXISTS, so an ACTIVE rule citing
    workflow:proposals/ci/... passed even though that workflow is inert.
    """
    import yaml as _yaml
    d = tempfile.mkdtemp()
    # live workflow with the job an ACTIVE rule may cite
    os.makedirs(os.path.join(d, ".github", "workflows"), exist_ok=True)
    open(os.path.join(d, ".github", "workflows", "gate.yml"), "w").write(
        "name: gate\non:\n  push:\njobs:\n  secret-scan:\n    name: Secret Scan\n"
        "    runs-on: ubuntu-latest\n")
    # inert draft under proposals/ (what an ACTIVE rule must NOT cite)
    os.makedirs(os.path.join(d, "proposals", "ci"), exist_ok=True)
    open(os.path.join(d, "proposals", "ci", "draft.yml"), "w").write(
        "name: draft\non:\n  push:\njobs:\n  claimed-shas:\n    runs-on: ubuntu-latest\n")
    reg = {"tests": {}}
    open(os.path.join(d, "registry.json"), "w").write(json.dumps(reg))
    os.makedirs(os.path.join(d, "docs", "directives"), exist_ok=True)

    def idx(rows):
        body = "\n".join(f"| {rid} | {st} | {enf} | NONE |" for rid, st, enf in rows)
        text = ( "# index\n<!-- BEGIN RULE INVENTORY -->\n"
                 "reading-note-path: docs/note.md\n"
                 "pinned-blobs: b6bb35e84956d84a8c9d9e65734c212baf3ed075\n"
                 "| Rule ID | Status | Enforcement point(s) | Acceptance test(s) |\n"
                 "|---|---|---|---|\n" + body + "\n<!-- END RULE INVENTORY -->\n" )
        open(os.path.join(d, "docs", "directives", "INDEX.md"), "w").write(text)
        note_lines = ["blob b6bb35e84956d84a8c9d9e65734c212baf3ed075 read"]
        note_lines += [f"| {rid} | -> mapped |" for rid, _, _ in rows]
        open(os.path.join(d, "docs", "note.md"), "w").write("\n".join(note_lines) + "\n")

    def run_idx():
        return run("reading_note_map_check.py", "--index", "docs/directives/INDEX.md",
                   "--repo-root", d, "--registry", "registry.json", "--include-proposed")

    # 1. ACTIVE rule citing an inert proposals/ draft MUST FAIL
    idx([("EC-X", "ACTIVE", "workflow:proposals/ci/draft.yml#claimed-shas")])
    r = run_idx()
    check("AT-17", "ACTIVE rule citing an inert proposals/ workflow fails",
          r.returncode == 1 and "outside the live CI path" in r.stdout, r.stdout)
    # 2. ACTIVE rule citing a ref-pinned workflow MUST FAIL (pin can diverge from live)
    idx([("EC-X", "ACTIVE", "workflow:.github/workflows/gate.yml@origin/main#secret-scan")])
    r = run_idx()
    check("AT-17", "ACTIVE rule citing an @ref-pinned workflow fails",
          r.returncode == 1 and "pins its workflow" in r.stdout, r.stdout)
    # 3. ACTIVE rule citing a live working-tree workflow with the named job PASSES
    idx([("EC-X", "ACTIVE", "workflow:.github/workflows/gate.yml#secret-scan")])
    r = run_idx()
    check("AT-17", "ACTIVE rule citing the live unpinned workflow passes",
          r.returncode == 0, r.stdout)
    # 4. PROPOSED rule citing an inert draft is still ALLOWED (activation path)
    idx([("EC-X", "PROPOSED", "workflow:proposals/ci/draft.yml#claimed-shas")])
    r = run_idx()
    check("AT-17", "PROPOSED rule citing an inert proposals/ workflow is allowed",
          r.returncode == 0, r.stdout)
    # 5. ACTIVE rule citing a live-path workflow that is missing from the tree fails
    idx([("EC-X", "ACTIVE", "workflow:.github/workflows/not-there.yml#some-job")])
    r = run_idx()
    check("AT-17", "ACTIVE rule citing a live-path file absent from the tree fails",
          r.returncode == 1 and "absent from the working tree" in r.stdout, r.stdout)


# ------------------------------------ AT-20 (EC-18/24 schema hardening, Round-19)

def test_at20_pushlog_schema_hardening():
    """Round-19 direction 4: the validation must reject duplicate IDs, malformed
    confirmations, missing required fields, and invalid parent-prefix history;
    an unpaired INTENT is reported PENDING without implying completion."""
    sha = lambda c: c * 40
    # 1. duplicate ENTRY-ID across a merged history (the P-015-CONFIRM collision class)
    dup = tmpfile("## 2026-10-20T10:00:00Z | P-015-INTENT | push-intent | branch=b credential=OWNER-PAT-1\n"
                  "## 2026-10-20T11:00:00Z | P-015-CONFIRM | push-confirm | branch=b remote_sha=" + sha("1") + " credential=OWNER-PAT-1 evidence=\"ls-remote b\"\n"
                  "## 2026-10-20T12:00:00Z | P-015-CONFIRM | push-confirm | branch=c remote_sha=" + sha("2") + " credential=OWNER-PAT-1 evidence=\"ls-remote c\"\n")
    r = run("evidence_sequence_check.py", dup, "--pushlog-schema")
    check("AT-20", "duplicate ENTRY-ID rejected", r.returncode == 1 and "duplicate ENTRY-ID P-015-CONFIRM" in r.stdout, r.stdout)
    # 2. CONFIRM missing a required field (evidence=) fails
    miss = tmpfile("## 2026-10-20T10:00:00Z | P-1-INTENT | push-intent | branch=b credential=OWNER-PAT-1\n"
                   "## 2026-10-20T11:00:00Z | P-1-CONFIRM | push-confirm | branch=b remote_sha=" + sha("1") + " credential=OWNER-PAT-1\n")
    r = run("evidence_sequence_check.py", miss, "--pushlog-schema")
    check("AT-20", "CONFIRM missing evidence= field fails", r.returncode == 1 and "required field 'evidence'" in r.stdout, r.stdout)
    # 3. CONFIRM without any preceding INTENT fails
    lone = tmpfile("## 2026-10-20T11:00:00Z | P-2-CONFIRM | push-confirm | branch=b remote_sha=" + sha("3") + " credential=OWNER-PAT-1 evidence=\"ls-remote\"\n")
    r = run("evidence_sequence_check.py", lone, "--pushlog-schema")
    check("AT-20", "CONFIRM without a preceding INTENT fails", r.returncode == 1 and "no preceding" in r.stdout, r.stdout)
    # 4. INTENT without CONFIRM is reported PENDING and does not fail
    pend = tmpfile("## 2026-10-20T10:00:00Z | P-3-INTENT | push-intent | branch=b credential=OWNER-PAT-1\n")
    r = run("evidence_sequence_check.py", pend, "--pushlog-schema")
    check("AT-20", "unpaired INTENT reported PENDING without failing",
          r.returncode == 0 and "PENDING" in r.stdout and "P-3-INTENT" in r.stdout, r.stdout)
    # 5. INTENT missing branch= fails
    bad = tmpfile("## 2026-10-20T10:00:00Z | P-4-INTENT | push-intent | credential=OWNER-PAT-1\n")
    r = run("evidence_sequence_check.py", bad, "--pushlog-schema")
    check("AT-20", "INTENT missing branch= fails", r.returncode == 1 and "required field 'branch'" in r.stdout, r.stdout)
    # 6. ENTRY-ID outside the P- grammar rejected in schema mode
    weird = tmpfile("## 2026-10-20T10:00:00Z | E1 | push-intent | branch=b credential=OWNER-PAT-1\n")
    r = run("evidence_sequence_check.py", weird, "--pushlog-schema")
    check("AT-20", "ENTRY-ID outside the P- grammar rejected", r.returncode == 1 and "grammar P-<slug>" in r.stdout, r.stdout)
    # 7. the real mainline push log passes the hardened schema
    real = os.path.join(REPO, "docs", "evidence", "push-log", "PUSH_LOG.md")
    r = run("evidence_sequence_check.py", real, "--pushlog-schema")
    check("AT-20", "the real mainline push log passes the hardened schema", r.returncode == 0, r.stdout)


# ------------------------------------ AT-19 (EC-22 Status-summary, Round-19)

def _empty_manifest():
    f = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
    f.write('{"runs": {}}')
    f.close()
    return f.name


def test_at19_status_summary_grammar():
    """Round-19: PR-body Status-summary lines carry the EC-22 vocabulary.

    The PR template declares the vocabulary for Status-summary; the lint now
    enforces the same grammar as STATUS: lines, so off-vocabulary tokens
    (e.g. ACTIVATION_PENDING) in a PR body or report fail machine-side.
    """
    # 1. off-vocabulary token in a Status-summary line is flagged
    bad = tmpfile("Status-summary: ACTIVATION_PENDING(branch=x, sha=" + "1" * 40 + ")\n")
    r = run("status_vocabulary_lint.py", bad, "--runs-manifest", _empty_manifest())
    check("AT-19", "off-vocabulary Status-summary token flagged",
          r.returncode == 1 and "ACTIVATION_PENDING" in r.stdout, r.stdout)
    # 2. vocabulary-compliant Status-summary passes
    good = tmpfile("Status-summary: NOT_RUN (activation workflow not yet rerun on the fixed candidate)\n")
    r = run("status_vocabulary_lint.py", good, "--runs-manifest", _empty_manifest())
    check("AT-19", "compliant Status-summary passes", r.returncode == 0, r.stdout)
    # 3. RUN_PASS inside Status-summary needs run_id + sha and a manifest hit
    rp = tmpfile("Status-summary: RUN_PASS (run_id=123, sha=" + "a" * 40 + ")\n")
    r = run("status_vocabulary_lint.py", rp, "--runs-manifest", _empty_manifest())
    check("AT-19", "RUN_PASS in Status-summary is validated against the runs manifest",
          r.returncode == 1 and "not found in the runs manifest" in r.stdout, r.stdout)
    # 4. PUSHED / MERGED tokens pass without a manifest record
    pm = tmpfile("Status-summary: PUSHED\nStatus-summary: MERGED\n")
    r = run("status_vocabulary_lint.py", pm, "--runs-manifest", _empty_manifest())
    check("AT-19", "PUSHED/MERGED Status-summary tokens pass", r.returncode == 0, r.stdout)
