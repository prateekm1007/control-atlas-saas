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
    # --runs-manifest takes a FILE path (the lint fail-closes on an unusable
    # manifest, so the inline-JSON form this test once passed was itself a
    # hidden check() failure -- the lint correctly refused it).
    empty_manifest = tmpfile('{"runs": {}}', ".json")
    bad = tmpfile("STATUS: RUNNING\nSTATUS: DONE\nSTATUS: RUN_PASS\n"
                  "Gitleaks came back clean scan\nFailures are pre-existing\nProgress: 92%\n")
    r = run("status_vocabulary_lint.py", bad, "--runs-manifest", empty_manifest)
    check("AT-04", "unknown status tokens flagged", "RUNNING" in r.stdout and "DONE" in r.stdout, r.stdout)
    check("AT-04", "RUN_PASS without id/args flagged", "without run_id" in r.stdout, r.stdout)
    check("AT-04", "clean scan without citation flagged", "'clean scan'" in r.stdout, r.stdout)
    check("AT-04", "pre-existing without citation flagged", "'pre-existing'" in r.stdout, r.stdout)
    check("AT-04", "percentage without run records flagged", "not supported by run records" in r.stdout, r.stdout)
    good = tmpfile("STATUS: NOT_RUN\nSTATUS: RUN_FAIL\nSTATUS: PUSHED\nSTATUS: MERGED\n")
    r = run("status_vocabulary_lint.py", good, "--runs-manifest", empty_manifest)
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
    # 3. passing named tests produce a complete evidence record (the body must
    #    also carry the fail-closed EC-16/EC-19 fields: Claimed-SHAs 'none' and a
    #    Directive-ID -- without them the checker correctly exits 1)
    body = tmpfile("Directive-ID: TEST\nClaimed-SHAs: none\n"
                   "Acceptance-Tests:\n  test_ok.py::test_first_passes\n  test_ok.py::test_second_passes\n")
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
    # each repo salts its file content: two synth repos created within the same
    # second would otherwise carry identical commit timestamps and produce the
    # SAME commit SHA (the timestamp-collision class the Round-19 fix removed
    # and the R20/R21 rewrites reintroduced -- it made the "unrelated base"
    # scenario pass because base3 literally was base2)
    salt = os.urandom(8).hex()
    proj = tempfile.mkdtemp()
    g = lambda *a: subprocess.run(["git", *a], cwd=proj, capture_output=True, text=True)
    g("init", "-q"); g("config", "user.email", "t@t"); g("config", "user.name", "t")
    open(os.path.join(proj, "f.txt"), "w").write(f"base {salt}\n")
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
    nf404 = tmpfile('{"status": 404}')
    fd403 = tmpfile('{"status": 403}')
    fmism = tmpfile('{"status": 200, "body": {"sha": "%s"}}' % ("f" * 40))
    # 1. valid base: local + remote 200 exact + ancestral + correct stack
    r = run("check_declared_base.py", "--base", base, "--head", head, "--stack-parent", "main-line",
            "--repo-slug", "x/y", "--api-json", ok200, "--stack-refs", refs, cwd=proj)
    check("AT-10", "valid base passes (local+remote+ancestry+stack)",
          r.returncode == 0, r.stdout)
    # 2. locally present but remote 404 -> NOT_ON_REMOTE (local is not remote proof)
    r = run("check_declared_base.py", "--base", base, "--head", head, "--stack-parent", "main-line",
            "--repo-slug", "x/y", "--api-json", nf404, "--stack-refs", refs, cwd=proj)
    check("AT-10", "local object + remote 404 => NOT_ON_REMOTE",
          r.returncode == 1 and "NOT_ON_REMOTE" in r.stdout, r.stdout)
    # 3. remote 403 => VERIFICATION_FAILURE (never misreported as absent)
    r = run("check_declared_base.py", "--base", base, "--head", head, "--stack-parent", "main-line",
            "--repo-slug", "x/y", "--api-json", fd403, "--stack-refs", refs, cwd=proj)
    check("AT-10", "403 => REMOTE_VERIFICATION_FAILURE (fail closed)",
          r.returncode == 1 and "REMOTE_VERIFICATION_FAILURE" in r.stdout and "NOT_ON_REMOTE" not in r.stdout, r.stdout)
    # 4. returned SHA mismatch => error
    r = run("check_declared_base.py", "--base", base, "--head", head, "--stack-parent", "main-line",
            "--repo-slug", "x/y", "--api-json", fmism,
            "--stack-refs", refs, cwd=proj)
    check("AT-10", "returned SHA mismatch fails closed",
          r.returncode == 1 and "does not match requested" in r.stdout, r.stdout)
    # 5. wrong stack-parent: an ORPHAN science line whose history does not
    #    contain the base (the previous fixture mapped every stack class to the
    #    same ref, so the base was trivially valid and no mismatch could occur)
    gs = lambda *a: subprocess.run(["git", *a], cwd=proj, capture_output=True, text=True)
    gs("checkout", "-q", "--orphan", "science")
    open(os.path.join(proj, "s.txt"), "w").write("science line\n")
    gs("add", "."); gs("commit", "-qm", "science root")
    gs("checkout", "-q", "main")
    refs_sc = tmpfile(json.dumps({"main-line": "main", "science-line": "science", "policy-line": "main"}))
    r = run("check_declared_base.py", "--base", base, "--head", head, "--stack-parent", "science-line",
            "--repo-slug", "x/y", "--api-json", ok200, "--stack-refs", refs_sc, cwd=proj)
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
    # 8. unrelated but locally-present base (fixture lies about remote to isolate
    #    ancestry): the salted repos never share objects, so proj2 must FETCH
    #    proj3's history for the base to be locally present at all -- then the
    #    ancestry check is what must fail it
    proj2, base2, head2 = _synth_repo()
    proj3, base3, _ = _synth_repo()
    subprocess.run(["git", "-C", proj2, "fetch", "-q", proj3], check=True)
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
        f"STATUS: RUN_FAIL\n"
        f"Scan claim: clean scan (scan=gitleaks run 39376279686, sha={od})\n"
        f"Failure is pre-existing (base={'2' * 40})\n"
        f"STATUS: RUN_PASS(run_id=39376279686, sha={od}) pct=50\n")
    # 4 STATUS lines, 2 accepted RUN_PASS -> computed progress 50% == the pct= claim
    r = run("status_vocabulary_lint.py", good, "--runs-manifest", manifest)
    check("AT-11", "manifest-backed statuses pass with computed pct", r.returncode == 0, r.stdout)
    # without a manifest, RUN_PASS must fail closed
    r = run("status_vocabulary_lint.py", good)
    check("AT-11", "RUN_PASS without manifest fails closed",
          r.returncode == 1 and "no unverified progress" in r.stdout, r.stdout)


# ---------------------------------------------------------------- AT-12 (EC-16 strictness)

def test_at12_strict_sha_lists():
    head = rev(["rev-parse", "HEAD"])
    # --api-json takes a FILE path (inline JSON made the checker fail closed on
    # an unusable fixture file -- the hidden-failure class this round fixes)
    ok_head = tmpfile('{"status": 200, "body": {"sha": "%s"}}' % head)
    nf404 = tmpfile('{"status": 404}')
    fmism = tmpfile('{"status": 200, "body": {"sha": "%s"}}' % ("e" * 40))
    codes = {c: tmpfile('{"status": %d}' % c) for c in (401, 403, 429, 500)}
    # 1. one malformed entry rejects the WHOLE list (the valid SHA is not checked)
    mixed = tmpfile(f"Claimed-SHAs:\n  {head}\n  abc123\n  {'a' * 41}\n")
    r = run("verify_claimed_shas.py", "--body", mixed, "--repo", "x/y",
            "--api-json", ok_head)
    check("AT-12", "any malformed entry rejects the whole list",
          r.returncode == 1 and "whole list is rejected" in r.stdout and "abc123" in r.stdout, r.stdout)
    check("AT-12", "no partial verification happened (valid SHA not reported OK)",
          "verified on remote" not in r.stdout, r.stdout)
    # 2. remote 404 = absent (distinct)
    body404 = tmpfile(f"Claimed-SHAs:\n  {head}\n")
    r = run("verify_claimed_shas.py", "--body", body404, "--repo", "x/y", "--api-json", nf404)
    check("AT-12", "404 reported as not on remote",
          r.returncode == 1 and "not on remote" in r.stdout, r.stdout)
    # 3. remote 403/429/500 = VERIFICATION_FAILURE (never misreported as absent)
    for code in (401, 403, 429, 500):
        r = run("verify_claimed_shas.py", "--body", body404, "--repo", "x/y",
                "--api-json", codes[code])
        check("AT-12", f"HTTP {code} => VERIFICATION_FAILURE (fail closed)",
              r.returncode == 1 and "VERIFICATION_FAILURE" in r.stdout and "not on remote" not in r.stdout, r.stdout)
    # 4. returned SHA mismatch = fail closed
    r = run("verify_claimed_shas.py", "--body", body404, "--repo", "x/y",
            "--api-json", fmism)
    check("AT-12", "returned-SHA mismatch fails closed",
          r.returncode == 1 and "does not match requested" in r.stdout, r.stdout)
    # 5. short and uppercase entries rejected by name
    ugly = tmpfile("Claimed-SHAs:\n  " + head[:12] + "\n  " + head.upper() + "\n")
    r = run("verify_claimed_shas.py", "--body", ugly, "--repo", "x/y",
            "--api-json", ok_head)
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
    import re as _re
    body_lines = [l for l in open(os.path.join(rem, "toscanini-gate.branch-triggers.yml"))
                  if not l.lstrip().startswith("#")]
    # match the DEFECT TOKEN 'ain]' precisely: a naive substring also matches the
    # tail of the legitimate token 'main]' (e.g. a flow-style branches: [main]
    # line), which is a false positive, not the defect
    check("AT-14", "proposed gate body contains no 'ain]' defect",
          not any(_re.search(r"(?<!m)ain\]", l) for l in body_lines), "found in non-comment lines")
    fxg = yaml.safe_load(open(os.path.join(fx, "toscanini_gate_AIN_DEFECT.yml")))
    fx_on = fxg.get(True, fxg.get("on"))
    check("AT-14", "fixture reproduces main's defects (push '*', pull_request 'ain]')",
          fx_on["push"]["branches"] == ["*"] and fx_on["pull_request"]["branches"] == "ain]",
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
    check("AT-16", "eleven checker scripts exist",
          len(checkers) == 11, str(checkers))
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
    check("AT-16", "receipt carries the corrected counts (7 workflows / 11 checkers)",
          "7 workflow" in receipt and "11 checker" in receipt, "correction entry missing")


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
    #    named test (this is the exact regression: it used to fail the EC-20 job;
    #    the body carries the fail-closed EC-16/EC-19 fields so a clean run can exit 0)
    body = tmpfile("Directive-ID: TEST\nClaimed-SHAs: none\n"
                   "Acceptance-Tests:\n  legacytests/test_good.py::test_named_passes\n")
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
    #    (the citation must declare the events it needs AND the changes it
    #    governs: the fixture triggers on push and governs the whole tree)
    idx([("EC-X", "ACTIVE", "workflow:.github/workflows/gate.yml#secret-scan~events=push~governs=**")])
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


# --------------------- AT-21 (EC-23 fingerprint allowlist, Round-20)

def _at21_baseline_entry(fp, rule=None, file=None, cls=None):
    # the entry's rule/file/line are DERIVED from the fingerprint itself: a
    # baseline entry that contradicts the fingerprint's own embedded location is
    # exactly the contradiction class the checker must (and does) reject
    parts = fp.split(":")
    return {
        "fingerprint": fp,
        "rule": rule or (parts[2] if len(parts) > 2 else "generic-api-key"),
        "file": file or (parts[1] if len(parts) > 1 else "docs/x.md"),
        "line": int(parts[3]) if len(parts) > 3 and parts[3].isdigit() else 5,
        "commit": parts[0],
        "classification": cls or "FALSE-POSITIVE-PROVENANCE-DIGEST",
        "source": "test scan fixture",
        "basis": "test basis: recorded content digest, not a credential",
    }


def _at21_baseline(entries, a02=False):
    doc = {"metadata": {"scanner_version": "8.30.1", "policy": "explicit reviewed allowlist"},
           "entries": entries}
    if a02 or any(e["classification"] == "DOCUMENTED-HISTORICAL-A-02" for e in entries):
        doc["a02_note"] = "A-02 stays open; issuer-side revocation unverified."
    return tmpfile(json.dumps(doc), ".json")


def _at21_report(rows):
    return tmpfile(json.dumps(rows), ".json")


def test_at21_secret_baseline_allowlist():
    """Round-20 direction 3: the all-refs classification must be an explicit,
    reviewed fingerprint allowlist -- a new finding under docs/ or evidence/
    FAILS (the Round-19 category exemption admitted that class)."""
    fp = "a" * 40 + ":docs/digests.json:generic-api-key:12"
    # 1. exact documented allowlisted fingerprint passes
    r = run("secret_baseline_check.py", "--report",
            _at21_report([{"Fingerprint": fp, "RuleID": "generic-api-key",
                           "File": "docs/digests.json", "StartLine": 12}]),
            "--baseline", _at21_baseline([_at21_baseline_entry(fp)]))
    check("AT-21", "exact documented allowlisted fingerprint passes",
          r.returncode == 0 and "PASS" in r.stdout, r.stdout)
    # 2. a NEW finding under docs/ (unknown fingerprint) FAILS -- the core fix
    new_fp = "b" * 40 + ":docs/new_evidence.md:generic-api-key:3"
    r = run("secret_baseline_check.py", "--report",
            _at21_report([{"Fingerprint": new_fp, "RuleID": "generic-api-key",
                           "File": "docs/new_evidence.md", "StartLine": 3}]),
            "--baseline", _at21_baseline([_at21_baseline_entry(fp)]))
    check("AT-21", "new finding under docs/ fails (no category exemption)",
          r.returncode == 1 and "UNKNOWN FINGERPRINT" in r.stdout
          and "docs/new_evidence.md" in r.stdout, r.stdout)
    # 3. same failure class under evidence/
    ev_fp = "c" * 40 + ":evidence/snap.json:generic-api-key:9"
    r = run("secret_baseline_check.py", "--report",
            _at21_report([{"Fingerprint": ev_fp, "RuleID": "generic-api-key",
                           "File": "evidence/snap.json", "StartLine": 9}]),
            "--baseline", _at21_baseline([_at21_baseline_entry(fp)]))
    check("AT-21", "new finding under evidence/ fails",
          r.returncode == 1 and "UNKNOWN FINGERPRINT" in r.stdout, r.stdout)
    # 4. fingerprint matches but rule/file contradict the reviewed entry -> FAIL
    r = run("secret_baseline_check.py", "--report",
            _at21_report([{"Fingerprint": fp, "RuleID": "aws-access-token",
                           "File": "docs/digests.json", "StartLine": 12}]),
            "--baseline", _at21_baseline([_at21_baseline_entry(fp)]))
    check("AT-21", "matched fingerprint with contradicting rule fails",
          r.returncode == 1 and "contradiction" in r.stdout, r.stdout)
    # 5. baseline integrity: an entry without basis/source/classification fails
    broken = _at21_baseline_entry(fp); del broken["basis"]
    r = run("secret_baseline_check.py", "--report",
            _at21_report([{"Fingerprint": fp, "RuleID": "generic-api-key",
                           "File": "docs/digests.json", "StartLine": 12}]),
            "--baseline", _at21_baseline([broken]))
    check("AT-21", "undocumented baseline entry fails (allowlist is review material)",
          r.returncode == 1 and "missing required field" in r.stdout, r.stdout)
    # 6. H-1 class disclosed with the A-02 note, not silently suppressed
    h1_fp = "d" * 40 + ":system_diagnosis.py:gcp-api-key:8"
    r = run("secret_baseline_check.py", "--report",
            _at21_report([{"Fingerprint": h1_fp, "RuleID": "gcp-api-key",
                           "File": "system_diagnosis.py", "StartLine": 8}]),
            "--baseline", _at21_baseline(
                [_at21_baseline_entry(h1_fp, rule="gcp-api-key",
                                      file="system_diagnosis.py",
                                      cls="DOCUMENTED-HISTORICAL-A-02")]))
    check("AT-21", "H-1 fingerprint disclosed with the open A-02 blocker",
          r.returncode == 0 and "A-02 open" in r.stdout, r.stdout)
    # 7. empty report passes (allowlist is not a must-appear set)
    r = run("secret_baseline_check.py", "--report", _at21_report([]),
            "--baseline", _at21_baseline([_at21_baseline_entry(fp)]))
    check("AT-21", "empty report passes", r.returncode == 0, r.stdout)
    # 8. missing report fails closed (scanner produced nothing)
    r = run("secret_baseline_check.py", "--report", "/nonexistent-report.json",
            "--baseline", _at21_baseline([_at21_baseline_entry(fp)]))
    check("AT-21", "missing report fails closed", r.returncode == 1
          and "does not exist" in r.stdout, r.stdout)
    # 9. duplicate report rows for one fingerprint collapse honestly
    row = {"Fingerprint": fp, "RuleID": "generic-api-key",
           "File": "docs/digests.json", "StartLine": 12}
    r = run("secret_baseline_check.py", "--report", _at21_report([row, dict(row)]),
            "--baseline", _at21_baseline([_at21_baseline_entry(fp)]))
    check("AT-21", "duplicate report rows for one fingerprint collapse",
          r.returncode == 0 and "2 report row(s) / 1 distinct" in r.stdout, r.stdout)
    # 10. the real tracked baseline is well-formed and self-consistent
    real = os.path.join(REPO, "docs", "evidence", "ec16-ec26", "secret_scan_baseline.json")
    doc = json.load(open(real))
    check("AT-21", "real baseline has 47 documented entries with required fields",
          len(doc["entries"]) == 47 and all(
              all(k in e for k in ("fingerprint", "rule", "file", "classification",
                                   "source", "basis")) for e in doc["entries"]),
          str(len(doc.get("entries", []))))
    h1 = [e for e in doc["entries"] if e["classification"] == "DOCUMENTED-HISTORICAL-A-02"]
    check("AT-21", "real baseline keeps the two H-1 fingerprints separate with the A-02 note",
          len(h1) == 2 and "A-02" in doc.get("a02_note", ""), str(len(h1)))


# --------------------- AT-22 (EC-23 controls, Round-20)

def test_at22_secret_scan_controls():
    """Round-20 direction 3: a scanner crash or invocation error must not
    satisfy the positive control; the control must verify the expected finding,
    the fired rule, and key redaction; the clean negative control must pass."""
    key = "AKIA" + "T7W2QZ5XN4PRM3KV"
    os.environ["SEEDED_KEY_AT22"] = key
    good_row = {"RuleID": "aws-access-token", "File": "fixture.txt",
                "StartLine": 1, "Secret": "REDACTED", "Match": "REDACTED"}
    rep = lambda rows: tmpfile(json.dumps(rows), ".json")
    # 1. genuine detection: exit 101 + expected finding + redacted key -> PASS
    r = run("secret_control_check.py", "positive", "--report", rep([good_row]),
            "--exit-code", "101", "--expect-rule", "aws-access-token",
            "--expect-file", "fixture.txt", "--key-env", "SEEDED_KEY_AT22")
    check("AT-22", "genuine detection with redaction passes",
          r.returncode == 0 and "redaction verified" in r.stdout, r.stdout)
    # 2. scanner exit 0 (nothing found) fails the control
    r = run("secret_control_check.py", "positive", "--report", rep([good_row]),
            "--exit-code", "0", "--expect-rule", "aws-access-token",
            "--expect-file", "fixture.txt", "--key-env", "SEEDED_KEY_AT22")
    check("AT-22", "exit 0 (key NOT detected) fails", r.returncode == 1
          and "NOT detected" in r.stdout, r.stdout)
    # 3. scanner execution error (exit 1) must not satisfy the control
    r = run("secret_control_check.py", "positive", "--report", rep([good_row]),
            "--exit-code", "1", "--expect-rule", "aws-access-token",
            "--expect-file", "fixture.txt", "--key-env", "SEEDED_KEY_AT22")
    check("AT-22", "scanner execution error does not satisfy the positive control",
          r.returncode == 1 and "execution error" in r.stdout.lower()
          and "masquerading" in r.stdout, r.stdout)
    # 4. missing report fails closed
    r = run("secret_control_check.py", "positive", "--report", "/nonexistent-at22.json",
            "--exit-code", "101", "--expect-rule", "aws-access-token",
            "--expect-file", "fixture.txt", "--key-env", "SEEDED_KEY_AT22")
    check("AT-22", "missing report fails closed", r.returncode == 1
          and "no report" in r.stdout, r.stdout)
    # 5. malformed report fails closed
    bad = tmpfile("{not json", ".json")
    r = run("secret_control_check.py", "positive", "--report", bad,
            "--exit-code", "101", "--expect-rule", "aws-access-token",
            "--expect-file", "fixture.txt", "--key-env", "SEEDED_KEY_AT22")
    check("AT-22", "malformed report fails closed", r.returncode == 1
          and "not parseable" in r.stdout, r.stdout)
    # 6. unredacted key in the report fails the control
    leaked = dict(good_row, Secret=key)
    r = run("secret_control_check.py", "positive", "--report", rep([leaked]),
            "--exit-code", "101", "--expect-rule", "aws-access-token",
            "--expect-file", "fixture.txt", "--key-env", "SEEDED_KEY_AT22")
    check("AT-22", "unredacted key in the report fails",
          r.returncode == 1 and "redaction FAILED" in r.stdout, r.stdout)
    # 7. wrong rule fired fails
    wrong = dict(good_row, RuleID="gcp-api-key")
    r = run("secret_control_check.py", "positive", "--report", rep([wrong]),
            "--exit-code", "101", "--expect-rule", "aws-access-token",
            "--expect-file", "fixture.txt", "--key-env", "SEEDED_KEY_AT22")
    check("AT-22", "expected rule not firing fails", r.returncode == 1
          and "expected rule did not fire" in r.stdout, r.stdout)
    # 8. unexpected additional finding fails
    extra = dict(good_row, RuleID="gcp-api-key", File="other.txt")
    r = run("secret_control_check.py", "positive", "--report", rep([good_row, extra]),
            "--exit-code", "101", "--expect-rule", "aws-access-token",
            "--expect-file", "fixture.txt", "--key-env", "SEEDED_KEY_AT22")
    check("AT-22", "unexpected additional finding fails", r.returncode == 1
          and "unexpected additional" in r.stdout, r.stdout)
    # 9. the raw key must not appear anywhere in the report bytes
    withkey = tmpfile(json.dumps([dict(good_row)]) + "\n# " + key + "\n", ".json")
    r = run("secret_control_check.py", "positive", "--report", withkey,
            "--exit-code", "101", "--expect-rule", "aws-access-token",
            "--expect-file", "fixture.txt", "--key-env", "SEEDED_KEY_AT22")
    check("AT-22", "raw key appearing anywhere in the report fails",
          r.returncode == 1 and "raw seeded key appears" in r.stdout, r.stdout)
    # 10. clean repository: exit 0 + zero findings passes; anything else fails
    r = run("secret_control_check.py", "negative", "--report", rep([]), "--exit-code", "0")
    check("AT-22", "clean repository passes the negative control",
          r.returncode == 0, r.stdout)
    r = run("secret_control_check.py", "negative", "--report", rep([good_row]),
            "--exit-code", "101")
    check("AT-22", "clean repository with a finding fails",
          r.returncode == 1 and "false positive" in r.stdout, r.stdout)
    r = run("secret_control_check.py", "negative", "--report", rep([]), "--exit-code", "1")
    check("AT-22", "clean repository with scanner error fails",
          r.returncode == 1 and "expected 0" in r.stdout, r.stdout)


# --------------------- AT-23 (EC-26 trigger liveness, Round-20)

def test_at23_workflow_trigger_liveness():
    """Round-20 direction 5: a workflow under .github/workflows/ with the named
    job but only a manual-dispatch trigger must NOT satisfy an ACTIVE rule;
    liveness requires the declared events to actually trigger, with branch
    coverage of the enforcement branch."""
    d = tempfile.mkdtemp()
    os.makedirs(os.path.join(d, ".github", "workflows"), exist_ok=True)
    wf = os.path.join(d, ".github", "workflows", "gate.yml")

    def write_wf(text):
        open(wf, "w").write(text)

    reg = {"tests": {}}
    open(os.path.join(d, "registry.json"), "w").write(json.dumps(reg))
    os.makedirs(os.path.join(d, "docs", "directives"), exist_ok=True)

    def idx(enf):
        text = ("# index\n<!-- BEGIN RULE INVENTORY -->\n"
                "reading-note-path: docs/note.md\n"
                "pinned-blobs: b6bb35e84956d84a8c9d9e65734c212baf3ed075\n"
                "| Rule ID | Status | Enforcement point(s) | Acceptance test(s) |\n"
                "|---|---|---|---|\n"
                "| EC-X | ACTIVE | " + enf + " | NONE |\n"
                "<!-- END RULE INVENTORY -->\n")
        open(os.path.join(d, "docs", "directives", "INDEX.md"), "w").write(text)
        open(os.path.join(d, "docs", "note.md"), "w").write(
            "blob b6bb35e84956d84a8c9d9e65734c212baf3ed075 read\n| EC-X | -> mapped |\n")

    def run_idx():
        return run("reading_note_map_check.py", "--index", "docs/directives/INDEX.md",
                   "--repo-root", d, "--registry", "registry.json")

    job = "jobs:\n  secret-scan:\n    name: Secret Scan\n    runs-on: ubuntu-latest\n"
    # 1. valid push/PR workflow with matching annotation passes
    write_wf("name: gate\non:\n  push:\n    branches: ['**']\n  pull_request:\n"
             "    branches: ['**']\n" + job)
    idx("workflow:.github/workflows/gate.yml#secret-scan~events=push+pull_request~governs=**")
    r = run_idx()
    check("AT-23", "valid push/PR workflow with declared events passes",
          r.returncode == 0, r.stdout)
    # 2. path-correct, job exists, but MANUAL-ONLY -> fails
    write_wf("name: gate\non:\n  workflow_dispatch:\n" + job)
    idx("workflow:.github/workflows/gate.yml#secret-scan~events=push+pull_request~governs=**")
    r = run_idx()
    check("AT-23", "manual-only (workflow_dispatch) workflow fails",
          r.returncode == 1 and "NOT triggered by" in r.stdout, r.stdout)
    # 3. push trigger with a branch filter excluding main -> fails
    write_wf("name: gate\non:\n  push:\n    branches: [feature-x]\n" + job)
    idx("workflow:.github/workflows/gate.yml#secret-scan~events=push~governs=**")
    r = run_idx()
    check("AT-23", "branch filter excluding the enforcement branch fails",
          r.returncode == 1 and "exclude the enforcement branch" in r.stdout, r.stdout)
    # 4. branches-ignore covering main -> fails
    write_wf("name: gate\non:\n  push:\n    branches-ignore: [main]\n" + job)
    idx("workflow:.github/workflows/gate.yml#secret-scan~events=push~governs=**")
    r = run_idx()
    check("AT-23", "branches-ignore covering the enforcement branch fails",
          r.returncode == 1 and "ignores the enforcement branch" in r.stdout, r.stdout)
    # 5. ACTIVE reference without the events annotation -> fails
    write_wf("name: gate\non:\n  push:\n" + job)
    idx("workflow:.github/workflows/gate.yml#secret-scan")
    r = run_idx()
    check("AT-23", "ACTIVE workflow reference without ~events declaration fails",
          r.returncode == 1 and "must declare the events" in r.stdout, r.stdout)
    # 6. unknown event token in the annotation -> fails
    idx("workflow:.github/workflows/gate.yml#secret-scan~events=banana")
    r = run_idx()
    check("AT-23", "unknown event token in the annotation fails",
          r.returncode == 1 and "unknown event token" in r.stdout, r.stdout)
    # 7. list-form on: [push, pull_request] passes
    write_wf("name: gate\non: [push, pull_request]\n" + job)
    idx("workflow:.github/workflows/gate.yml#secret-scan~events=push+pull_request~governs=**")
    r = run_idx()
    check("AT-23", "list-form on: [push, pull_request] passes",
          r.returncode == 0, r.stdout)
    # 8. no on: block at all -> fails
    write_wf("name: gate\n" + job)
    idx("workflow:.github/workflows/gate.yml#secret-scan~events=push~governs=**")
    r = run_idx()
    check("AT-23", "workflow with no on: block fails",
          r.returncode == 1 and "no 'on:' triggers" in r.stdout, r.stdout)
    # 9. pull_request declared but workflow only has push -> fails
    write_wf("name: gate\non:\n  push:\n" + job)
    idx("workflow:.github/workflows/gate.yml#secret-scan~events=pull_request~governs=**")
    r = run_idx()
    check("AT-23", "declared pull_request missing from the workflow fails",
          r.returncode == 1 and "NOT triggered by 'pull_request'" in r.stdout, r.stdout)
    # 10. PROPOSED rule citing a manual-only draft WITHOUT an events annotation
    #     is allowed (activation path: the draft is inert; events are declared
    #     at activation time, when the rule goes ACTIVE)
    write_wf("name: gate\non:\n  workflow_dispatch:\n" + job)
    text = open(os.path.join(d, "docs", "directives", "INDEX.md")).read()
    open(os.path.join(d, "docs", "directives", "INDEX.md"), "w").write(
        text.replace("~events=pull_request", "").replace("| EC-X | ACTIVE |", "| EC-X | PROPOSED |"))
    r = run("reading_note_map_check.py", "--index", "docs/directives/INDEX.md",
            "--repo-root", d, "--registry", "registry.json", "--include-proposed")
    check("AT-23", "PROPOSED rule citing a manual-only draft is allowed (activation path)",
          r.returncode == 0, r.stdout)
    # 10b. a PROPOSED row that DOES declare an event the draft lacks still fails
    #      (whenever an ~events declaration is present it is validated -- found by
    #      the Round-21-followup soundness audit, now an explicit check)
    open(os.path.join(d, "docs", "directives", "INDEX.md"), "w").write(
        text.replace("| EC-X | ACTIVE |", "| EC-X | PROPOSED |"))
    r = run("reading_note_map_check.py", "--index", "docs/directives/INDEX.md",
            "--repo-root", d, "--registry", "registry.json", "--include-proposed")
    check("AT-23", "PROPOSED rule with a declared event the draft lacks fails",
          r.returncode == 1 and "NOT triggered by" in r.stdout, r.stdout)


# --------------------- AT-24 (EC-24 remote-ref verification, Round-20)

def test_at24_pushlog_remote_verification():
    """Round-20 direction 4: every confirmation must be verified against the
    live remote ref; stale confirmations are HISTORICAL, not current; the
    branch under test's tip must be declared; fabricated evidence, mismatched
    branches, unpaired confirmations, and missing remote verification fail."""
    import subprocess as sp
    sha = lambda c: c * 40
    d = tempfile.mkdtemp()
    bare = os.path.join(d, "remote.git")
    sp.run(["git", "init", "-q", "--bare", bare], check=True)
    src = os.path.join(d, "src")
    sp.run(["git", "init", "-q", src], check=True)

    def commit_and_push(tag):
        sp.run(["git", "-C", src, "-c", "user.email=t@t", "-c", "user.name=t",
                "commit", "-q", "--allow-empty", "-m", tag], check=True)
        sp.run(["git", "-C", src, "push", "-q", bare, "HEAD:refs/heads/topic/b"], check=True)
        return sp.run(["git", "-C", src, "rev-parse", "HEAD"],
                      capture_output=True, text=True).stdout.strip()

    s1 = commit_and_push("first")
    s2 = commit_and_push("second")

    def log_with(entries):
        return tmpfile(entries, ".md")

    def rv(log, branch="topic/b", tip=s2, remote=bare, extra=()):
        return run("evidence_sequence_check.py", log, "--pushlog-schema",
                   "--remote-verify", "--remote", remote,
                   "--branch", branch, "--expect-tip", tip, *extra)

    # 1. CURRENT confirmation: live tip matches the recorded remote_sha -> PASS
    cur = log_with("## 2026-10-20T10:00:00Z | P-1-INTENT | push-intent | branch=topic/b credential=OWNER-PAT-1\n"
                   "## 2026-10-20T11:00:00Z | P-1-CONFIRM | push-confirm | branch=topic/b remote_sha=%s "
                   "credential=OWNER-PAT-1 evidence=\"git ls-remote origin refs/heads/topic/b -> %s\"\n" % (s2, s2))
    r = rv(cur)
    check("AT-24", "CURRENT confirmation (live tip == remote_sha) passes",
          r.returncode == 0 and "CURRENT" in r.stdout, r.stdout)
    # 2. valid-but-stale SHA: CONFIRM for the older tip; current tip undeclared -> FAIL
    stale = log_with("## 2026-10-20T10:00:00Z | P-2-INTENT | push-intent | branch=topic/b credential=OWNER-PAT-1\n"
                     "## 2026-10-20T11:00:00Z | P-2-CONFIRM | push-confirm | branch=topic/b remote_sha=%s "
                     "credential=OWNER-PAT-1 evidence=\"git ls-remote origin refs/heads/topic/b -> %s\"\n" % (s1, s1))
    r = rv(stale)
    check("AT-24", "stale confirmation is HISTORICAL and the undeclared current tip FAILS",
          r.returncode == 1 and "HISTORICAL" in r.stdout and "UNDECLARED" in r.stdout,
          r.stdout)
    # 3. (Round-21 direction 4) the current tip carries an unpaired terminal
    #    INTENT introduced by the exact candidate, with a verified pre-push tip
    #    -> PENDING pass. The log must live in a git repository so the checker
    #    can bind the intent to the commit that introduced it.
    logdir = os.path.join(src, "docs", "evidence", "push-log")
    os.makedirs(logdir, exist_ok=True)
    logpath = os.path.join(logdir, "PUSH_LOG.md")

    def commit_log(text, tag):
        open(logpath, "w").write(text)
        sp.run(["git", "-C", src, "add", "docs"], check=True)
        sp.run(["git", "-C", src, "-c", "user.email=t@t", "-c", "user.name=t",
                "commit", "-q", "--allow-empty", "-m", tag], check=True)
        return sp.run(["git", "-C", src, "rev-parse", "HEAD"],
                      capture_output=True, text=True).stdout.strip()

    def push_b():
        sp.run(["git", "-C", src, "push", "-q", bare, "HEAD:refs/heads/topic/b"], check=True)

    cA = commit_log("## 2026-10-20T10:00:00Z | P-3-INTENT | push-intent | "
                    "branch=topic/b expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n", "log A")
    push_b()
    cB = commit_log("## 2026-10-20T10:00:00Z | P-3-INTENT | push-intent | "
                    "branch=topic/b expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
                    "## 2026-10-20T11:00:00Z | P-3-CONFIRM | push-confirm | branch=topic/b remote_sha=%s "
                    "credential=OWNER-PAT-1 evidence=\"git ls-remote origin refs/heads/topic/b -> %s\"\n"
                    "## 2026-10-20T12:00:00Z | P-4-INTENT | push-intent | branch=topic/b "
                    "expected_remote_tip_before_push=%s credential=OWNER-PAT-1\n" % (cA, cA, cA),
                    "log B")
    push_b()
    r = run("evidence_sequence_check.py", logpath, "--pushlog-schema", "--remote-verify",
            "--remote", bare, "--branch", "topic/b", "--expect-tip", cB)
    check("AT-24", "tip declared by a terminal INTENT bound to the exact candidate is PENDING, not a failure",
          r.returncode == 0 and "PENDING" in r.stdout and "P-4-INTENT" in r.stdout
          and "BOUND" in r.stdout, r.stdout)
    # 4. RECONCILE for the branch does not consume the terminal INTENT
    cC = commit_log("## 2026-10-20T10:00:00Z | P-3-INTENT | push-intent | "
                    "branch=topic/b expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
                    "## 2026-10-20T11:00:00Z | P-3-CONFIRM | push-confirm | branch=topic/b remote_sha=%s "
                    "credential=OWNER-PAT-1 evidence=\"git ls-remote origin refs/heads/topic/b -> %s\"\n"
                    "## 2026-10-20T11:30:00Z | P-6 | push-reconcile | branch=topic/b remote_sha=%s "
                    "credential=OWNER-PAT-1 evidence=\"git ls-remote origin refs/heads/topic/b -> %s verified now\"\n"
                    "## 2026-10-20T12:00:00Z | P-7-INTENT | push-intent | branch=topic/b "
                    "expected_remote_tip_before_push=%s credential=OWNER-PAT-1\n"
                    % (cA, cA, cB, cB, cB), "log C")
    push_b()
    r = run("evidence_sequence_check.py", logpath, "--pushlog-schema", "--remote-verify",
            "--remote", bare, "--branch", "topic/b", "--expect-tip", cC)
    check("AT-24", "a RECONCILE never consumes the terminal INTENT (PENDING preserved)",
          r.returncode == 0 and "PENDING" in r.stdout and "P-7-INTENT" in r.stdout
          and "P-6" in r.stdout, r.stdout)
    # 5. mismatched branch: CONFIRM claims branch=topic/b but the evidence names only branch=c
    mism = log_with("## 2026-10-20T10:00:00Z | P-8-INTENT | push-intent | branch=topic/b credential=OWNER-PAT-1\n"
                    "## 2026-10-20T11:00:00Z | P-8-CONFIRM | push-confirm | branch=topic/b remote_sha=%s "
                    "credential=OWNER-PAT-1 evidence=\"git ls-remote origin refs/heads/c -> %s\"\n" % (s2, s2))
    r = rv(mism)
    check("AT-24", "confirmation whose evidence names a different branch fails",
          r.returncode == 1 and "mismatched branch" in r.stdout, r.stdout)
    # 6. fabricated evidence: neither field nor body contains the recorded SHA
    fab = log_with("## 2026-10-20T10:00:00Z | P-9-INTENT | push-intent | branch=topic/b credential=OWNER-PAT-1\n"
                   "## 2026-10-20T11:00:00Z | P-9-CONFIRM | push-confirm | branch=topic/b remote_sha=%s "
                   "credential=OWNER-PAT-1 evidence=\"git ls-remote origin refs/heads/topic/b (ok)\"\n" % sha("e"))
    r = rv(fab)
    check("AT-24", "fabricated evidence string fails",
          r.returncode == 1 and "fabricated or unsupported evidence" in r.stdout, r.stdout)
    # 7. unpaired confirmation fails under remote verification too
    unpaired = log_with("## 2026-10-20T11:00:00Z | P-10-CONFIRM | push-confirm | branch=topic/b "
                        "remote_sha=%s credential=OWNER-PAT-1 evidence=\"git ls-remote origin "
                        "refs/heads/topic/b -> %s\"\n" % (s2, s2))
    r = rv(unpaired)
    check("AT-24", "unpaired confirmation fails",
          r.returncode == 1 and "no preceding" in r.stdout, r.stdout)
    # 8. missing remote verification: unreachable remote fails closed
    r = rv(cur, remote=os.path.join(d, "nonexistent.git"))
    check("AT-24", "unreachable remote fails closed (missing remote verification)",
          r.returncode == 1 and "could not be verified" in r.stdout
          and "fail closed" in r.stdout, r.stdout)
    # 9. runs-manifest contradiction: cited run's head_sha differs -> FAIL.
    #    The confirmation records the CURRENT live tip (cC) so the classification
    #    is CURRENT and the assertion isolates the citation semantics
    man = tmpfile(json.dumps({"runs": {"999999": {"head_sha": sha("9"),
                                                  "conclusion": "success"}}}), ".json")
    contra = log_with("## 2026-10-20T10:00:00Z | P-11-INTENT | push-intent | branch=topic/b credential=OWNER-PAT-1\n"
                      "## 2026-10-20T11:00:00Z | P-11-CONFIRM | push-confirm | branch=topic/b remote_sha=%s "
                      "credential=OWNER-PAT-1 evidence=\"git ls-remote origin refs/heads/topic/b -> %s; gate run 999999\"\n"
                      % (cC, cC))
    r = rv(contra, tip=cC, extra=("--runs-manifest", man))
    check("AT-24", "cited run whose head_sha contradicts the recorded SHA fails",
          r.returncode == 1 and "evidence does not support" in r.stdout, r.stdout)
    # 10. runs-manifest corroboration passes
    man2 = tmpfile(json.dumps({"runs": {"999999": {"head_sha": cC,
                                                   "conclusion": "success"}}}), ".json")
    r = rv(contra, tip=cC, extra=("--runs-manifest", man2))
    check("AT-24", "cited run corroborating the recorded SHA passes",
          r.returncode == 0, r.stdout)
    # 11. expect-tip mismatch (checked-out candidate is not the live tip) fails
    r = rv(cur, tip=s1)
    check("AT-24", "candidate that is not the live tip fails",
          r.returncode == 1 and "not the expected candidate" in r.stdout, r.stdout)
    # 12. the real mainline log classifies every confirmation without fabricating.
    #     The off-lineage historical SHAs must exist in the local object store for
    #     the evidence-back proof (a CI checkout of all branch tips does not fetch
    #     unreachable commits; there the full verification -- including the commits
    #     API -- is performed by the push-log-integrity job, which always runs with
    #     --repo and a token in the same workflow).
    real = os.path.join(REPO, "docs", "evidence", "push-log", "PUSH_LOG.md")
    r = run("evidence_sequence_check.py", real, "--pushlog-schema")
    check("AT-24", "the real mainline log passes the schema (identifier pairing)",
          r.returncode == 0, r.stdout)
    have_hist = sp.run(["git", "-C", REPO, "cat-file", "-e",
                        "17eb88184c9f461c6bc9ee592f916f7ff4c921c4^{commit}"],
                       capture_output=True).returncode == 0
    if have_hist:
        r = run("evidence_sequence_check.py", real, "--pushlog-schema", "--remote-verify")
        check("AT-24", "the real mainline log passes remote classification (no fabrication)",
              r.returncode == 0, r.stdout)
    else:
        print("AT-24 NOTE: off-lineage historical objects absent from this checkout; "
              "the full remote classification runs in the push-log-integrity CI job "
              "(--repo + token + manifest)")


# ------------------------------------------- AT-25 (EC-23 disclosure, Round-21)

def test_at25_secret_control_disclosure():
    """Round-21 direction 3: the EC-23 control checker must NEVER print
    Secret/Match field values or raw key material on ANY path (success or
    failure); a requested capture log must exist and be readable (fail
    closed); missing or malformed reports fail; the negative mode requires a
    parsed empty-finding report rather than a trusted exit status."""
    import subprocess as sp
    # synthetic test key, never a real credential; assembled from two
    # non-adjacent parts so the tracked source bytes cannot match the
    # generic-api-key/aws-access-token rules (same discipline as the seeded
    # workflow control)
    KEY = "AKIA" + "ZZR21SYNTHETICKEYQQ7X"
    d = tempfile.mkdtemp()

    def write_json(name, obj):
        p = os.path.join(d, name)
        json.dump(obj, open(p, "w"))
        return p

    def write_raw(name, text):
        p = os.path.join(d, name)
        open(p, "w").write(text)
        return p

    def ctrl(mode, *args):
        env = dict(os.environ)
        env["SEEDED_KEY"] = KEY
        return subprocess.run(
            [PY, os.path.join(CI, "secret_control_check.py"), mode, *args],
            capture_output=True, text=True, env=env, shell=False)

    finding_ok = [{"RuleID": "aws-access-token", "File": "fixture.txt",
                   "StartLine": 1, "Secret": "REDACTED", "Match": "REDACTED"}]
    log_ok = write_raw("scan.log", "scanner output with no key material\n")

    # 1. success path: PASS and the key never appears in stdout/stderr
    r = ctrl("positive", "--report", write_json("ok.json", finding_ok),
             "--exit-code", "101", "--log-file", log_ok)
    check("AT-25", "redacted positive control passes",
          r.returncode == 0 and "PASS" in r.stdout, r.stdout + r.stderr)
    check("AT-25", "the synthetic key never appears in the checker output (success path)",
          KEY not in r.stdout and KEY not in r.stderr, r.stdout + r.stderr)
    # 2. report Secret field carries the raw key -> FAIL, key never printed
    leak_secret = [{"RuleID": "aws-access-token", "File": "fixture.txt",
                    "StartLine": 1, "Secret": KEY, "Match": "REDACTED"}]
    r = ctrl("positive", "--report", write_json("leak_secret.json", leak_secret),
             "--exit-code", "101")
    check("AT-25", "Secret field carrying the raw key fails redaction",
          r.returncode == 1 and "redaction FAILED" in r.stdout, r.stdout)
    check("AT-25", "the leaked Secret value is never printed",
          KEY not in r.stdout and KEY not in r.stderr, r.stdout + r.stderr)
    # 3. Match field carries the raw key -> FAIL, key never printed
    leak_match = [{"RuleID": "aws-access-token", "File": "fixture.txt",
                   "StartLine": 1, "Secret": "REDACTED", "Match": KEY}]
    r = ctrl("positive", "--report", write_json("leak_match.json", leak_match),
             "--exit-code", "101")
    check("AT-25", "Match field carrying the raw key fails redaction",
          r.returncode == 1 and "redaction FAILED" in r.stdout, r.stdout)
    check("AT-25", "the leaked Match value is never printed",
          KEY not in r.stdout and KEY not in r.stderr, r.stdout + r.stderr)
    # 4. raw key present in the report BYTES (extra field) -> FAIL, not printed
    r = ctrl("positive", "--report",
             write_raw("leak_bytes.json", json.dumps(finding_ok)[:-1] + ', "Extra": "%s"}' % KEY),
             "--exit-code", "101")
    check("AT-25", "raw key in the report bytes fails redaction without being printed",
          r.returncode == 1 and KEY not in r.stdout and KEY not in r.stderr,
          r.stdout + r.stderr)
    # 5. raw key in the captured scanner log -> FAIL, not printed
    leaky_log = write_raw("leaky.log", "scanning...\nfound key %s in fixture\n" % KEY)
    r = ctrl("positive", "--report", write_json("ok2.json", finding_ok),
             "--exit-code", "101", "--log-file", leaky_log)
    check("AT-25", "raw key in the captured log fails redaction without being printed",
          r.returncode == 1 and KEY not in r.stdout and KEY not in r.stderr,
          r.stdout + r.stderr)
    # 6. requested capture log MISSING -> fail closed
    r = ctrl("positive", "--report", write_json("ok3.json", finding_ok),
             "--exit-code", "101", "--log-file", os.path.join(d, "missing.log"))
    check("AT-25", "a requested but missing capture log fails closed",
          r.returncode == 1 and "MISSING" in r.stdout, r.stdout)
    # 7. malformed report (not JSON) -> FAIL
    r = ctrl("positive", "--report", write_raw("bad.json", "not json at all"),
             "--exit-code", "101")
    check("AT-25", "malformed positive report fails",
          r.returncode == 1 and "parseable JSON" in r.stdout, r.stdout)
    # 8. missing report -> FAIL
    r = ctrl("positive", "--report", os.path.join(d, "absent.json"), "--exit-code", "101")
    check("AT-25", "missing positive report fails",
          r.returncode == 1 and "absent" in r.stdout, r.stdout)
    # 9. negative mode: findings present -> FAIL; missing report -> FAIL;
    #    malformed -> FAIL; requested log missing -> FAIL
    r = ctrl("negative", "--report", write_json("neg_find.json", finding_ok),
             "--exit-code", "0")
    check("AT-25", "negative control with findings fails",
          r.returncode == 1 and "finding(s)" in r.stdout, r.stdout)
    r = ctrl("negative", "--report", os.path.join(d, "absent2.json"), "--exit-code", "0")
    check("AT-25", "negative control with a missing report fails (exit status alone is not trusted)",
          r.returncode == 1 and "absent" in r.stdout, r.stdout)
    r = ctrl("negative", "--report", write_raw("neg_bad.json", "{oops"),
             "--exit-code", "0")
    check("AT-25", "negative control with a malformed report fails",
          r.returncode == 1, r.stdout)
    r = ctrl("negative", "--report", write_json("neg_ok.json", []),
             "--exit-code", "0", "--log-file", os.path.join(d, "missing2.log"))
    check("AT-25", "negative control with a requested-but-missing log fails closed",
          r.returncode == 1 and "MISSING" in r.stdout, r.stdout)
    # 10. negative success: parsed empty report + existing log -> PASS, no key
    r = ctrl("negative", "--report", write_json("neg_clean.json", []),
             "--exit-code", "0", "--log-file", log_ok)
    check("AT-25", "verified-empty negative control passes without printing the key",
          r.returncode == 0 and KEY not in r.stdout and KEY not in r.stderr,
          r.stdout + r.stderr)


# --------------------------------- AT-26 (EC-24 intent binding, Round-21)

def test_at26_intent_binding_and_pairing():
    """Round-21 direction 4: confirmations pair their OWN identifier's intent;
    the single accepted terminal intent must be introduced by the exact
    candidate and declare the verified pre-push tip; excess unmatched intents
    fail unless voided; a stale intent can never declare later commits."""
    import subprocess as sp
    d = tempfile.mkdtemp()
    bare = os.path.join(d, "remote.git")
    sp.run(["git", "init", "-q", "--bare", bare], check=True)
    src = os.path.join(d, "src")
    sp.run(["git", "init", "-q", "-b", "main", src], check=True)
    logdir = os.path.join(src, "docs", "evidence", "push-log")
    os.makedirs(logdir, exist_ok=True)
    logpath = os.path.join(logdir, "PUSH_LOG.md")

    def commit_log(text, tag):
        open(logpath, "w").write(text)
        sp.run(["git", "-C", src, "add", "docs"], check=True)
        sp.run(["git", "-C", src, "-c", "user.email=t@t", "-c", "user.name=t",
                "commit", "-q", "--allow-empty", "-m", tag], check=True)
        return sp.run(["git", "-C", src, "rev-parse", "HEAD"],
                      capture_output=True, text=True).stdout.strip()

    def push(refspec="HEAD:refs/heads/b"):
        sp.run(["git", "-C", src, "push", "-q", bare, refspec], check=True)

    def chk(tip):
        return run("evidence_sequence_check.py", logpath, "--pushlog-schema",
                   "--remote-verify", "--remote", bare, "--branch", "b",
                   "--expect-tip", tip)

    # 1. identifier pairing: a CONFIRM with a DIFFERENT id must not consume the intent
    t = ("## 2026-10-20T10:00:00Z | P-1-INTENT | push-intent | branch=b credential=OWNER-PAT-1\n"
         "## 2026-10-20T11:00:00Z | P-9-CONFIRM | push-confirm | branch=b remote_sha="
         + "a" * 40 + " credential=OWNER-PAT-1 evidence=\"ls-remote b -> " + "a" * 40 + "\"\n")
    f = tmpfile(t, ".md")
    r = run("evidence_sequence_check.py", f, "--pushlog-schema")
    check("AT-26", "a confirmation must match its own intent identifier (no cross-consumption)",
          r.returncode == 1 and "identifier P-9" in r.stdout, r.stdout)
    # 2. correct identifier pairing passes schema
    t2 = t.replace("P-9-CONFIRM", "P-1-CONFIRM")
    f2 = tmpfile(t2, ".md")
    r = run("evidence_sequence_check.py", f2, "--pushlog-schema")
    check("AT-26", "matching identifier pairing passes",
          r.returncode == 0, r.stdout)
    # 3. legitimate terminal flow: intent -> confirm -> terminal intent on the exact candidate
    a = commit_log("## 2026-10-20T10:00:00Z | P-2-INTENT | push-intent | branch=b "
                   "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n", "A")
    push()
    b = commit_log("## 2026-10-20T10:00:00Z | P-2-INTENT | push-intent | branch=b "
                   "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
                   "## 2026-10-20T11:00:00Z | P-2-CONFIRM | push-confirm | branch=b remote_sha=%s "
                   "credential=OWNER-PAT-1 evidence=\"git ls-remote r refs/heads/b -> %s\"\n"
                   "## 2026-10-20T12:00:00Z | P-3-INTENT | push-intent | branch=b "
                   "expected_remote_tip_before_push=%s credential=OWNER-PAT-1\n" % (a, a, a), "B")
    push()
    r = chk(b)
    check("AT-26", "terminal intent bound to the exact candidate is accepted as PENDING",
          r.returncode == 0 and "PENDING" in r.stdout and "P-3-INTENT" in r.stdout
          and "BOUND" in r.stdout, r.stdout)
    # 4. a LATER commit without a new intent: the old intent must not declare it
    c = commit_log(open(logpath).read(), "C")  # same log, new commit on top
    push()
    r = chk(c)
    check("AT-26", "a stale intent cannot declare later commits (candidate advanced)",
          r.returncode == 1 and "does not bind to the candidate" in r.stdout
          and "undeclared" in r.stdout, r.stdout)
    # 5. expected_remote_tip_before_push must match the newest recorded branch state
    bad_exp = ("## 2026-10-20T10:00:00Z | P-4-INTENT | push-intent | branch=b "
               "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
               "## 2026-10-20T11:00:00Z | P-4-CONFIRM | push-confirm | branch=b remote_sha=%s "
               "credential=OWNER-PAT-1 evidence=\"git ls-remote r refs/heads/b -> %s\"\n"
               "## 2026-10-20T12:00:00Z | P-5-INTENT | push-intent | branch=b "
               "expected_remote_tip_before_push=%s credential=OWNER-PAT-1\n"
               % (a, a, "f" * 40))
    e = commit_log(bad_exp, "E")
    push()
    r = chk(e)
    check("AT-26", "a declared pre-push tip that mismatches the recorded state fails",
          r.returncode == 1 and "intervening ref update or a stale declaration" in r.stdout,
          r.stdout)
    # 6. multiple unmatched intents fail...
    two = ("## 2026-10-20T10:00:00Z | P-6-INTENT | push-intent | branch=b "
           "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
           "## 2026-10-20T10:30:00Z | P-7-INTENT | push-intent | branch=b "
           "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n")
    g = commit_log(two, "G")
    push()
    r = chk(g)
    check("AT-26", "multiple unmatched intents for the branch under test fail",
          r.returncode == 1 and "multiple unmatched intents" in r.stdout, r.stdout)
    # 7. ...unless the excess intents are voided by intent-reconcile records.
    #    Fresh identifiers: the voided intents must be INTRODUCED by the same
    #    candidate commit that carries the reconcile record (an intent introduced
    #    by an earlier commit can never declare a later push)
    voided = ("## 2026-10-20T10:00:00Z | P-9-INTENT | push-intent | branch=b "
              "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
              "## 2026-10-20T10:30:00Z | P-10-INTENT | push-intent | branch=b "
              "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
              "## 2026-10-20T10:40:00Z | P-11 | intent-reconcile | voids=P-9-INTENT "
              "branch=b credential=OWNER-PAT-1 evidence=\"superseded by P-10 on the same push\"\n")
    h = commit_log(voided, "H")
    push()
    r = chk(h)
    check("AT-26", "excess unmatched intents voided by intent-reconcile are tolerated",
          r.returncode == 0 and "PENDING" in r.stdout and "P-10-INTENT" in r.stdout,
          r.stdout)
    # 8. ABSENT declared while a prior record exists fails
    absent_bad = ("## 2026-10-20T10:00:00Z | P-10-INTENT | push-intent | branch=b "
                  "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
                  "## 2026-10-20T11:00:00Z | P-10-CONFIRM | push-confirm | branch=b remote_sha=%s "
                  "credential=OWNER-PAT-1 evidence=\"git ls-remote r refs/heads/b -> %s\"\n"
                  "## 2026-10-20T12:00:00Z | P-11-INTENT | push-intent | branch=b "
                  "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n" % (a, a))
    i2 = commit_log(absent_bad, "I")
    push()
    r = chk(i2)
    check("AT-26", "ABSENT pre-push tip contradicted by a prior record fails",
          r.returncode == 1 and "already records branch state" in r.stdout, r.stdout)


# ------------------------------ AT-27 (EC-24 historical proof, Round-21)

def test_at27_historical_proof():
    """Round-21 direction 5: a HISTORICAL confirmation must be proven real --
    ancestry from the live tip, or a reviewed evidence-back record whose SHA
    exists; every cited run must verify (manifest, else fail closed); failed
    runs are disclosed and cannot prove a push; dangling backing fails."""
    import subprocess as sp
    d = tempfile.mkdtemp()
    bare = os.path.join(d, "remote.git")
    sp.run(["git", "init", "-q", "--bare", bare], check=True)
    src = os.path.join(d, "src")
    sp.run(["git", "init", "-q", "-b", "main", src], check=True)

    def commit(repo, tag, path=None, text=None):
        if path and text is not None:
            os.makedirs(os.path.dirname(os.path.join(repo, path)), exist_ok=True)
            open(os.path.join(repo, path), "w").write(text)
            sp.run(["git", "-C", repo, "add", "docs"], check=True)
        sp.run(["git", "-C", repo, "-c", "user.email=t@t", "-c", "user.name=t",
                "commit", "-q", "--allow-empty", "-m", tag], check=True)
        return sp.run(["git", "-C", repo, "rev-parse", "HEAD"],
                      capture_output=True, text=True).stdout.strip()

    logrel = "docs/evidence/push-log/PUSH_LOG.md"
    logpath = os.path.join(src, logrel)
    # build an off-lineage history: root -> A pushed to b, then b rewritten to B
    # (root exists so the reset below always has a parent to fall back to)
    commit(src, "root")
    a = commit(src, "A")
    sp.run(["git", "-C", src, "push", "-q", bare, a + ":refs/heads/b"], check=True)
    sp.run(["git", "-C", src, "reset", "-q", "--hard", "HEAD~1"], check=True)
    b = commit(src, "B")
    sp.run(["git", "-C", src, "push", "-q", "--force", bare, b + ":refs/heads/b"], check=True)
    # A is off b's lineage but still in src's object store

    def logtext(extra):
        return ("## 2026-10-20T10:00:00Z | P-1-INTENT | push-intent | branch=b "
                "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
                "## 2026-10-20T11:00:00Z | P-1-CONFIRM | push-confirm | branch=b remote_sha=%s "
                "credential=OWNER-PAT-1 evidence=\"git ls-remote r refs/heads/b -> %s\"\n" % (a, a)
                + extra)

    def run_with(extra_text, *extra_args):
        c = commit(src, "log", path=logrel, text=logtext(extra_text))
        sp.run(["git", "-C", src, "push", "-q", "--force", bare, c + ":refs/heads/b"],
               check=True)
        return run("evidence_sequence_check.py", logpath, "--pushlog-schema",
                   "--remote-verify", "--remote", bare, "--branch", "b",
                   "--expect-tip", c, *extra_args)

    # 1. off-lineage HISTORICAL without an evidence-back -> FAIL
    r = run_with("")
    check("AT-27", "an off-lineage historical confirmation without backing fails",
          r.returncode == 1 and "no reviewed evidence-back record supports it" in r.stdout,
          r.stdout)
    # 2. with a matching evidence-back: Round-22 requires INDEPENDENT push
    #    provenance (commit existence alone is not push evidence). The back
    #    carries a recorded push_event= citation (id/ref/before/after) plus a
    #    review_ref= that resolves into a tracked adjudicating register whose
    #    section mentions the backed entry.
    register = os.path.join(src, "docs", "decisions", "REGISTER.md")
    os.makedirs(os.path.dirname(register), exist_ok=True)
    open(register, "w").write(
        "# Test register\n\n## Part 99 -- Round-22 test adjudication\n\n"
        "The adjudicator reviewed the evidence-back P-2 for the off-lineage "
        "confirmation P-1-CONFIRM (SHA %s) and accepted the recorded "
        "PushEvent ref-transition evidence.\n" % a[:12])
    back = ("## 2026-10-20T12:00:00Z | P-2 | evidence-back | backs=P-1-CONFIRM branch=b "
            "remote_sha=%s credential=OWNER-PAT-1 review_ref=docs/decisions/REGISTER.md#Part-99 "
            "evidence=\"branch rebuilt per reviewed direction; commit %s preserved in the "
            "object store; push_event=12345678 ref=refs/heads/b before=%s after=%s "
            "created=2026-10-20T09:00:00Z\"\n" % (a, a, "0" * 40, a)
            + "## 2026-10-20T13:00:00Z | P-3-INTENT | push-intent | branch=b "
            "expected_remote_tip_before_push=%s credential=OWNER-PAT-1\n" % a)
    r = run_with(back)
    check("AT-27", "an evidence-backed off-lineage confirmation passes only with adjudicated push_event provenance",
          r.returncode == 0 and "PUSH-PROVEN" in r.stdout
          and "recorded PushEvent 12345678" in r.stdout, r.stdout)
    # 3. evidence-back for a SHA that exists nowhere -> FAIL (fail closed): the
    #    backing BINDS to the backed entry (same declared SHA) but that SHA is
    #    not a real commit -- the object-store/API proof is what must reject it
    ghost_sha = "e" * 40
    ghost_log = ("## 2026-10-20T10:00:00Z | P-4-INTENT | push-intent | branch=b "
                 "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
                 "## 2026-10-20T11:00:00Z | P-4-CONFIRM | push-confirm | branch=b remote_sha=%s "
                 "credential=OWNER-PAT-1 evidence=\"git ls-remote r refs/heads/b -> %s\"\n"
                 "## 2026-10-20T12:00:00Z | P-5 | evidence-back | backs=P-4-CONFIRM branch=b "
                 "remote_sha=%s credential=OWNER-PAT-1 evidence=\"fabricated backing\"\n"
                 % (ghost_sha, ghost_sha, ghost_sha))
    c3 = commit(src, "log3", path=logrel, text=ghost_log)
    sp.run(["git", "-C", src, "push", "-q", "--force", bare, c3 + ":refs/heads/b"], check=True)
    r = run("evidence_sequence_check.py", logpath, "--pushlog-schema", "--remote-verify",
            "--remote", bare, "--branch", "b", "--expect-tip", c3)
    check("AT-27", "an evidence-back for a nonexistent SHA fails (fail closed)",
          r.returncode == 1 and "not present in the local object store" in r.stdout,
          r.stdout)
    # 4. dangling evidence-back (backs a nonexistent entry) -> FAIL
    dangling = ("## 2026-10-20T12:00:00Z | P-4 | evidence-back | backs=P-99-CONFIRM "
                "branch=b remote_sha=%s credential=OWNER-PAT-1 evidence=\"x\"\n" % a)
    r = run_with(dangling)
    check("AT-27", "a dangling evidence-back fails",
          r.returncode == 1 and "dangling backing" in r.stdout, r.stdout)
    # 5. a cited run OUTSIDE the manifest with no API configured -> FAIL (5.6)
    man = tmpfile(json.dumps({"runs": {"111111": {"head_sha": a, "conclusion": "success"}}}),
                  ".json")
    cite = ("## 2026-10-20T12:00:00Z | P-5 | evidence-back | backs=P-1-CONFIRM branch=b "
            "remote_sha=%s credential=OWNER-PAT-1 evidence=\"run 222222 on this head\"\n" % a)
    r = run_with(cite, "--runs-manifest", man)
    check("AT-27", "a cited run outside the manifest window with no API fails the gate",
          r.returncode == 1 and "could not be verified" in r.stdout
          and "manifest omission is never an approval" in r.stdout, r.stdout)
    # 6. a cited run whose head_sha contradicts the recorded SHA -> FAIL
    bad_man = tmpfile(json.dumps({"runs": {"222222": {"head_sha": "b" * 40,
                                                     "conclusion": "success"}}}), ".json")
    r = run_with(cite, "--runs-manifest", bad_man)
    check("AT-27", "a cited run with a contradicting head_sha fails",
          r.returncode == 1 and "evidence does not support the recorded SHA" in r.stdout,
          r.stdout)
    # 7. a cited run on the exact head with a FAILED conclusion and NO push
    #    event (pull_request) is disclosed, proves neither CI success nor the
    #    push, and cannot replace the backing requirement
    fail_man = tmpfile(json.dumps({"runs": {"222222": {"head_sha": a,
                                                      "head_branch": "b",
                                                      "event": "pull_request",
                                                      "conclusion": "failure"}}}), ".json")
    cite_fail = ("## 2026-10-20T10:00:00Z | P-6-INTENT | push-intent | branch=b "
                 "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
                 "## 2026-10-20T11:00:00Z | P-6-CONFIRM | push-confirm | branch=b remote_sha=%s "
                 "credential=OWNER-PAT-1 evidence=\"git ls-remote r refs/heads/b -> %s; "
                 "gate run 222222\"\n" % (a, a))
    c6 = commit(src, "log6", path=logrel, text=cite_fail)
    sp.run(["git", "-C", src, "push", "-q", "--force", bare, c6 + ":refs/heads/b"], check=True)
    r = run("evidence_sequence_check.py", logpath, "--pushlog-schema", "--remote-verify",
            "--remote", bare, "--branch", "b", "--expect-tip", c6,
            "--runs-manifest", fail_man)
    check("AT-27", "a failed pull_request-event run proves neither CI success nor the push",
          r.returncode == 1 and "successful gate evidence" in r.stdout
          and "commit existence is not push provenance" in r.stdout, r.stdout)
    # 8. a push-EVENT run (event=push, head_branch, exact head) with a FAILED
    #    conclusion proves the PUSH but never CI success -- both facts are
    #    reported separately (Round-22 claim separation)
    pushfail_man = tmpfile(json.dumps({"runs": {"222222": {"head_sha": a,
                                                        "head_branch": "b",
                                                        "event": "push",
                                                        "conclusion": "failure"}}}), ".json")
    both = cite_fail + ("## 2026-10-20T12:00:00Z | P-7 | evidence-back | backs=P-6-CONFIRM "
                        "branch=b remote_sha=%s credential=OWNER-PAT-1 "
                        "review_ref=docs/decisions/REGISTER.md#Part-99 evidence=\"run 222222 "
                        "connects the head to branch b (event=push; conclusion=failure "
                        "disclosed separately)\"\n" % a
                        + "## 2026-10-20T13:00:00Z | P-8-INTENT | push-intent | branch=b "
                        "expected_remote_tip_before_push=%s credential=OWNER-PAT-1\n" % a)
    c7 = commit(src, "log7", path=logrel, text=both)
    sp.run(["git", "-C", src, "push", "-q", "--force", bare, c7 + ":refs/heads/b"], check=True)
    r = run("evidence_sequence_check.py", logpath, "--pushlog-schema", "--remote-verify",
            "--remote", bare, "--branch", "b", "--expect-tip", c7,
            "--runs-manifest", pushfail_man)
    check("AT-27", "a failed push-event run proves the push event but never CI success (claims separated)",
          r.returncode == 0 and "PUSH-PROVEN by push-event run 222222" in r.stdout
          and "none of the cited runs succeeded" in r.stdout, r.stdout)
    # 9. a SUCCESSFUL push-event run on the exact head -> PUSH-PROVEN with
    #    CI-success support
    ok_man = tmpfile(json.dumps({"runs": {"222222": {"head_sha": a, "head_branch": "b",
                                                    "event": "push",
                                                    "conclusion": "success"}}}), ".json")
    r = run("evidence_sequence_check.py", logpath, "--pushlog-schema", "--remote-verify",
            "--remote", bare, "--branch", "b", "--expect-tip", c7,
            "--runs-manifest", ok_man)
    check("AT-27", "a verified successful push-event run on the exact head proves the branch push",
          r.returncode == 0 and "PUSH-PROVEN by push-event run 222222" in r.stdout
          and "CI-success: 1 verified successful run" in r.stdout, r.stdout)


# ------------------------------ AT-28 (protected-main merge, Round-21)

def test_at28_merge_transition():
    """Round-21 direction 6: a push to the protected main is declared ONLY by a
    verified merge-intent -- the candidate must be a merge commit of the
    declared pre-merge base and the declared approved PR head, the intent must
    ride the PR branch, and the PR must be merged with an approving non-author
    review. Raced/unrelated updates and unapproved merges FAIL."""
    import subprocess as sp
    d = tempfile.mkdtemp()
    bare = os.path.join(d, "remote.git")
    sp.run(["git", "init", "-q", "--bare", bare], check=True)
    src = os.path.join(d, "src")
    sp.run(["git", "init", "-q", "-b", "main", src], check=True)

    def commit(tag, path=None, text=None):
        if path and text is not None:
            os.makedirs(os.path.dirname(os.path.join(src, path)), exist_ok=True)
            open(os.path.join(src, path), "w").write(text)
            sp.run(["git", "-C", src, "add", "docs"], check=True)
        sp.run(["git", "-C", src, "-c", "user.email=t@t", "-c", "user.name=t",
                "commit", "-q", "--allow-empty", "-m", tag], check=True)
        return sp.run(["git", "-C", src, "rev-parse", "HEAD"],
                      capture_output=True, text=True).stdout.strip()

    def push(refspec, force=False):
        sp.run(["git", "-C", src, "push", "-q"] + (["--force"] if force else []) +
               [bare, refspec], check=True)

    logrel = "docs/evidence/push-log/PUSH_LOG.md"
    logpath = os.path.join(src, logrel)

    # ---- build: main@base, feat@head2 carrying P-1-INTENT + P-2 merge-intent
    base = commit("base")
    push("HEAD:refs/heads/main")
    sp.run(["git", "-C", src, "checkout", "-q", "-b", "feat"], check=True)
    # expected_head=DEFERRED: the head is the commit that introduces the
    # merge-intent itself (a record cannot embed its own commit's SHA; the
    # checker resolves DEFERRED via the introducing commit)
    mi = ("## 2026-10-20T10:00:00Z | P-1-INTENT | push-intent | branch=feat "
          "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
          "## 2026-10-20T11:00:00Z | P-2 | merge-intent | target=main pr=7 "
          "source_branch=feat expected_pre_merge_base=%s expected_head=DEFERRED "
          "credential=OWNER-PAT-1 authorized_by=\"auditor direction 6\"\n" % base)
    head2 = commit("merge-intent", path=logrel, text=mi)
    push("HEAD:refs/heads/feat")

    def chk(branch, tip, pr_file=None, rev_file=None):
        extra = []
        if pr_file:
            extra += ["--pr-json", pr_file]
        if rev_file:
            extra += ["--reviews-json", rev_file]
        return run("evidence_sequence_check.py", logpath, "--pushlog-schema",
                   "--remote-verify", "--remote", bare, "--branch", branch,
                   "--expect-tip", tip, *extra)

    # ---- phase 1: PRE-MERGE, on the source branch (main still at base)
    r = chk("feat", head2)
    check("AT-28", "on the source branch the merge-intent binds the exact candidate and base",
          r.returncode == 0 and "future reviewed merge" in r.stdout, r.stdout)
    # a stale merge-intent (candidate advanced past the DEFERRED head) fails
    head3 = commit("after-intent", path=logrel, text=mi)
    push("HEAD:refs/heads/feat", force=True)
    r = chk("feat", head3)
    check("AT-28", "a stale merge-intent on an advanced candidate fails (redeclare before merging)",
          r.returncode == 1 and "redeclare the merge-intent" in r.stdout, r.stdout)
    sp.run(["git", "-C", src, "reset", "-q", "--hard", head2], check=True)
    push("HEAD:refs/heads/feat", force=True)

    # ---- phase 2: the reviewed merge lands on main
    sp.run(["git", "-C", src, "checkout", "-q", "main"], check=True)
    sp.run(["git", "-C", src, "-c", "user.email=t@t", "-c", "user.name=t",
            "merge", "-q", "--no-ff", "-m", "merge pr 7", "feat"], check=True)
    m = sp.run(["git", "-C", src, "rev-parse", "HEAD"],
               capture_output=True, text=True).stdout.strip()
    push("HEAD:refs/heads/main")
    pr_ok = tmpfile(json.dumps({"merged": True, "merge_commit_sha": m,
                                "merged_at": "2026-10-20T12:00:00Z",
                                "user": {"login": "author"},
                                "head": {"sha": head2, "ref": "feat"},
                                "base": {"ref": "main"}}), ".json")
    rev_ok = tmpfile(json.dumps([{"state": "APPROVED",
                                  "user": {"login": "reviewer"},
                                  "commit_id": head2,
                                  "submitted_at": "2026-10-20T11:30:00Z"}]), ".json")
    rev_author = tmpfile(json.dumps([{"state": "APPROVED",
                                      "user": {"login": "author"},
                                      "commit_id": head2,
                                      "submitted_at": "2026-10-20T11:30:00Z"}]), ".json")
    rev_none = tmpfile(json.dumps([{"state": "COMMENTED",
                                    "user": {"login": "reviewer"},
                                    "commit_id": head2,
                                    "submitted_at": "2026-10-20T11:30:00Z"}]), ".json")
    pr_open = tmpfile(json.dumps({"merged": False, "merge_commit_sha": None,
                                  "merged_at": None,
                                  "user": {"login": "author"},
                                  "head": {"sha": head2, "ref": "feat"},
                                  "base": {"ref": "main"}}), ".json")
    r = chk("main", m, pr_ok, rev_ok)
    check("AT-28", "a legitimate reviewed-PR merge of the declared base+head passes on main",
          r.returncode == 0 and "MERGE-VERIFIED" in r.stdout
          and "DECLARED" in r.stdout, r.stdout)
    r = chk("main", m, pr_ok, rev_none)
    check("AT-28", "a merge without an approving review fails",
          r.returncode == 1 and "no approving review" in r.stdout, r.stdout)
    r = chk("main", m, pr_ok, rev_author)
    check("AT-28", "an approval from the PR author does not satisfy independence",
          r.returncode == 1 and "no approving review" in r.stdout, r.stdout)
    r = chk("main", m, pr_open, rev_ok)
    check("AT-28", "a merge-intent whose PR is not merged fails",
          r.returncode == 1 and "is not merged" in r.stdout, r.stdout)

    # ---- phase 3: raced main (a commit lands after the declared base)
    sp.run(["git", "-C", src, "reset", "-q", "--hard", base], check=True)
    commit("racer")
    push("HEAD:refs/heads/main", force=True)
    sp.run(["git", "-C", src, "-c", "user.email=t@t", "-c", "user.name=t",
            "merge", "-q", "--no-ff", "-m", "merge pr 7 (raced)", "feat"], check=True)
    m2 = sp.run(["git", "-C", src, "rev-parse", "HEAD"],
                capture_output=True, text=True).stdout.strip()
    push("HEAD:refs/heads/main", force=True)
    r = chk("main", m2, pr_ok, rev_ok)
    check("AT-28", "a merge onto a raced main (first parent != declared base) fails",
          r.returncode == 1 and "not the declared pre-merge base" in r.stdout, r.stdout)

    # ---- phase 4: an unrelated direct (non-merge) push to main. The pushed log
    #    carries the merge-intent content (a protected main always carries the
    #    log; the earlier fixture dropped it, crashing the checker on a missing
    #    file instead of exercising the non-merge rejection)
    sp.run(["git", "-C", src, "reset", "-q", "--hard", base], check=True)
    direct = commit("direct push", path=logrel, text=mi)
    push("HEAD:refs/heads/main", force=True)
    r = chk("main", direct, pr_ok, rev_ok)
    check("AT-28", "a direct non-merge push to main fails",
          r.returncode == 1 and "not a merge commit" in r.stdout, r.stdout)

def test_at29_path_filter_liveness():
    """Round-21 direction 7: an ACTIVE workflow reference must declare ~governs=,
    and the live workflow's paths/paths-ignore filters must not exclude the
    governed changes; filter values of an invalid type FAIL."""
    d = tempfile.mkdtemp()
    os.makedirs(os.path.join(d, ".github", "workflows"), exist_ok=True)
    wf = os.path.join(d, ".github", "workflows", "gate.yml")
    reg = {"tests": {}}
    open(os.path.join(d, "registry.json"), "w").write(json.dumps(reg))
    os.makedirs(os.path.join(d, "docs", "directives"), exist_ok=True)

    def write_wf(text):
        open(wf, "w").write(text)

    def idx(enf):
        open(os.path.join(d, "docs", "directives", "INDEX.md"), "w").write(
            "# index\n<!-- BEGIN RULE INVENTORY -->\n"
            "reading-note-path: docs/note.md\n"
            "pinned-blobs: b6bb35e84956d84a8c9d9e65734c212baf3ed075\n"
            "| Rule ID | Status | Enforcement point(s) | Acceptance test(s) |\n"
            "|---|---|---|---|\n"
            "| EC-X | ACTIVE | " + enf + " | NONE |\n"
            "<!-- END RULE INVENTORY -->\n")
        open(os.path.join(d, "docs", "note.md"), "w").write(
            "blob b6bb35e84956d84a8c9d9e65734c212baf3ed075 read\n| EC-X | -> mapped |\n")

    def run_idx():
        return run("reading_note_map_check.py", "--index", "docs/directives/INDEX.md",
                   "--repo-root", d, "--registry", "registry.json")

    job = "jobs:\n  gate-job:\n    name: Gate Job\n    runs-on: ubuntu-latest\n"
    base_on = "on:\n  push:\n    branches: ['**']\n  pull_request:\n    branches: ['**']\n"
    # 1. ACTIVE reference without ~governs= fails
    write_wf("name: gate\n" + base_on + job)
    idx("workflow:.github/workflows/gate.yml#gate-job~events=push+pull_request")
    r = run_idx()
    check("AT-29", "ACTIVE workflow reference without ~governs fails",
          r.returncode == 1 and "must declare the changes it governs" in r.stdout, r.stdout)
    # 2. paths filter excluding the governed path fails
    write_wf("name: gate\non:\n  push:\n    branches: ['**']\n    paths: [docs/other/**]\n"
             "  pull_request:\n    branches: ['**']\n" + job)
    idx("workflow:.github/workflows/gate.yml#gate-job~events=push+pull_request"
        "~governs=docs/evidence/push-log/**")
    r = run_idx()
    check("AT-29", "a paths filter excluding the governed path fails",
          r.returncode == 1 and "exclude the governed path" in r.stdout, r.stdout)
    # 3. paths-ignore covering the governed path fails
    write_wf("name: gate\non:\n  push:\n    branches: ['**']\n"
             "    paths-ignore: [docs/evidence/**]\n  pull_request:\n    branches: ['**']\n" + job)
    idx("workflow:.github/workflows/gate.yml#gate-job~events=push+pull_request"
        "~governs=docs/evidence/push-log/PUSH_LOG.md")
    r = run_idx()
    check("AT-29", "a paths-ignore filter covering the governed path fails",
          r.returncode == 1 and "paths-ignore" in r.stdout
          and "excludes the governed path" in r.stdout, r.stdout)
    # 4. covering paths filter passes
    write_wf("name: gate\non:\n  push:\n    branches: ['**']\n"
             "    paths: [docs/evidence/**, scripts/**]\n  pull_request:\n    branches: ['**']\n"
             + job)
    idx("workflow:.github/workflows/gate.yml#gate-job~events=push+pull_request"
        "~governs=docs/evidence/push-log/**+scripts/ci/evidence_sequence_check.py")
    r = run_idx()
    check("AT-29", "a paths filter covering every governed path passes",
          r.returncode == 0, r.stdout)
    # 5. no path filters at all passes
    write_wf("name: gate\n" + base_on + job)
    idx("workflow:.github/workflows/gate.yml#gate-job~events=push+pull_request"
        "~governs=docs/evidence/push-log/**")
    r = run_idx()
    check("AT-29", "an unfiltered workflow passes for any governed set",
          r.returncode == 0, r.stdout)
    # 6. tree-wide governance (~governs=**) with a restrictive paths filter fails
    write_wf("name: gate\non:\n  push:\n    branches: ['**']\n    paths: [docs/**]\n"
             "  pull_request:\n    branches: ['**']\n" + job)
    idx("workflow:.github/workflows/gate.yml#gate-job~events=push+pull_request~governs=**")
    r = run_idx()
    check("AT-29", "tree-wide governance cannot be satisfied by a restrictive paths filter",
          r.returncode == 1 and "governs the whole tree" in r.stdout, r.stdout)
    # 7. INVALID filter types fail even when they would otherwise pass
    write_wf("name: gate\non:\n  push:\n    branches: {wrong: type}\n"
             "  pull_request:\n    branches: ['**']\n" + job)
    idx("workflow:.github/workflows/gate.yml#gate-job~events=push+pull_request~governs=**")
    r = run_idx()
    check("AT-29", "an invalid branch-filter type fails",
          r.returncode == 1 and "INVALID filter type" in r.stdout, r.stdout)
    write_wf("name: gate\non:\n  push:\n    branches: ['**']\n    paths: 123\n"
             "  pull_request:\n    branches: ['**']\n" + job)
    idx("workflow:.github/workflows/gate.yml#gate-job~events=push+pull_request~governs=**")
    r = run_idx()
    check("AT-29", "an invalid paths-filter type fails",
          r.returncode == 1 and "INVALID filter type" in r.stdout, r.stdout)
    # 8. unknown annotation key fails closed
    idx("workflow:.github/workflows/gate.yml#gate-job~events=push+pull_request"
        "~governs=**~unexpected=xyz")
    r = run_idx()
    check("AT-29", "an unknown annotation key fails closed",
          r.returncode == 1 and "unknown annotation key" in r.stdout, r.stdout)


# ------------------- AT-30 (EC-23 report-field leakage, Round-22 direction 2)

def test_at30_secret_control_field_leakage():
    """Round-22 direction 2: NO report-derived value may ever reach stdout or
    stderr. A malicious or malformed report can place the synthetic key in
    RuleID, File, StartLine, an unexpected field, a non-object array entry, or
    malformed JSON bytes; every such report must fail closed WITHOUT printing
    the key, on positive and negative paths alike."""
    KEY = "AKIA" + "ZZR22FIELDPROBEQQ9"
    d = tempfile.mkdtemp()

    def write_raw(name, text):
        p = os.path.join(d, name)
        open(p, "w").write(text)
        return p

    def ctrl(mode, *args):
        env = dict(os.environ)
        env["SEEDED_KEY"] = KEY
        return subprocess.run(
            [PY, os.path.join(CI, "secret_control_check.py"), mode, *args],
            capture_output=True, text=True, env=env, shell=False)

    good = [{"RuleID": "aws-access-token", "File": "fixture.txt",
             "StartLine": 1, "Secret": "REDACTED", "Match": "REDACTED"}]
    clean_log = write_raw("clean.log", "no key material here\n")

    def no_leak(r):
        return KEY not in r.stdout and KEY not in r.stderr

    # 1. an UNEXPECTED finding whose RuleID carries the key: count-only
    #    diagnostics; the key never reaches the output
    r1 = write_raw("r1.json", json.dumps(
        good + [{"RuleID": KEY, "File": "other.txt", "StartLine": 3}]))
    r = ctrl("positive", "--report", r1, "--exit-code", "101")
    check("AT-30", "unexpected finding with the key in RuleID fails without printing it",
          r.returncode == 1 and "unexpected additional finding(s): 1" in r.stdout
          and no_leak(r), r.stdout + r.stderr)
    # 2. the key in the File field of an unexpected finding
    r2 = write_raw("r2.json", json.dumps(
        good + [{"RuleID": "generic-api-key", "File": KEY, "StartLine": 4}]))
    r = ctrl("positive", "--report", r2, "--exit-code", "101")
    check("AT-30", "unexpected finding with the key in File fails without printing it",
          r.returncode == 1 and no_leak(r), r.stdout + r.stderr)
    # 3. the key as a STRING StartLine: type validation fails closed
    r3 = write_raw("r3.json", json.dumps(
        [{"RuleID": "aws-access-token", "File": "fixture.txt",
          "StartLine": KEY, "Secret": "REDACTED", "Match": "REDACTED"}]))
    r = ctrl("positive", "--report", r3, "--exit-code", "101")
    check("AT-30", "a string StartLine carrying the key fails type validation without printing it",
          r.returncode == 1 and "'StartLine' is not a" in r.stdout
          and "contents withheld" in r.stdout and no_leak(r), r.stdout + r.stderr)
    # 4. non-object array entries (the key embedded in a string row)
    r4 = write_raw("r4.json", json.dumps(
        good + ["row carrying " + KEY, 5, None]))
    r = ctrl("positive", "--report", r4, "--exit-code", "101")
    check("AT-30", "non-object report rows fail closed without printing row contents",
          r.returncode == 1 and "is not a JSON object" in r.stdout
          and "contents withheld" in r.stdout and no_leak(r), r.stdout + r.stderr)
    # 5. malformed JSON with the key in the raw bytes
    r5 = write_raw("r5.json", "[{\"RuleID\": \"" + KEY + "\", oops")
    r = ctrl("positive", "--report", r5, "--exit-code", "101")
    check("AT-30", "malformed JSON carrying the key fails closed without printing it",
          r.returncode == 1 and "not parseable JSON" in r.stdout
          and "redaction FAILED" in r.stdout and no_leak(r), r.stdout + r.stderr)
    # 6. non-array JSON (dict) with the key inside a value
    r6 = write_raw("r6.json", json.dumps({"RuleID": KEY}))
    r = ctrl("positive", "--report", r6, "--exit-code", "101")
    check("AT-30", "a non-array report fails closed without printing values",
          r.returncode == 1 and "not a JSON array" in r.stdout and no_leak(r),
          r.stdout + r.stderr)
    # 7. missing report fails; requested-but-missing log fails closed
    r = ctrl("positive", "--report", os.path.join(d, "absent.json"),
             "--exit-code", "101")
    check("AT-30", "a missing report fails without printing anything derived",
          r.returncode == 1 and "absent" in r.stdout and no_leak(r),
          r.stdout + r.stderr)
    rok = write_raw("rok.json", json.dumps(good))
    r = ctrl("positive", "--report", rok, "--exit-code", "101",
             "--log-file", os.path.join(d, "missing.log"))
    check("AT-30", "a requested-but-missing capture log fails closed",
          r.returncode == 1 and "MISSING" in r.stdout, r.stdout)
    # 8. negative mode: findings present -> count-only message, no values
    r = ctrl("negative", "--report", r2, "--exit-code", "0")
    check("AT-30", "negative mode reports only the finding count (never rule/file values)",
          r.returncode == 1 and "produced 2 finding(s)" in r.stdout
          and "never printed" in r.stdout and no_leak(r), r.stdout + r.stderr)
    # 9. negative mode with a malformed row value fails closed
    r = ctrl("negative", "--report", r3, "--exit-code", "0")
    check("AT-30", "negative mode type-validates rows and fails closed",
          r.returncode == 1 and "'StartLine' is not a" in r.stdout
          and no_leak(r), r.stdout + r.stderr)
    # 10. clean paths still pass with the key nowhere in the output
    r = ctrl("positive", "--report", rok, "--exit-code", "101",
             "--log-file", clean_log)
    check("AT-30", "the hardened positive control passes without printing the key",
          r.returncode == 0 and no_leak(r), r.stdout + r.stderr)
    rclean = write_raw("rclean.json", "[]")
    r = ctrl("negative", "--report", rclean, "--exit-code", "0",
             "--log-file", clean_log)
    check("AT-30", "the hardened negative control passes without printing the key",
          r.returncode == 0 and no_leak(r), r.stdout + r.stderr)
    # 11. scanner execution error masquerading as detection still fails
    r = ctrl("positive", "--report", rok, "--exit-code", "1")
    check("AT-30", "a scanner exit 1 is rejected (execution error, not detection)",
          r.returncode == 1 and "masquerading" in r.stdout, r.stdout)


# --------------- AT-31 (EC-24 push provenance + pre-push binding, Round-22)

def test_at31_push_provenance_and_binding():
    """Round-22 directions 3-4: an off-lineage historical confirmation must
    prove the BRANCH PUSH (push-event run / PushEvent / adjudicated citation),
    never merely commit existence; a push-not-proven reclassification is the
    honest disclosure; and the terminal intent's declared pre-push state is
    bound to the trusted PushEvent before/after state (false ABSENT, prior
    unlogged ref updates, and raced updates all fail)."""
    import subprocess as sp
    d = tempfile.mkdtemp()
    bare = os.path.join(d, "remote.git")
    sp.run(["git", "init", "-q", "--bare", bare], check=True)
    src = os.path.join(d, "src")
    sp.run(["git", "init", "-q", "-b", "main", src], check=True)

    def commit(repo, tag, path=None, text=None, add="docs"):
        if path and text is not None:
            os.makedirs(os.path.dirname(os.path.join(repo, path)), exist_ok=True)
            open(os.path.join(repo, path), "w").write(text)
            sp.run(["git", "-C", repo, "add", add], check=True)
        sp.run(["git", "-C", repo, "-c", "user.email=t@t", "-c", "user.name=t",
                "commit", "-q", "--allow-empty", "-m", tag], check=True)
        return sp.run(["git", "-C", repo, "rev-parse", "HEAD"],
                      capture_output=True, text=True).stdout.strip()

    logrel = "docs/evidence/push-log/PUSH_LOG.md"
    logpath = os.path.join(src, logrel)
    commit(src, "root")
    a = commit(src, "A")
    sp.run(["git", "-C", src, "push", "-q", bare, a + ":refs/heads/b"], check=True)
    sp.run(["git", "-C", src, "reset", "-q", "--hard", "HEAD~1"], check=True)
    b2 = commit(src, "B")
    sp.run(["git", "-C", src, "push", "-q", "--force", bare, b2 + ":refs/heads/b"],
           check=True)
    # adjudication-register fixture for review_ref resolution
    reg = os.path.join(src, "docs", "decisions", "REGISTER.md")
    os.makedirs(os.path.dirname(reg), exist_ok=True)
    open(reg, "w").write(
        "# Test register\n\n## Part 98 -- Round-22 provenance adjudication\n\n"
        "The adjudicator reviewed the evidence-backs P-2/P-20 backing "
        "P-1-CONFIRM (SHA %s) and the reclassification record P-21.\n" % a[:12])
    sp.run(["git", "-C", src, "add", "docs"], check=True)

    base_log = ("## 2026-10-20T10:00:00Z | P-1-INTENT | push-intent | branch=b "
                "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
                "## 2026-10-20T11:00:00Z | P-1-CONFIRM | push-confirm | branch=b remote_sha=%s "
                "credential=OWNER-PAT-1 evidence=\"git ls-remote r refs/heads/b -> %s\"\n" % (a, a))
    empty_man = tmpfile(json.dumps({"runs": {}}), ".json")
    _seq = [0]

    def scenario(extra, *args, cand_events=None, extra_events=()):
        # every scenario's candidate is declared by a fresh terminal intent
        # (the two-phase protocol) with a UNIQUE id per scenario (the
        # introducing-commit binding uses git log -S, so a reused id would
        # resolve to the first scenario that introduced it); the pre-push tip
        # is the newest recorded state for b (P-1-CONFIRM's remote_sha = a)
        _seq[0] += 1
        text = base_log + extra + (
            "## 2026-10-20T14:00:00Z | P-9%d-INTENT | push-intent | branch=b "
            "expected_remote_tip_before_push=%s credential=OWNER-PAT-1\n"
            % (_seq[0], a))
        c = commit(src, "log", path=logrel, text=text)
        sp.run(["git", "-C", src, "push", "-q", "--force", bare, c + ":refs/heads/b"],
               check=True)
        args = list(args)
        if cand_events is not None or extra_events:
            events = list(extra_events) + [dict(e, after=c) for e in (cand_events or [])]
            args += ["--push-events-json", tmpfile(json.dumps(events), ".json")]
        return run("evidence_sequence_check.py", logpath, "--pushlog-schema",
                   "--remote-verify", "--remote", bare, "--branch", "b",
                   "--expect-tip", c, *args)

    # 1. a real but never-proven-pushed commit: back proves existence only
    narrative = ("## 2026-10-20T12:00:00Z | P-2 | evidence-back | backs=P-1-CONFIRM "
                 "branch=b remote_sha=%s credential=OWNER-PAT-1 evidence=\"branch was "
                 "rebuilt per a reviewed direction; the commit is preserved in the "
                 "object store and the rebuild was documented in transcripts\"\n" % a)
    r = scenario(narrative, "--runs-manifest", empty_man)
    check("AT-31", "a rebuild narrative without ref-transition evidence is not push provenance",
          r.returncode == 1 and "proves only COMMIT EXISTENCE" in r.stdout, r.stdout)
    # 2. a push-event run on the WRONG branch contradicts the entry
    wrong_man = tmpfile(json.dumps({"runs": {"111111": {
        "head_sha": a, "head_branch": "other", "event": "push",
        "conclusion": "success"}}}), ".json")
    cite_other = ("## 2026-10-20T12:00:00Z | P-2 | evidence-back | backs=P-1-CONFIRM "
                  "branch=b remote_sha=%s credential=OWNER-PAT-1 evidence=\"run 111111 "
                  "on this head\"\n" % a)
    r = scenario(cite_other, "--runs-manifest", wrong_man)
    check("AT-31", "a push-event run on another branch fails the entry",
          r.returncode == 1 and "push event on another branch" in r.stdout, r.stdout)
    # 3. a non-push event (pull_request) never proves the push
    pr_man = tmpfile(json.dumps({"runs": {"111111": {
        "head_sha": a, "head_branch": "b", "event": "pull_request",
        "conclusion": "success"}}}), ".json")
    r = scenario(cite_other, "--runs-manifest", pr_man)
    check("AT-31", "a pull_request-event run does not prove the branch push",
          r.returncode == 1 and "proves only COMMIT EXISTENCE" in r.stdout, r.stdout)
    # 4. a cited run absent from the manifest with no API fails the gate
    r = scenario(cite_other.replace("111111", "999999"), "--runs-manifest", empty_man)
    check("AT-31", "a cited run that cannot be verified fails (missing API proof)",
          r.returncode == 1 and "could not be verified" in r.stdout
          and "manifest omission is never an approval" in r.stdout, r.stdout)
    # 5. push-not-proven reclassification with a VALID review_ref discloses honestly
    reclass = ("## 2026-10-20T12:00:00Z | P-21 | push-not-proven | reclassifies=P-1-CONFIRM "
               "branch=b remote_sha=%s credential=OWNER-PAT-1 "
               "review_ref=docs/decisions/REGISTER.md#Part-98 evidence=\"no "
               "push-event run exists for this head; the Events-API window no "
               "longer retains the transition; the branch push is NOT PROVEN and "
               "is never counted as evidence\"\n" % a)
    r = scenario(reclass, "--runs-manifest", empty_man)
    check("AT-31", "a push-not-proven reclassification with a valid review_ref is accepted",
          r.returncode == 0 and "NOT-PROVEN" in r.stdout
          and "RECLASSIFIED by P-21" in r.stdout, r.stdout)
    # 6. a reclassification whose review_ref does not mention the entry fails
    bad_reg = os.path.join(src, "docs", "decisions", "OTHER.md")
    open(bad_reg, "w").write("# Other\n\n## Part 50 -- unrelated\n\nnothing here\n")
    reclass_bad = reclass.replace("REGISTER.md#Part-98", "OTHER.md#Part-50")
    r = scenario(reclass_bad, "--runs-manifest", empty_man)
    check("AT-31", "a review_ref whose section does not mention the entry fails",
          r.returncode == 1 and "does not mention the backed entry" in r.stdout,
          r.stdout)
    # 7. duplicate reclassification of the same confirmation fails
    r = scenario(reclass + reclass.replace("P-21", "P-22"), "--runs-manifest", empty_man)
    check("AT-31", "a second reclassification of the same confirmation fails",
          r.returncode == 1 and "a second time" in r.stdout, r.stdout)
    # 8. a live Events-API PushEvent + review_ref proves the push (the
    #    candidate's own push event is provided alongside the historical one)
    back_ev = ("## 2026-10-20T12:00:00Z | P-20 | evidence-back | backs=P-1-CONFIRM "
               "branch=b remote_sha=%s credential=OWNER-PAT-1 "
               "review_ref=docs/decisions/REGISTER.md#Part-98 evidence=\"PushEvent "
               "observed live\"\n" % a)
    ev_for_a8 = {"id": "777001", "ref": "refs/heads/b", "before": "0" * 40,
                 "after": a, "created_at": "2026-10-20T09:00:00Z"}
    r = scenario(back_ev, "--runs-manifest", empty_man,
                 cand_events=[{"id": "777040", "ref": "refs/heads/b",
                               "before": a,
                               "created_at": "2026-10-20T14:00:00Z"}],
                 extra_events=[ev_for_a8])
    check("AT-31", "a live PushEvent with an adjudicating review_ref proves the push",
          r.returncode == 0 and "live Events-API PushEvent 777001" in r.stdout
          and "PUSH-PROVEN" in r.stdout, r.stdout)
    # 9. a recorded push_event citation on the WRONG ref fails
    bad_ref = back_ev.replace(
        "evidence=\"PushEvent observed live\"",
        "evidence=\"push_event=777002 ref=refs/heads/c before=%s after=%s "
        "created=2026-10-20T09:00:00Z\"" % ("0" * 40, a))
    r = scenario(bad_ref, "--runs-manifest", empty_man)
    check("AT-31", "a recorded push_event citation on the wrong ref fails",
          r.returncode == 1 and "mismatched ref-transition evidence" in r.stdout,
          r.stdout)
    # 10. a recorded push_event citation WITHOUT review_ref fails
    no_rr = bad_ref.replace("ref=refs/heads/c", "ref=refs/heads/b").replace(
        "review_ref=docs/decisions/REGISTER.md#Part-98 ", "")
    r = scenario(no_rr, "--runs-manifest", empty_man)
    check("AT-31", "recorded push_event evidence without review_ref is not adjudicated",
          r.returncode == 1 and "carries no review_ref" in r.stdout, r.stdout)

    # ---- pre-push-state binding (direction 4): PushEvent before/after.
    # Every binding scenario carries the adjudicated back for the off-lineage
    # P-1-CONFIRM (otherwise its missing provenance would fail the run for an
    # unrelated reason); the candidate's own event is injected with the real
    # candidate SHA, which cannot be known before the commit is made.
    back_ev_rec = ("## 2026-10-20T12:00:00Z | P-20 | evidence-back | backs=P-1-CONFIRM "
                   "branch=b remote_sha=%s credential=OWNER-PAT-1 "
                   "review_ref=docs/decisions/REGISTER.md#Part-98 evidence=\"push_event="
                   "777000 ref=refs/heads/b before=%s after=%s "
                   "created=2026-10-20T09:00:00Z\"\n" % (a, "0" * 40, a))
    ev_for_a = {"id": "777000", "ref": "refs/heads/b", "before": "0" * 40,
                "after": a, "created_at": "2026-10-20T09:00:00Z"}

    def binding_scenario(intent_text, cand_events, branch="b", prelude=None,
                         extra_events=()):
        _seq[0] += 1
        intent_text = intent_text.replace(
            "P-30-INTENT", "P-3%d-INTENT" % _seq[0]).replace(
            "P-31-INTENT", "P-3%d-INTENT" % _seq[0])
        text = base_log + (prelude if prelude is not None else back_ev_rec) + intent_text
        c = commit(src, "logB", path=logrel, text=text)
        sp.run(["git", "-C", src, "push", "-q", "--force", bare,
                c + ":refs/heads/" + branch], check=True)
        events = list(extra_events) + [dict(e, after=c) for e in cand_events]
        ev_file = tmpfile(json.dumps(events), ".json")
        return run("evidence_sequence_check.py", logpath, "--pushlog-schema",
                   "--remote-verify", "--remote", bare, "--branch", branch,
                   "--expect-tip", c, "--push-events-json", ev_file)

    term = ("## 2026-10-20T13:00:00Z | P-30-INTENT | push-intent | branch=b "
            "expected_remote_tip_before_push=%s credential=OWNER-PAT-1\n" % a)
    # 11. declared tip matches the PushEvent before -> BOUND
    r = binding_scenario(term, [{"id": "777010", "ref": "refs/heads/b",
                                 "before": a,
                                 "created_at": "2026-10-20T14:00:00Z"}],
                         extra_events=[ev_for_a])
    check("AT-31", "the terminal intent's declared pre-push tip is bound to the PushEvent before",
          r.returncode == 0 and "BOUND: the PushEvent verifies the ref transition" in r.stdout,
          r.stdout)
    # 12. a prior UNLOGGED ref update (before != declared) fails
    r = binding_scenario(term, [{"id": "777011", "ref": "refs/heads/b",
                                 "before": "f" * 40,
                                 "created_at": "2026-10-20T14:00:00Z"}],
                         extra_events=[ev_for_a])
    check("AT-31", "a prior unlogged ref update (PushEvent before != declared tip) fails",
          r.returncode == 1 and "prior unlogged ref update" in r.stdout, r.stdout)
    # 13. no PushEvent for the candidate fails closed (raced or absent)
    r = binding_scenario(term, [], extra_events=[ev_for_a])
    check("AT-31", "absent push-event evidence for the candidate fails closed",
          r.returncode == 1 and "no push-event evidence exists for the candidate" in r.stdout,
          r.stdout)
    # 14. false ABSENT: on a FRESH branch the intent declares ABSENT but the
    #     PushEvent records a prior (nonzero) tip
    term_absent = ("## 2026-10-20T13:00:00Z | P-31-INTENT | push-intent | branch=fresh "
                   "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n")
    r = binding_scenario(term_absent, [{"id": "777012", "ref": "refs/heads/fresh",
                                        "before": "e" * 40,
                                        "created_at": "2026-10-20T14:00:00Z"}],
                         branch="fresh")
    check("AT-31", "a false ABSENT declaration (nonzero PushEvent before) fails",
          r.returncode == 1 and "false ABSENT" in r.stdout, r.stdout)
    # 15. ABSENT with the all-zero before (newly created branch) is bound
    r = binding_scenario(term_absent, [{"id": "777013", "ref": "refs/heads/fresh",
                                        "before": "0" * 40,
                                        "created_at": "2026-10-20T14:00:00Z"}],
                         branch="fresh")
    check("AT-31", "an ABSENT declaration backed by the all-zero before is bound",
          r.returncode == 0
          and "the PushEvent verifies the ref transition" in r.stdout
          and "before=000000000000" in r.stdout, r.stdout)

    # 16-18. PRIOR-TIP citations (Round-22): a later push on the same ref
    # whose before field records the backed SHA proves the ref pointed there
    # (used when the push's own event rolled off and the API elides after)
    def prior_tip_scenario(back_text, events, prelude=None):
        _seq[0] += 1
        text = base_log + (prelude if prelude is not None else "") + back_text + (
            "## 2026-10-20T14:00:00Z | P-8%d-INTENT | push-intent | branch=b "
            "expected_remote_tip_before_push=%s credential=OWNER-PAT-1\n"
            % (_seq[0], a))
        c = commit(src, "logC", path=logrel, text=text)
        sp.run(["git", "-C", src, "push", "-q", "--force", bare, c + ":refs/heads/b"],
               check=True)
        # the candidate's own push event carries after == the candidate SHA,
        # which is only known after the commit; the sentinel is replaced here
        prepared = [dict(e, after=c) if e.get("after") == "@@CAND@@" else e
                    for e in events]
        ev_file = tmpfile(json.dumps(prepared), ".json")
        return run("evidence_sequence_check.py", logpath, "--pushlog-schema",
                   "--remote-verify", "--remote", bare, "--branch", "b",
                   "--expect-tip", c, "--runs-manifest", empty_man,
                   "--push-events-json", ev_file)

    prior_ev = {"id": "777100", "ref": "refs/heads/b", "before": a,
                "after": None, "created_at": "2026-10-20T09:30:00Z"}
    cand_ev = {"id": "777101", "ref": "refs/heads/b", "before": a,
               "after": "@@CAND@@", "created_at": "2026-10-20T14:00:00Z"}
    pt_back = ("## 2026-10-20T12:30:00Z | P-40 | evidence-back | backs=P-1-CONFIRM "
               "branch=b remote_sha=%s credential=OWNER-PAT-1 "
               "review_ref=docs/decisions/REGISTER.md#Part-98 evidence=\"a later "
               "push on the ref records the prior tip; citation: push_event=777100 "
               "ref=refs/heads/b prior_tip=%s; the event's after is elided\"\n" % (a, a))
    r = prior_tip_scenario(pt_back, [prior_ev, cand_ev])
    check("AT-31", "a prior_tip citation (later push's before == backed SHA) proves the push",
          r.returncode == 0 and "PUSH-PROVEN" in r.stdout
          and "prior_tip=" in r.stdout
          and "still observable live" in r.stdout, r.stdout)
    pt_wrong_ref = pt_back.replace("ref=refs/heads/b prior_tip", "ref=refs/heads/c prior_tip")
    r = prior_tip_scenario(pt_wrong_ref, [prior_ev, cand_ev])
    check("AT-31", "a prior_tip citation on the wrong ref fails",
          r.returncode == 1 and "mismatched ref-transition evidence" in r.stdout, r.stdout)
    pt_no_rr = pt_back.replace("review_ref=docs/decisions/REGISTER.md#Part-98 ", "")
    r = prior_tip_scenario(pt_no_rr, [prior_ev, cand_ev])
    check("AT-31", "a prior_tip citation without review_ref is not adjudicated",
          r.returncode == 1 and "no review_ref" in r.stdout, r.stdout)

    # 19-20. ADAPTED binding modes for this repository's feed (verified live
    # 2026-10-20): PushEvent.after is ELIDED in the events feed (only before
    # is reliable), and ref updates made through the git-data API generate NO
    # PushEvent at all. Mode 2: a PushEvent whose before == the declared
    # pre-push tip (after elided) + a push-event run on the candidate. Mode
    # 3: no PushEvent; the push-event run on the exact candidate alone.
    def adapted_binding_scenario(events, run_in_manifest=True):
        _seq[0] += 1
        intent = term.replace("P-30-INTENT", "P-3%d-INTENT" % _seq[0])
        text = base_log + back_ev_rec + intent
        c = commit(src, "logD", path=logrel, text=text)
        sp.run(["git", "-C", src, "push", "-q", "--force", bare, c + ":refs/heads/b"],
               check=True)
        ev_file = tmpfile(json.dumps(events), ".json")
        runs_map = {}
        if run_in_manifest:
            runs_map["555001"] = {"head_sha": c, "head_branch": "b",
                                  "event": "push", "conclusion": "success"}
        man = tmpfile(json.dumps({"runs": runs_map}), ".json")
        return run("evidence_sequence_check.py", logpath, "--pushlog-schema",
                   "--remote-verify", "--remote", bare, "--branch", "b",
                   "--expect-tip", c, "--runs-manifest", man,
                   "--push-events-json", ev_file)

    # 19. mode 2: the candidate's own push event carries before == declared
    #     tip with after ELIDED; the push-event run ties the candidate
    r = adapted_binding_scenario(
        [{"id": "777200", "ref": "refs/heads/b", "before": a, "after": None,
          "created_at": "2026-10-20T14:00:00Z"}])
    check("AT-31", "a before-matching PushEvent (after elided) + push-run binds the transition",
          r.returncode == 0 and "records the transition from the declared "
          "pre-push tip" in r.stdout and "feed elides after" in r.stdout,
          r.stdout)
    # 20. mode 3: no PushEvent at all (API-driven ref update class); the
    #     push-event run on the exact candidate is the binding evidence
    r = adapted_binding_scenario([])
    check("AT-31", "a push-event run on the exact candidate binds when no PushEvent exists",
          r.returncode == 0 and "push-event run 555001" in r.stdout
          and "ties the candidate to a genuine push event" in r.stdout,
          r.stdout)
    # 21. neither a PushEvent nor a push-event run -> fail closed
    r = adapted_binding_scenario([], run_in_manifest=False)
    check("AT-31", "no push-event evidence at all fails closed",
          r.returncode == 1 and "no push-event evidence exists" in r.stdout,
          r.stdout)


# --------------- AT-32 (EC-24 intent-reconciliation invariants, Round-22)

def test_at32_intent_reconciliation_invariants():
    """Round-22 direction 4: every intent-reconcile must target a PRECEDING
    eligible intent, match its branch semantics (push-intent: branch;
    merge-intent: source_branch), void it exactly once, and never void a
    consumed intent; the unmatched-intent invariant holds GLOBALLY (no branch
    may hold more than one unvoided unmatched intent, including branches
    other than the one under test)."""
    import subprocess as sp
    d = tempfile.mkdtemp()

    def schema(logtext):
        f = tmpfile(logtext, ".md")
        return run("evidence_sequence_check.py", f, "--pushlog-schema")

    def log(entry_ids):
        """Build a log with push-intents for branches b/c from ENTRY_IDS."""
        out = []
        for n_, (eid, br) in enumerate(entry_ids, 1):
            out.append("## 2026-10-20T1%d:00:00Z | %s | push-intent | branch=%s "
                       "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
                       % (n_, eid, br))
        return "".join(out)

    # 1. cross-branch void: intent for branch c, reconcile declares branch=b
    t = log([("P-1-INTENT", "c")]) + \
        "## 2026-10-20T19:00:00Z | P-9 | intent-reconcile | voids=P-1-INTENT " \
        "branch=b credential=OWNER-PAT-1 evidence=\"x\"\n"
    r = schema(t)
    check("AT-32", "a cross-branch void is rejected",
          r.returncode == 1 and "cross-branch voiding rejected" in r.stdout, r.stdout)
    # 2. voiding a confirmed (consumed) intent is rejected
    t = ("## 2026-10-20T10:00:00Z | P-1-INTENT | push-intent | branch=b "
         "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
         "## 2026-10-20T11:00:00Z | P-1-CONFIRM | push-confirm | branch=b remote_sha="
         + "a" * 40 + " credential=OWNER-PAT-1 evidence=\"ls-remote b -> " + "a" * 40 + "\"\n"
         "## 2026-10-20T12:00:00Z | P-9 | intent-reconcile | voids=P-1-INTENT "
         "branch=b credential=OWNER-PAT-1 evidence=\"x\"\n")
    r = schema(t)
    check("AT-32", "voiding an intent its CONFIRM already completed is rejected",
          r.returncode == 1 and "already completed by its CONFIRM" in r.stdout, r.stdout)
    # 3. duplicate void of the same intent is rejected
    t = log([("P-2-INTENT", "b"), ("P-3-INTENT", "b")]) + \
        "## 2026-10-20T19:00:00Z | P-9 | intent-reconcile | voids=P-2-INTENT branch=b " \
        "credential=OWNER-PAT-1 evidence=\"x\"\n" \
        "## 2026-10-20T19:30:00Z | P-10 | intent-reconcile | voids=P-2-INTENT branch=b " \
        "credential=OWNER-PAT-1 evidence=\"y\"\n"
    r = schema(t)
    check("AT-32", "a duplicate void of the same intent is rejected",
          r.returncode == 1 and "a second time" in r.stdout, r.stdout)
    # 4. a void that targets a FUTURE intent is rejected
    t = "## 2026-10-20T10:00:00Z | P-9 | intent-reconcile | voids=P-2-INTENT " \
        "branch=b credential=OWNER-PAT-1 evidence=\"x\"\n" + log([("P-2-INTENT", "b")])
    r = schema(t)
    check("AT-32", "a void of a non-preceding intent is rejected",
          r.returncode == 1 and "does not PRECEDE" in r.stdout, r.stdout)
    # 5. the unmatched-intent invariant is GLOBAL: two unmatched intents on
    #    branch c fail even though the branch under test is b
    bare = os.path.join(d, "remote.git")
    sp.run(["git", "init", "-q", "--bare", bare], check=True)
    src = os.path.join(d, "src")
    sp.run(["git", "init", "-q", "-b", "main", src], check=True)
    logrel = "docs/evidence/push-log/PUSH_LOG.md"
    logpath = os.path.join(src, logrel)

    def commit_log(text, tag):
        os.makedirs(os.path.dirname(os.path.join(src, logrel)), exist_ok=True)
        open(logpath, "w").write(text)
        sp.run(["git", "-C", src, "add", "docs"], check=True)
        sp.run(["git", "-C", src, "-c", "user.email=t@t", "-c", "user.name=t",
                "commit", "-q", "--allow-empty", "-m", tag], check=True)
        return sp.run(["git", "-C", src, "rev-parse", "HEAD"],
                      capture_output=True, text=True).stdout.strip()

    _seq = [0]

    def remote_check(text, *args):
        # unique terminal-intent id per call (the introducing-commit binding
        # uses git log -S; a reused id would resolve to an earlier candidate)
        _seq[0] += 1
        text = text.replace("P-4-INTENT", "P-4%d-INTENT" % _seq[0])
        c = commit_log(text, "log")
        sp.run(["git", "-C", src, "push", "-q", bare, c + ":refs/heads/b"], check=True)
        return run("evidence_sequence_check.py", logpath, "--pushlog-schema",
                   "--remote-verify", "--remote", bare, "--branch", "b",
                   "--expect-tip", c, *args)

    two_c = log([("P-4-INTENT", "b"), ("P-5-INTENT", "c"), ("P-6-INTENT", "c")])
    r = remote_check(two_c)
    check("AT-32", "multiple unmatched intents on a branch OTHER than the one under test fail",
          r.returncode == 1 and "multiple unmatched intents for branch c" in r.stdout,
          r.stdout)
    # 6. the excess intent on the other branch, validly voided, passes
    fixed = two_c + \
        "## 2026-10-20T19:00:00Z | P-11 | intent-reconcile | voids=P-5-INTENT " \
        "branch=c credential=OWNER-PAT-1 evidence=\"superseded by P-6\"\n"
    r = remote_check(fixed)
    check("AT-32", "an excess intent on another branch, validly voided, is tolerated",
          r.returncode == 0 and "PENDING" in r.stdout, r.stdout)
    # 7. merge-intent voids are scoped to the SOURCE branch
    mi = ("## 2026-10-20T10:00:00Z | P-4-INTENT | push-intent | branch=feat "
          "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
          "## 2026-10-20T11:00:00Z | P-5 | merge-intent | target=main pr=7 "
          "source_branch=feat expected_pre_merge_base=" + "d" * 40 + " "
          "expected_head=DEFERRED credential=OWNER-PAT-1 authorized_by=\"t\"\n")
    ok_void = mi + ("## 2026-10-20T12:00:00Z | P-12 | intent-reconcile | "
                    "voids=P-5 branch=feat credential=OWNER-PAT-1 evidence=\"stale\"\n"
                    "## 2026-10-20T13:00:00Z | P-6-INTENT | push-intent | branch=b "
                    "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n")
    r = schema(ok_void)
    check("AT-32", "a merge-intent void scoped to its source branch is valid",
          r.returncode == 0, r.stdout)
    bad_scope = mi + ("## 2026-10-20T12:00:00Z | P-12 | intent-reconcile | "
                      "voids=P-5 branch=main credential=OWNER-PAT-1 evidence=\"stale\"\n")
    r = schema(bad_scope)
    check("AT-32", "a merge-intent void scoped to the target branch is rejected",
          r.returncode == 1 and "cross-branch voiding rejected" in r.stdout, r.stdout)
    # 8. a CONFIRM naming a VOIDED intent is rejected (a void must not
    #    substitute for a confirmation)
    t = ("## 2026-10-20T10:00:00Z | P-7-INTENT | push-intent | branch=b "
         "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
         "## 2026-10-20T11:00:00Z | P-13 | intent-reconcile | voids=P-7-INTENT "
         "branch=b credential=OWNER-PAT-1 evidence=\"stale\"\n"
         "## 2026-10-20T12:00:00Z | P-7-CONFIRM | push-confirm | branch=b remote_sha="
         + "a" * 40 + " credential=OWNER-PAT-1 evidence=\"ls-remote b -> " + "a" * 40 + "\"\n")
    r = schema(t)
    check("AT-32", "a confirmation pairing a voided intent is rejected",
          r.returncode == 1 and "no longer eligible" in r.stdout, r.stdout)
    # 9. a void naming a non-intent entry is rejected
    t = log([("P-8-INTENT", "b")]) + \
        "## 2026-10-20T19:00:00Z | P-14 | intent-reconcile | voids=P-8-CONFIRM " \
        "branch=b credential=OWNER-PAT-1 evidence=\"x\"\n"
    r = schema(t)
    check("AT-32", "a void naming a non-intent entry is rejected (dangling void)",
          r.returncode == 1 and "not an intent entry" in r.stdout, r.stdout)


# ------------- AT-33 (merge-transition proof + pre-merge gate, Round-22)

def test_at33_merge_proof_and_premerge_gate():
    """Round-22 direction 5: the merge transition is proven by the PR API's
    merge_commit_sha (not parents alone), the exact source ref, merged_at,
    and a FRESH non-author approval on the exact head (commit_id, submitted
    before the merge, latest state per reviewer); the merge PushEvent's
    before must equal the declared base; a second sequential merge follows
    the first; and the PRE-MERGE gate fails while no approval exists on the
    exact head."""
    import subprocess as sp
    d = tempfile.mkdtemp()
    bare = os.path.join(d, "remote.git")
    sp.run(["git", "init", "-q", "--bare", bare], check=True)
    src = os.path.join(d, "src")
    sp.run(["git", "init", "-q", "-b", "main", src], check=True)

    def commit(tag, path=None, text=None):
        if path and text is not None:
            os.makedirs(os.path.dirname(os.path.join(src, path)), exist_ok=True)
            open(os.path.join(src, path), "w").write(text)
            sp.run(["git", "-C", src, "add", "docs"], check=True)
        sp.run(["git", "-C", src, "-c", "user.email=t@t", "-c", "user.name=t",
                "commit", "-q", "--allow-empty", "-m", tag], check=True)
        return sp.run(["git", "-C", src, "rev-parse", "HEAD"],
                      capture_output=True, text=True).stdout.strip()

    def push(refspec, force=False):
        sp.run(["git", "-C", src, "push", "-q"] + (["--force"] if force else []) +
               [bare, refspec], check=True)

    logrel = "docs/evidence/push-log/PUSH_LOG.md"
    logpath = os.path.join(src, logrel)
    base = commit("base")
    push("HEAD:refs/heads/main")
    sp.run(["git", "-C", src, "checkout", "-q", "-b", "feat"], check=True)
    mi = ("## 2026-10-20T10:00:00Z | P-1-INTENT | push-intent | branch=feat "
          "expected_remote_tip_before_push=ABSENT credential=OWNER-PAT-1\n"
          "## 2026-10-20T11:00:00Z | P-2 | merge-intent | target=main pr=7 "
          "source_branch=feat expected_pre_merge_base=%s expected_head=DEFERRED "
          "credential=OWNER-PAT-1 authorized_by=\"auditor direction\"\n" % base)
    head2 = commit("merge-intent", path=logrel, text=mi)
    push("HEAD:refs/heads/feat")

    def chk(branch, tip, pr_file=None, rev_file=None, ev_file=None, pre=False,
            pr_number=None):
        extra = []
        if pr_file:
            extra += ["--pr-json", pr_file]
        if rev_file:
            extra += ["--reviews-json", rev_file]
        if ev_file:
            extra += ["--push-events-json", ev_file]
        if pre:
            extra += ["--pre-merge-verify", "--pr-number", str(pr_number or 7)]
        return run("evidence_sequence_check.py", logpath, "--pushlog-schema",
                   "--remote-verify", "--remote", bare, "--branch", branch,
                   "--expect-tip", tip, *extra)

    # ---- phase 1: PRE-MERGE gate on the source branch
    pr_open = tmpfile(json.dumps({"merged": False, "merge_commit_sha": None,
                                  "merged_at": None, "user": {"login": "author"},
                                  "head": {"sha": head2, "ref": "feat"},
                                  "base": {"ref": "main"}}), ".json")
    rev_none = tmpfile(json.dumps([]), ".json")
    rev_ok = tmpfile(json.dumps([{"state": "APPROVED", "user": {"login": "reviewer"},
                                  "commit_id": head2,
                                  "submitted_at": "2026-10-20T11:30:00Z"}]), ".json")
    r = chk("feat", head2, pr_open, rev_none, pre=True)
    check("AT-33", "the pre-merge gate FAILS while no approval exists on the exact head",
          r.returncode == 1 and "PRE-MERGE BLOCKED" in r.stdout, r.stdout)
    r = chk("feat", head2, pr_open, rev_ok, pre=True)
    check("AT-33", "the pre-merge gate passes with a fresh non-author approval on the exact head",
          r.returncode == 0 and "PRE-MERGE-VERIFIED" in r.stdout, r.stdout)
    pr_other_head = tmpfile(json.dumps({"merged": False, "merge_commit_sha": None,
                                        "merged_at": None, "user": {"login": "author"},
                                        "head": {"sha": "9" * 40, "ref": "feat"},
                                        "base": {"ref": "main"}}), ".json")
    r = chk("feat", head2, pr_other_head, rev_ok, pre=True)
    check("AT-33", "the pre-merge gate rejects a PR whose head is not the candidate",
          r.returncode == 1 and "not the candidate" in r.stdout, r.stdout)
    pr_wrong_base = tmpfile(json.dumps({"merged": False, "merge_commit_sha": None,
                                        "merged_at": None, "user": {"login": "author"},
                                        "head": {"sha": head2, "ref": "feat"},
                                        "base": {"ref": "dev"}}), ".json")
    r = chk("feat", head2, pr_wrong_base, rev_ok, pre=True)
    check("AT-33", "the pre-merge gate rejects a PR based on the wrong target",
          r.returncode == 1 and "not the merge-intent's declared target" in r.stdout,
          r.stdout)
    r = chk("feat", head2, pre=True)  # no fixtures, no API
    check("AT-33", "the pre-merge gate fails closed when the PR cannot be verified",
          r.returncode == 1 and "could not be verified pre-merge" in r.stdout, r.stdout)

    # ---- phase 2: DEFERRED supersession + a fresh terminal intent. The
    # candidate advances on the SAME source branch (the real-world
    # P-061 -> P-068 -> P-072 pattern); the prior merge-intent is superseded
    # via the introducing commits, and the two-phase protocol delivers
    # P-1-CONFIRM plus a fresh terminal intent for the new candidate.
    sp.run(["git", "-C", src, "checkout", "-q", "feat"], check=True)
    mi_chain = (mi
                + "## 2026-10-20T12:30:00Z | P-1-CONFIRM | push-confirm | branch=feat "
                  "remote_sha=%s credential=OWNER-PAT-1 evidence=\"git ls-remote r "
                  "refs/heads/feat -> %s\"\n" % (head2, head2)
                + "## 2026-10-20T13:00:00Z | P-4 | merge-intent | target=main pr=7 "
                  "source_branch=feat expected_pre_merge_base=%s expected_head=DEFERRED "
                  "credential=OWNER-PAT-1 authorized_by=\"auditor direction\"\n" % base
                + "## 2026-10-20T13:30:00Z | P-5-INTENT | push-intent | branch=feat "
                  "expected_remote_tip_before_push=%s credential=OWNER-PAT-1\n" % head2)
    head6 = commit("redeclared merge-intent", path=logrel, text=mi_chain)
    push("HEAD:refs/heads/feat")
    pr_red = tmpfile(json.dumps({"merged": False, "merge_commit_sha": None,
                                 "merged_at": None, "user": {"login": "author"},
                                 "head": {"sha": head6, "ref": "feat"},
                                 "base": {"ref": "main"}}), ".json")
    rev_red = tmpfile(json.dumps([{"state": "APPROVED",
                                   "user": {"login": "reviewer"},
                                   "commit_id": head6,
                                   "submitted_at": "2026-10-20T13:45:00Z"}]), ".json")
    r = chk("feat", head6, pr_red, rev_red, pre=True)
    check("AT-33", "a prior DEFERRED merge-intent on the same lineage is SUPERSEDED, not a failure",
          r.returncode == 0 and "SUPERSEDED" in r.stdout
          and "PRE-MERGE-VERIFIED" in r.stdout, r.stdout)

    # ---- phase 3: the reviewed merge lands; post-merge proof hardening
    # (merging the CURRENT feat head, which the active P-4 intent declares)
    sp.run(["git", "-C", src, "checkout", "-q", "main"], check=True)
    sp.run(["git", "-C", src, "-c", "user.email=t@t", "-c", "user.name=t",
            "merge", "-q", "--no-ff", "-m", "merge pr 7", "feat"], check=True)
    m = sp.run(["git", "-C", src, "rev-parse", "HEAD"],
               capture_output=True, text=True).stdout.strip()
    push("HEAD:refs/heads/main")

    def prf(**over):
        doc = {"merged": True, "merge_commit_sha": m,
               "merged_at": "2026-10-20T14:00:00Z", "user": {"login": "author"},
               "head": {"sha": head6, "ref": "feat"}, "base": {"ref": "main"}}
        doc.update(over)
        return tmpfile(json.dumps(doc), ".json")

    ev_merge = tmpfile(json.dumps([{"id": "888001", "ref": "refs/heads/main",
                                    "before": base, "after": m,
                                    "created_at": "2026-10-20T14:00:05Z"}]), ".json")
    r = chk("main", m, prf(), rev_red, ev_merge)
    check("AT-33", "a legitimate merge with complete PR + approval + PushEvent proof passes",
          r.returncode == 0 and "MERGE-VERIFIED" in r.stdout
          and "merge PushEvent verifies" in r.stdout, r.stdout)
    # mismatched merge_commit_sha fails (parents alone are insufficient)
    r = chk("main", m, prf(merge_commit_sha="5" * 40), rev_ok, ev_merge)
    check("AT-33", "a mismatched merge_commit_sha fails",
          r.returncode == 1 and "merge_commit_sha is" in r.stdout, r.stdout)
    # wrong source branch fails
    r = chk("main", m, prf(head={"sha": head2, "ref": "other"}), rev_ok, ev_merge)
    check("AT-33", "a PR merged from the wrong source branch fails",
          r.returncode == 1 and "not the declared source branch" in r.stdout, r.stdout)
    # stale-SHA approval (approval bound to an older head) fails
    rev_stale = tmpfile(json.dumps([{"state": "APPROVED",
                                     "user": {"login": "reviewer"},
                                     "commit_id": base,
                                     "submitted_at": "2026-10-20T11:30:00Z"}]), ".json")
    r = chk("main", m, prf(), rev_stale, ev_merge)
    check("AT-33", "an approval bound to a stale commit does not satisfy the exact-head rule",
          r.returncode == 1 and "approval bound to commit" in r.stdout, r.stdout)
    # approval submitted AFTER the merge fails
    rev_late = tmpfile(json.dumps([{"state": "APPROVED", "user": {"login": "reviewer"},
                                    "commit_id": head6,
                                    "submitted_at": "2026-10-20T14:30:00Z"}]), ".json")
    r = chk("main", m, prf(), rev_late, ev_merge)
    check("AT-33", "an approval submitted after the merge fails",
          r.returncode == 1 and "submitted after the merge" in r.stdout, r.stdout)
    # dismissed / superseded approvals fail
    rev_dismissed = tmpfile(json.dumps([
        {"state": "APPROVED", "user": {"login": "reviewer"}, "commit_id": head6,
         "submitted_at": "2026-10-20T13:45:00Z"},
        {"state": "DISMISSED", "user": {"login": "reviewer"}, "commit_id": head6,
         "submitted_at": "2026-10-20T13:50:00Z"}]), ".json")
    r = chk("main", m, prf(), rev_dismissed, ev_merge)
    check("AT-33", "a later DISMISSED state supersedes the earlier approval",
          r.returncode == 1 and "latest state 'DISMISSED'" in r.stdout, r.stdout)
    rev_superseded = tmpfile(json.dumps([
        {"state": "APPROVED", "user": {"login": "reviewer"}, "commit_id": head6,
         "submitted_at": "2026-10-20T13:45:00Z"},
        {"state": "CHANGES_REQUESTED", "user": {"login": "reviewer"},
         "commit_id": head6, "submitted_at": "2026-10-20T13:50:00Z"}]), ".json")
    r = chk("main", m, prf(), rev_superseded, ev_merge)
    check("AT-33", "a later CHANGES_REQUESTED state supersedes the earlier approval",
          r.returncode == 1 and "latest state 'CHANGES_REQUESTED'" in r.stdout, r.stdout)
    # raced merge PushEvent (before != declared base) fails
    ev_raced = tmpfile(json.dumps([{"id": "888002", "ref": "refs/heads/main",
                                    "before": "c" * 40, "after": m,
                                    "created_at": "2026-10-20T12:00:05Z"}]), ".json")
    r = chk("main", m, prf(), rev_ok, ev_raced)
    check("AT-33", "a merge PushEvent whose before is not the declared base fails",
          r.returncode == 1
          and "the protected-branch transition is not the declared one" in r.stdout,
          r.stdout)

    # ---- phase 3: sequential second merge (the PR #30 -> #31 pattern)
    sp.run(["git", "-C", src, "checkout", "-q", "-b", "feat2", "main"], check=True)
    mi2 = mi + ("## 2026-10-20T13:00:00Z | P-3 | merge-intent | target=main pr=8 "
                "source_branch=feat2 expected_pre_merge_base=%s expected_head=DEFERRED "
                "credential=OWNER-PAT-1 authorized_by=\"auditor direction\"\n" % m)
    head5 = commit("second merge-intent", path=logrel, text=mi2)
    push("HEAD:refs/heads/feat2")
    sp.run(["git", "-C", src, "checkout", "-q", "main"], check=True)
    sp.run(["git", "-C", src, "-c", "user.email=t@t", "-c", "user.name=t",
            "merge", "-q", "--no-ff", "-m", "merge pr 8", "feat2"], check=True)
    m2 = sp.run(["git", "-C", src, "rev-parse", "HEAD"],
                capture_output=True, text=True).stdout.strip()
    push("HEAD:refs/heads/main")
    pr8 = tmpfile(json.dumps({"merged": True, "merge_commit_sha": m2,
                              "merged_at": "2026-10-20T14:00:00Z",
                              "user": {"login": "author"},
                              "head": {"sha": head5, "ref": "feat2"},
                              "base": {"ref": "main"}}), ".json")
    rev8 = tmpfile(json.dumps([{"state": "APPROVED", "user": {"login": "reviewer2"},
                                "commit_id": head5,
                                "submitted_at": "2026-10-20T13:30:00Z"}]), ".json")
    ev8 = tmpfile(json.dumps([{"id": "888003", "ref": "refs/heads/main",
                               "before": m, "after": m2,
                               "created_at": "2026-10-20T14:00:05Z"}]), ".json")
    r = chk("main", m2, pr8, rev8, ev8)
    check("AT-33", "a second sequential merge on the new base verifies (the #30 -> #31 pattern)",
          r.returncode == 0 and "MERGE-VERIFIED" in r.stdout
          and "P-3" in r.stdout, r.stdout)

    # ---- phase 4 (closed): the sequential second merge completed the scenario;
    # the DEFERRED-supersession proof ran in phase 2 above.


# ------------------- AT-34 (EC-26 conservative glob coverage, Round-22)

def test_at34_glob_coverage_and_ordered_filters():
    """Round-22 direction 6: representative-sample coverage is replaced by a
    sound conservative glob-coverage check; ordered negative patterns follow
    GitHub Actions filter semantics; empty filter lists and malformed event
    configurations fail closed."""
    d = tempfile.mkdtemp()
    os.makedirs(os.path.join(d, ".github", "workflows"), exist_ok=True)
    wf = os.path.join(d, ".github", "workflows", "gate.yml")
    os.makedirs(os.path.join(d, "docs", "directives"), exist_ok=True)
    open(os.path.join(d, "registry.json"), "w").write(json.dumps({"tests": {}}))

    def write_wf(text):
        open(wf, "w").write(text)

    def idx(enf):
        open(os.path.join(d, "docs", "directives", "INDEX.md"), "w").write(
            "# index\n<!-- BEGIN RULE INVENTORY -->\n"
            "reading-note-path: docs/note.md\n"
            "pinned-blobs: b6bb35e84956d84a8c9d9e65734c212baf3ed075\n"
            "| Rule ID | Status | Enforcement point(s) | Acceptance test(s) |\n"
            "|---|---|---|---|\n"
            "| EC-X | ACTIVE | " + enf + " | NONE |\n"
            "<!-- END RULE INVENTORY -->\n")
        open(os.path.join(d, "docs", "note.md"), "w").write(
            "blob b6bb35e84956d84a8c9d9e65734c212baf3ed075 read\n| EC-X | -> mapped |\n")

    def run_idx():
        return run("reading_note_map_check.py", "--index", "docs/directives/INDEX.md",
                   "--repo-root", d, "--registry", "registry.json")

    job = "jobs:\n  gate-job:\n    name: Gate Job\n    runs-on: ubuntu-latest\n"

    def on_with(push_cfg, pr_cfg="    branches: ['**']\n"):
        return ("name: gate\non:\n  push:\n" + push_cfg +
                "  pull_request:\n" + pr_cfg + job)

    gov = "~events=push+pull_request~governs=docs/**"
    # 1. docs/** governed by only docs/X: one representative sample would pass,
    #    the conservative witness set fails
    write_wf(on_with("    branches: ['**']\n    paths: [docs/X]\n"))
    idx("workflow:.github/workflows/gate.yml#gate-job" + gov)
    r = run_idx()
    check("AT-34", "docs/** governed by only docs/X fails (representative-sample hole closed)",
          r.returncode == 1 and "would exclude the governed path" in r.stdout, r.stdout)
    # 2. ordered negative pattern excludes the governed subtree
    write_wf(on_with("    branches: ['**']\n    paths: ['**', '!docs/evidence/**']\n"))
    idx("workflow:.github/workflows/gate.yml#gate-job"
        "~events=push+pull_request~governs=docs/evidence/**")
    r = run_idx()
    check("AT-34", "an ordered negative pattern excluding the governed subtree fails",
          r.returncode == 1 and "ordered negative patterns" in r.stdout, r.stdout)
    # 3. the SAME patterns in the OPPOSITE order re-include the subtree (order matters)
    write_wf(on_with("    branches: ['**']\n    paths: ['!docs/evidence/**', '**']\n"))
    r = run_idx()
    check("AT-34", "reversing the pattern order changes the outcome (ordered semantics)",
          r.returncode == 0, r.stdout)
    # 4. branch filters use ordered semantics too: ['**', '!main'] excludes main
    write_wf(on_with("    branches: ['**', '!main']\n"))
    idx("workflow:.github/workflows/gate.yml#gate-job"
        "~events=push+pull_request~governs=docs/**")
    r = run_idx()
    check("AT-34", "branches ['**', '!main'] excludes the enforcement branch",
          r.returncode == 1 and "ordered filter semantics" in r.stdout, r.stdout)
    # 5. empty filter lists fail closed (never read as unrestricted)
    write_wf(on_with("    branches: []\n"))
    r = run_idx()
    check("AT-34", "an empty branches list fails closed",
          r.returncode == 1 and "EMPTY branches list" in r.stdout, r.stdout)
    write_wf(on_with("    branches: ['**']\n    paths: []\n"))
    r = run_idx()
    check("AT-34", "an empty paths list fails closed",
          r.returncode == 1 and "EMPTY paths list" in r.stdout, r.stdout)
    # 6. malformed event configuration types fail closed
    write_wf("name: gate\non:\n  push:\n  - main\n  - dev\n  pull_request:\n"
             "    branches: ['**']\n" + job)
    r = run_idx()
    check("AT-34", "a list-valued event configuration fails closed",
          r.returncode == 1 and "MALFORMED event configuration" in r.stdout, r.stdout)
    write_wf("name: gate\non:\n  push: main\n  pull_request:\n"
             "    branches: ['**']\n" + job)
    r = run_idx()
    check("AT-34", "a string-valued event configuration fails closed",
          r.returncode == 1 and "MALFORMED event configuration" in r.stdout, r.stdout)
    write_wf("name: gate\non:\n  push: 5\n  pull_request:\n"
             "    branches: ['**']\n" + job)
    r = run_idx()
    check("AT-34", "an int-valued event configuration fails closed",
          r.returncode == 1 and "MALFORMED event configuration" in r.stdout, r.stdout)
    # 7. broad coverage passes: every witness of both globs is covered
    write_wf(on_with("    branches: ['**']\n    paths: [docs/**, scripts/**]\n"))
    idx("workflow:.github/workflows/gate.yml#gate-job"
        "~events=push+pull_request~governs=docs/**+scripts/ci/**")
    r = run_idx()
    check("AT-34", "a filter covering every structural witness of the globs passes",
          r.returncode == 0, r.stdout)
    # 8. a component-level wildcard still covers its own glob
    write_wf(on_with("    branches: ['**']\n    paths: [docs/**]\n"))
    idx("workflow:.github/workflows/gate.yml#gate-job"
        "~events=push+pull_request~governs=docs/*.md")
    r = run_idx()
    check("AT-34", "docs/*.md is covered by docs/** (witness shapes all match)",
          r.returncode == 0, r.stdout)
    # 9. paths-ignore over a witness SHAPE fails (the old representative
    #    sample docs/X would NOT have matched this ignore -- the witness set
    #    closes that hole)
    write_wf(on_with("    branches: ['**']\n    paths-ignore: [docs/*/w1]\n"))
    idx("workflow:.github/workflows/gate.yml#gate-job"
        "~events=push+pull_request~governs=docs/**")
    r = run_idx()
    check("AT-34", "a paths-ignore covering a governed witness shape fails",
          r.returncode == 1 and "paths-ignore" in r.stdout
          and "excludes the governed path" in r.stdout, r.stdout)
