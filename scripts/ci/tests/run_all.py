#!/usr/bin/env python3
"""Acceptance tests AT-01..AT-07 for the EC-16..EC-26 enforcement scaffolding.

EC-20 dogfood: this package names its acceptance tests and they must exist and pass.
Plain asserts, stdlib only; runs each checker as a subprocess. Exit 0 = all pass.
"""
import os
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
                          capture_output=True, text=True, cwd=cwd)


def check(at_id, name, cond, detail=""):
    (PASS if cond else FAIL).append(f"{at_id} {name}" + (f" -- {detail}" if not cond else ""))


def tmpfile(content: str) -> str:
    f = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
    f.write(content)
    f.close()
    return f.name


def at01_claimed_shas():
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                          cwd=REPO).stdout.strip()
    good = tmpfile(f"Claimed-SHAs:\n  {head}\nDirective-ID: TEST\nCloses-items: TEST\nStatus-summary: PUSHED\n")
    bad = tmpfile("Claimed-SHAs:\n  " + "0" * 40 + "\n")
    r = run("verify_claimed_shas.py", "--body", good, "--local")
    check("AT-01", "real HEAD passes offline verification", r.returncode == 0, r.stdout)
    r = run("verify_claimed_shas.py", "--body", bad, "--local")
    check("AT-01", "absent SHA fails", r.returncode == 1 and "not on remote" in r.stdout, r.stdout)
    none = tmpfile("no claimed shas here\n")
    r = run("verify_claimed_shas.py", "--body", none, "--local")
    check("AT-01", "missing field fails closed", r.returncode == 1 and "no full 40-hex" in r.stdout, r.stdout)
    r = run("verify_claimed_shas.py", "--body", bad, "--require-directive-id-only")
    check("AT-01", "directive-only mode ignores SHA absence, flags missing Directive-ID",
          r.returncode == 1 and "Directive-ID" in r.stdout and "not on remote" not in r.stdout, r.stdout)


def at02_declared_base():
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                          cwd=REPO).stdout.strip()
    base = subprocess.run(["git", "rev-parse", "HEAD~1"], capture_output=True, text=True,
                          cwd=REPO).stdout.strip()
    r = run("check_declared_base.py", "--base", base, "--head", head)
    check("AT-02", "ancestor base passes", r.returncode == 0, r.stdout)
    r = run("check_declared_base.py", "--base", "0" * 40, "--head", head)
    check("AT-02", "absent base fails closed", r.returncode == 1 and "not a commit" in r.stdout, r.stdout)


def at03_deploy_lint():
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
  push:
    branches: ["**"]
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
    r = run("workflow_deploy_trigger_lint.py", d)
    check("AT-03", "deploy-on-push flagged", r.returncode == 1 and "w0.yml" in r.stdout, r.stdout)
    check("AT-03", "deploy-on-PR and test jobs not flagged", "w1.yml" not in r.stdout
          and "w2.yml" not in r.stdout, r.stdout)
    r = run("workflow_deploy_trigger_lint.py", d, "--report-only")
    check("AT-03", "--report-only exits 0 with findings listed",
          r.returncode == 0 and "w0.yml" in r.stdout, r.stdout)


def at04_status_vocabulary():
    bad = tmpfile("Suite result: RUN_PASS\nProgress: 92%\nGitleaks came back clean scan\n"
                  "Failures are pre-existing\n")
    good = tmpfile("Suite result: RUN_PASS(run 118475032870)\nProgress: 92% (run 118475032870)\n"
                   "Failures are pre-existing (base 2195a04)\n")
    r = run("status_vocabulary_lint.py", bad)
    check("AT-04", "RUN_PASS without id flagged", "RUN_PASS without run id" in r.stdout, r.stdout)
    check("AT-04", "percentage without id flagged", "percentage claim without a run id" in r.stdout, r.stdout)
    check("AT-04", "clean scan without citation flagged", "'clean scan'" in r.stdout, r.stdout)
    check("AT-04", "pre-existing without citation flagged", "'pre-existing'" in r.stdout, r.stdout)
    r = run("status_vocabulary_lint.py", good)
    check("AT-04", "run-id-cited lines pass", r.returncode == 0, r.stdout)


def at05_sequence():
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


def at06_reading_note():
    note = tmpfile("blob b6bb35e8 read in full\nEC-16 -> proposals/ci/verify-claimed-shas.yml\n")
    r = run("reading_note_map_check.py", "--note", note, "--blobs", "b6bb35e8,211ad2c5",
            "--rules", "EC-16,EC-17")
    check("AT-06", "missing blob and unmapped rule flagged",
          r.returncode == 1 and "211ad2c5" in r.stdout and "EC-17" in r.stdout, r.stdout)
    good = tmpfile("blob b6bb35e8 ... blob 211ad2c5\n| EC-16 | yml |\n| EC-17 | practice |\n")
    r = run("reading_note_map_check.py", "--note", good, "--blobs", "b6bb35e8,211ad2c5",
            "--rules", "EC-16,EC-17")
    check("AT-06", "complete note passes", r.returncode == 0, r.stdout)


def at07_pr_template():
    tpl = open(os.path.join(REPO, "proposals", "ci", "PULL_REQUEST_TEMPLATE.md"),
               encoding="utf-8").read()
    for field in ("Directive-ID", "Claimed-SHAs", "Base-SHA", "Acceptance-Tests", "Status-summary"):
        check("AT-07", f"PR template carries {field}", field + ":" in tpl)


def main() -> int:
    at01_claimed_shas()
    at02_declared_base()
    at03_deploy_lint()
    at04_status_vocabulary()
    at05_sequence()
    at06_reading_note()
    at07_pr_template()
    print(f"PASS {len(PASS)}: " + "; ".join(PASS))
    if FAIL:
        print(f"FAIL {len(FAIL)}: " + "; ".join(FAIL))
        return 1
    print("ALL ACCEPTANCE TESTS PASS (AT-01..AT-07)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
