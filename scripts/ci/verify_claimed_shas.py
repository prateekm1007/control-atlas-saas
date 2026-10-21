#!/usr/bin/env python3
"""EC-16/19/20 claim verification (inert until adopted; called by proposals/ci/verify-claimed-shas.yml).

Auditor revisions 2026-10-20 implemented here:
  EC-16 (item 3.5):
    * STRICT parsing: the Claimed-SHAs field is parsed entry by entry; an entry that is not a
      full 40-character lowercase hex commit ID rejects the WHOLE list (mixed valid/invalid lists
      are never partially accepted; short or malformed IDs are never silently ignored).
    * Remote verification uses GET /repos/{repo}/commits/{sha}; the checker requires the returned
      commit SHA to EQUAL the requested SHA exactly.
    * 404/422 -> the commit is absent (a real verification result). Timeouts, permission errors
      (401/403), rate limiting (429), 5xx, and network failures are REMOTE VERIFICATION FAILURES
      (fail closed) and are never reported as "commit absent" or as success.
    * Least-privilege + redaction: the token is read from GH_TOKEN/GITHUB_TOKEN only, is never
      printed, and error output never includes request headers or credentials.
  EC-20 (item 3.2):
    * The former --check-tests-in-diff presence check is REMOVED (a path-substring match cannot
      validate pytest node IDs and does not execute anything).
    * --check-tests instead: parses Acceptance-Tests: from the PR body, verifies each identifier
      is collectable by pytest on the exact candidate (missing identifier -> fail), then EXECUTES
      the named tests via pytest. No user-controlled identifier is ever placed in an interpolated
      shell command: identifiers are validated against the collected set and passed as argument
      ARRAY elements (subprocess without shell). The command argv, candidate SHA, environment,
      run IDs, counts, and exit status are captured into an evidence JSON (--evidence-out).

Round-24 hardening (auditor direction 7 -- the gate passed while the PR body
  still described a DIFFERENT candidate):
    * --expect-head-sha: the EXACT current pull-request head SHA, passed as an
      explicit input by the workflow (never inferred from the body's own claim
      list). The body's Claimed-SHAs MUST contain it. A body that declares only
      stale -- but real -- SHAs FAILS.
    * The body's Status-summary MUST cite a run that (a) is real, (b) succeeded
      and (c) ran on that SAME exact head SHA. A real historical run id that
      describes another commit is REJECTED: "real but stale" is not evidence for
      the current candidate.
    * The cited run is verified through the Actions API (or --runs-manifest for
      offline tests); an unverifiable citation FAILS (fail closed).

Other modes:
  --local                     offline claimed-SHA verification against the local object store (tests).
  --require-directive-id-only EC-19 check only (Directive-ID / Closes-items / Status-summary present).

Exit codes: 0 = pass; 1 = fail (findings on stdout, one per line, `FAIL:`/`VERIFICATION_FAILURE:` prefixed);
2 = usage/internal error (also fail closed).
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

FULL_SHA = re.compile(r"^[0-9a-f]{40}$")
# The required checks of this repository. A Status-summary claim is only
# evidence when it cites a successful run of one of THESE workflows on the
# exact candidate head; a successful unrelated workflow on the same commit is
# not required-check evidence (Round-25 direction 3).
DEFAULT_REQUIRED_WORKFLOWS = ("Toscanini Stability Gate", "Governance Acceptance")
FIELDS = ("Directive-ID", "Closes-items", "Base-SHA", "Stack-parent",
          "Claimed-SHAs", "Acceptance-Tests", "Waiver", "Status-summary")


def field_values(body: str, field: str) -> list[str]:
    """Values of a template field: remainder of its line plus following indented lines."""
    out: list[str] = []
    pattern = re.compile(rf"^{re.escape(field)}:(.*)$")
    lines = body.splitlines()
    for i, line in enumerate(lines):
        m = pattern.match(line)
        if m:
            first = m.group(1).strip()
            if first and not first.startswith("<!--"):
                out.append(first)
            for j in range(i + 1, len(lines)):
                nxt = lines[j]
                if nxt.startswith(("  ", "\t")) and nxt.strip() and not nxt.strip().startswith("<!--"):
                    out.append(nxt.strip())
                else:
                    break
    return out


# ---------------------------------------------------------------- EC-16 ----

STATUS_LINE = re.compile(r"^\s*Status-summary:\s*(.+?)\s*$")
RUNPASS = re.compile(r"RUN_PASS\s*\(\s*run_id=(\d+)\s*,\s*sha=([0-9a-f]{40})\s*\)")


def api_run(repo: str, rid: str, token: str | None):
    """Verify one Actions run. Returns (status, detail) with status in
    {"present", "absent", "error"}; "present" carries the run record."""
    url = f"https://api.github.com/repos/{repo}/actions/runs/{rid}"
    req = urllib.request.Request(url, method="GET",
                                 headers={"Accept": "application/vnd.github+json",
                                          "X-GitHub-Api-Version": "2022-11-28"})
    if token:
        req.add_header("Authorization", f"token {token}")
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            body = json.loads(resp.read().decode("utf-8"))
            if body.get("id") is None:
                return ("error", "the Actions API returned no run record")
            return ("present", body)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return ("absent", f"HTTP 404: run {rid} not found on {repo}")
        return ("error", f"HTTP {e.code}: remote verification failure (fail closed)")
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        return ("error", f"network failure during remote verification (fail closed): "
                         f"{type(e).__name__}")
    except ValueError:
        return ("error", "unparseable API response (fail closed)")


def parse_claimed_shas(entries: list[str]):
    """Strict per-entry parse. Returns (shas, errors).

    ANY malformed/short/mixed-case entry rejects the whole list (auditor revision:
    'Reject mixed valid/invalid lists; never silently ignore short or malformed IDs').
    """
    shas, errors = [], []
    raw_tokens = [t for v in entries for t in v.replace(",", " ").split()]
    if not raw_tokens:
        return [], ["Claimed-SHAs field is empty (fail closed)"]
    for tok in raw_tokens:
        if not FULL_SHA.match(tok):
            errors.append(
                f"malformed Claimed-SHAs entry {tok!r}: entries must be full 40-character "
                f"lowercase hex commit IDs; the whole list is rejected")
            continue
        shas.append(tok)
    if errors:
        return [], errors  # reject the whole list on any malformed entry
    return shas, []


def verify_sha_local(sha: str):
    return ("present", sha) if subprocess.run(
        ["git", "cat-file", "-e", f"{sha}^{{commit}}"], capture_output=True).returncode == 0 \
        else ("absent", None)


def verify_sha_api(repo: str, sha: str, token: str | None, fixture: dict | None = None):
    """Verify one SHA against the intended repository via the GitHub commits API.

    Returns (status, returned_sha_or_reason):
      ("present", sha)  -- HTTP 200 and returned commit SHA equals the request
      ("absent", reason) -- HTTP 404/422: the commit does not exist on that remote
      ("error", reason)  -- permission/rate-limit/timeout/network/5xx or a returned
                            SHA that does not match the request: REMOTE VERIFICATION
                            FAILURE (fail closed; never reported as 'absent').
    Offline test mode: fixture = {"status": 200, "body": {"sha": ...}} | {"status": 404} ...
    """
    if fixture is not None:
        status = fixture.get("status")
        body = fixture.get("body") or {}
        if status == 200:
            returned = body.get("sha", "")
            if returned != sha:
                return ("error", f"returned commit SHA {returned!r} does not match requested {sha!r}")
            return ("present", returned)
        if status in (404, 422):
            return ("absent", f"HTTP {status}: commit not found on {repo}")
        return ("error", f"HTTP {status}: remote verification failure (fail closed)")
    url = f"https://api.github.com/repos/{repo}/commits/{sha}"
    req = urllib.request.Request(url, method="GET",
                                 headers={"Accept": "application/vnd.github+json",
                                          "X-GitHub-Api-Version": "2022-11-28"})
    if token:
        req.add_header("Authorization", f"token {token}")
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            body = json.loads(resp.read().decode("utf-8"))
            returned = body.get("sha", "")
            if returned != sha:
                return ("error", f"returned commit SHA {returned!r} does not match requested {sha!r}")
            return ("present", returned)
    except urllib.error.HTTPError as e:
        if e.code in (404, 422):
            return ("absent", f"HTTP {e.code}: commit not found on {repo}")
        return ("error", f"HTTP {e.code}: remote verification failure (fail closed)")
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        return ("error", f"network failure during remote verification (fail closed): {type(e).__name__}")
    except ValueError:
        return ("error", "unparseable API response (fail closed)")


# ---------------------------------------------------------------- EC-20 ----

def _pytest(repo: str, *argv: str):
    return subprocess.run([sys.executable, "-m", "pytest", *argv],
                          capture_output=True, text=True, cwd=repo, shell=False)


def collect_test_ids(repo: str, files: list[str] | None = None) -> tuple[set[str], str]:
    """Collect pytest node IDs, scoped to the named files when given.

    Round-19 EC-20 fix (runs 39417485854 / 39417488957 / 39417492355): repo-wide
    `pytest --collect-only` IMPORTS every test module under the root. Legacy
    script-style modules (gpu_worker/tests/*.py) execute their whole suite at
    import time and call sys.exit() at module level, which raises SystemExit
    inside pytest's collector -> INTERNALERROR -> collection aborts with
    "N test collected, M errors" and ZERO of the actually-named acceptance
    tests register as collectable, so the check failed closed with "named but
    not collectable" even though the named tests exist and pass.

    Collection is therefore scoped per-file: only the files that contain the
    named tests are imported. A named test whose file cannot be imported in
    this environment still fails closed (honest: the check cannot execute it
    here), but unrelated broken modules elsewhere in the tree cannot abort
    the collection of the named tests.

    Node IDs are normalized to repo-relative form: pytest displays node IDs
    relative to rootdir, so --rootdir is pinned to the repo and file args are
    passed absolute; a stray pyproject.toml/setup.cfg ABOVE the checkout can
    otherwise hijack rootdir and skew the displayed paths.
    """
    ids: set[str] = set()
    notes: list[str] = []
    repo_abs = os.path.abspath(repo)
    targets = files if files is not None else [None]
    for target in targets:
        arg = os.path.join(repo_abs, target) if target else None
        argv = ["--collect-only", "-q", "--rootdir", repo_abs] + ([arg] if arg else [])
        r = _pytest(repo, *argv)
        for line in r.stdout.splitlines():
            line = line.strip()
            if line.endswith("::") or not line or line.startswith(("=", "!", "no tests")):
                continue
            if "::" in line:
                ids.add(line.split()[0])
        if r.returncode != 0:
            scope = target or "<repo root>"
            notes.append(f"collection exit {r.returncode} for {scope}:\n" + r.stdout[-1500:])
    # defense in depth: if pytest still printed absolute paths, strip the repo prefix
    ids = {i[len(repo_abs) + 1:] if i.startswith(repo_abs + os.sep) else i for i in ids}
    return ids, "\n".join(notes)


def run_acceptance_tests(test_ids: list[str], repo: str, evidence_out: str | None):
    """Execute named acceptance tests on the exact candidate. Returns (findings, evidence)."""
    findings: list[str] = []
    started = datetime.now(timezone.utc)
    candidate = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                               cwd=repo).stdout.strip()
    # per-file scoped collection (see collect_test_ids docstring): bare paths
    # name their own file; node IDs name the file before the first '::'.
    named_files = sorted({(t.split("::", 1)[0] if "::" in t else t) for t in test_ids})
    absent_files = [f for f in named_files
                    if not os.path.isfile(os.path.join(repo, f))]
    for f in absent_files:
        findings.append(f"EC-20 named acceptance test file does not exist on the candidate: {f}")
    collectable_files = [f for f in named_files if f not in absent_files]
    collected, collect_note = collect_test_ids(repo, collectable_files)
    missing = [t for t in test_ids if t not in collected]
    for t in missing:
        findings.append(f"EC-20 acceptance test named but not collectable on the candidate: {t}")
    to_run = [t for t in test_ids if t not in missing]
    junit = os.path.join(os.path.dirname(os.path.abspath(evidence_out or "/tmp/e.json")),
                         f"junit_{int(time.time())}.xml") if evidence_out else \
        os.path.join("/tmp", f"junit_ec20_{int(time.time())}.xml")
    counts = {}
    exit_status = None
    executed_argv = []
    if to_run:
        executed_argv = [sys.executable, "-m", "pytest", *to_run, "--junitxml", junit, "-q", "--no-header"]
        r = subprocess.run(executed_argv, capture_output=True, text=True, cwd=repo, shell=False)
        exit_status = r.returncode
        try:
            root = ET.parse(junit).getroot()
            suite = root.find("testsuite")
            counts = {k: int(suite.get(k, 0)) for k in ("tests", "errors", "failures", "skipped")}
        except (OSError, ET.ParseError):
            counts = {}
            findings.append("EC-20 could not parse the junit report (fail closed)")
        if exit_status != 0:
            findings.append(f"EC-20 named acceptance tests FAILED on the candidate "
                            f"(pytest exit {exit_status}; counts={counts})")
    evidence = {
        "captured_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "candidate_sha": candidate,
        "command_argv": executed_argv,
        "shell_used": False,
        "python": sys.version.split()[0],
        "tests_requested": test_ids,
        "collection_scope_files": named_files,
        "tests_missing_from_collection": missing,
        "tests_executed": to_run,
        "junit_xml": junit,
        "counts": counts,
        "exit_status": exit_status,
        "collect_probe": collect_note,
    }
    if evidence_out:
        with open(evidence_out, "w", encoding="utf-8") as fh:
            json.dump(evidence, fh, indent=2)
    return findings, evidence


# ---------------------------------------------------------------- main ----

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--body", required=True, help="path to the PR body text")
    ap.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY", ""),
                    help="owner/name for API verification")
    ap.add_argument("--local", action="store_true", help="verify against local git objects (tests)")
    ap.add_argument("--require-directive-id-only", action="store_true")
    ap.add_argument("--check-tests", action="store_true",
                    help="EC-20: execute the named Acceptance-Tests on this candidate")
    ap.add_argument("--evidence-out", metavar="FILE", help="write the EC-20 evidence JSON here")
    ap.add_argument("--api-json", help="offline fixture: canned GitHub API response (tests)")
    ap.add_argument("--expect-head-sha", default="", metavar="SHA",
                    help="the EXACT current pull-request head SHA (passed explicitly by "
                         "the workflow). The PR body's Claimed-SHAs must contain it and "
                         "its Status-summary must cite a successful run on it.")
    ap.add_argument("--runs-manifest", default=None, metavar="FILE",
                    help="JSON {\"runs\": {\"<id>\": {head_sha, conclusion, event, "
                         "workflow_name}}} used to verify the cited Status-summary run "
                         "offline (tests)")
    ap.add_argument("--required-workflow", action="append", default=[],
                    metavar="NAME",
                    help="a workflow whose successful run on the exact head is REQUIRED "
                         "evidence (repeatable). Default: the two required checks of this "
                         "repository, the Toscanini Stability Gate and the Governance "
                         "Acceptance. A successful run of any other workflow is not "
                         "required-check evidence.")
    ap.set_defaults(required_workflow=[])
    args = ap.parse_args()
    if not args.required_workflow:
        args.required_workflow = list(DEFAULT_REQUIRED_WORKFLOWS)

    if not os.path.exists(args.body):
        print("FAIL: PR body file missing (fail closed)")
        return 1
    body = open(args.body, encoding="utf-8").read()
    findings: list[str] = []

    if args.require_directive_id_only:
        for f in ("Directive-ID", "Closes-items", "Status-summary"):
            if not any(v for v in field_values(body, f)):
                findings.append(f"EC-19 {f}: field missing or empty")
        print("\n".join(findings) if findings else "PASS: directive fields present")
        return 1 if findings else 0

    if not any(v for v in field_values(body, "Directive-ID")):
        findings.append("EC-19 Directive-ID: field missing or empty")

    if args.check_tests:
        tests = [t for v in field_values(body, "Acceptance-Tests") for t in v.split()]
        if not tests:
            findings.append("EC-20 Acceptance-Tests: field missing or empty (fail closed)")
        else:
            tf, ev = run_acceptance_tests(tests, os.getcwd(), args.evidence_out)
            findings.extend(tf)
            if not tf:
                print(f"OK  : EC-20 executed {len(ev['tests_executed'])} named test(s) on "
                      f"candidate {ev['candidate_sha'][:12]} counts={ev['counts']}")

    directive_only = args.require_directive_id_only
    if not directive_only:
        shas_field = field_values(body, "Claimed-SHAs")
        explicit_none = any(v.strip().lower() == "none" for v in shas_field)
        other_tokens = [t for v in shas_field for t in v.replace(",", " ").split()
                        if t.lower() != "none"]
        if explicit_none and other_tokens:
            findings.append("EC-16 Claimed-SHAs mixes the explicit 'none' declaration with "
                            "SHA entries (fail closed)")
            shas, parse_errors = [], []
        elif explicit_none:
            print("OK  : Claimed-SHAs: none (explicit no-claims declaration)")
            shas, parse_errors = [], []
        else:
            shas, parse_errors = parse_claimed_shas(shas_field)
        findings.extend(f"EC-16 {e}" for e in parse_errors)
        if not parse_errors and shas:
            token = os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
            fixture = json.load(open(args.api_json, encoding="utf-8")) if args.api_json else None
            for sha in sorted(set(shas)):
                if args.local:
                    status, detail = verify_sha_local(sha)
                else:
                    if not args.repo:
                        findings.append("EC-16 no --repo/GITHUB_REPOSITORY given for remote "
                                        "verification (fail closed)")
                        break
                    status, detail = verify_sha_api(args.repo, sha, token, fixture)
                if status == "present":
                    print(f"OK  : {sha} verified on remote ({detail})")
                elif status == "absent":
                    findings.append(f"EC-16 claimed SHA not on remote: {sha} ({detail})")
                else:
                    findings.append(f"EC-16 VERIFICATION_FAILURE for {sha} ({detail}); "
                                    f"fail closed -- not treated as absent")
        elif not shas_field:
            findings.append("EC-16 Claimed-SHAs: field missing (fail closed; use the literal "
                            "value 'none' to declare explicitly that a PR claims no results)")

    # ---- Round-24 direction 7 / Round-25 direction 3: bind the body to the
    #      EXACT current head AND to the REQUIRED WORKFLOW evidence.
    if args.expect_head_sha and not directive_only:
        head = args.expect_head_sha.strip()
        required = [w.strip() for w in args.required_workflow if w.strip()]
        if not FULL_SHA.match(head):
            findings.append("EC-16 --expect-head-sha must be a full 40-character lowercase "
                            "hex commit SHA (fail closed)")
        else:
            claimed = {t for v in field_values(body, "Claimed-SHAs")
                       for t in v.replace(",", " ").split()}
            if head not in claimed:
                findings.append(
                    f"EC-16 Claimed-SHAs does not declare the exact current pull-request "
                    f"head {head} (declared: {sorted(claimed) if claimed else 'none'}) -- "
                    f"the PR body must describe the candidate it is proposing (fail "
                    f"closed; a stale-but-real SHA is not evidence for this head)")
            cited = []
            for line in body.splitlines():
                m = STATUS_LINE.match(line)
                if not m:
                    continue
                for rid, sha in RUNPASS.findall(m.group(1)):
                    cited.append((rid, sha))
            token = os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
            manifest = None
            if args.runs_manifest:
                try:
                    manifest = json.load(open(args.runs_manifest, encoding="utf-8"))
                except (OSError, ValueError) as e:
                    print(f"FAIL: runs manifest unusable ({e}); fail closed")
                    return 1
            good = []
            satisfied = set()
            for rid, sha in cited:
                if sha != head:
                    # a citation for another commit is a STALE citation: it can
                    # never be evidence for this candidate
                    findings.append(
                        f"EC-22 Status-summary cites run_id={rid} for {sha[:12]}..., which "
                        f"is not the exact current head {head[:12]}... (a real but STALE "
                        f"run is not evidence for the current head)")
                    continue
                if manifest is not None:
                    rec = (manifest.get("runs") or {}).get(rid)
                    if rec is None:
                        findings.append(
                            f"EC-22 Status-summary cites run_id={rid}, which is not in the "
                            f"runs manifest (unverifiable -- fail closed)")
                        continue
                    if rec.get("head_sha") != head:
                        findings.append(
                            f"EC-22 Status-summary cites run_id={rid} for {sha[:12]}... but "
                            f"that run's head_sha is {str(rec.get('head_sha'))[:12]}... "
                            f"(wrong commit)")
                        continue
                    if rec.get("conclusion") != "success":
                        findings.append(
                            f"EC-22 Status-summary cites run_id={rid} whose conclusion is "
                            f"{rec.get('conclusion')!r}, not success")
                        continue
                    name = rec.get("workflow_name") or rec.get("name")
                else:
                    if not args.repo:
                        findings.append("EC-22 no --repo/GITHUB_REPOSITORY given to verify the "
                                        "cited Status-summary run (fail closed)")
                        break
                    status, rec = api_run(args.repo, rid, token)
                    if status != "present":
                        findings.append(f"EC-22 Status-summary cites run_id={rid}, which "
                                        f"could not be verified ({rec}) -- fail closed")
                        continue
                    if rec.get("head_sha") != head:
                        findings.append(
                            f"EC-22 Status-summary cites run_id={rid} for {sha[:12]}... but "
                            f"that run's head_sha is {str(rec.get('head_sha'))[:12]}... "
                            f"(a real but STALE run is not evidence for the current head)")
                        continue
                    if rec.get("conclusion") != "success":
                        findings.append(
                            f"EC-22 Status-summary cites run_id={rid} whose conclusion is "
                            f"{rec.get('conclusion')!r}, not success")
                        continue
                    name = rec.get("name")
                if not name:
                    findings.append(
                        f"EC-22 Status-summary cites run_id={rid}, whose workflow identity is "
                        f"unknown (the run record carries no workflow name) -- the run's "
                        f"identity cannot be established, so it cannot be required-check "
                        f"evidence (fail closed)")
                    continue
                if required and name not in required:
                    findings.append(
                        f"EC-22 Status-summary cites run_id={rid} of workflow {name!r}, which "
                        f"is not one of the required workflows {required} -- a successful "
                        f"unrelated workflow on the right commit is not required-check "
                        f"evidence")
                    continue
                good.append(rid)
                satisfied.add(name)
            if not good:
                findings.append(
                    f"EC-22 Status-summary cites no successful run on the exact current head "
                    f"{head} (cited RUN_PASS lines: {cited or 'none'}) -- a stale but real "
                    f"historical run id is not evidence for this candidate (fail closed)")
            for wf in required:
                if wf not in satisfied:
                    findings.append(
                        f"EC-22 Status-summary cites no successful run of the required "
                        f"workflow {wf!r} on the exact current head {head} -- the required "
                        f"check evidence must be identified explicitly, not inferred from a "
                        f"generic successful run (fail closed)")
            if good:
                print(f"OK  : Status-summary cites successful run(s) {sorted(set(good))} of "
                      f"the required workflow(s) {sorted(satisfied)} on the exact current head "
                      f"{head}")

    print("\n".join(findings) if findings else "")
    if not findings:
        print("PASS: EC-16/19/20 claim verification clean")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
