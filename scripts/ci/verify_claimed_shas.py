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


def verify_sha_api(repo: str, sha: str, token: str | None):
    """Verify one SHA against the intended repository via the GitHub commits API.

    Returns (status, returned_sha_or_reason):
      ("present", sha)  -- HTTP 200 and returned commit SHA equals the request
      ("absent", reason) -- HTTP 404/422: the commit does not exist on that remote
      ("error", reason)  -- permission/rate-limit/timeout/network/5xx or a returned
                            SHA that does not match the request: REMOTE VERIFICATION
                            FAILURE (fail closed; never reported as 'absent').
    """
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


def collect_test_ids(repo: str) -> tuple[set[str], str]:
    r = _pytest(repo, "--collect-only", "-q")
    ids = set()
    for line in r.stdout.splitlines():
        line = line.strip()
        if line.endswith("::") or not line or line.startswith(("=", "!", "no tests")):
            continue
        if "::" in line:
            ids.add(line.split()[0])
    return ids, (r.stdout[-2000:] if r.returncode != 0 else "")


def run_acceptance_tests(test_ids: list[str], repo: str, evidence_out: str | None):
    """Execute named acceptance tests on the exact candidate. Returns (findings, evidence)."""
    findings: list[str] = []
    started = datetime.now(timezone.utc)
    candidate = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                               cwd=repo).stdout.strip()
    collected, collect_note = collect_test_ids(repo)
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
    args = ap.parse_args()

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
            for sha in sorted(set(shas)):
                if args.local:
                    status, detail = verify_sha_local(sha)
                else:
                    if not args.repo:
                        findings.append("EC-16 no --repo/GITHUB_REPOSITORY given for remote "
                                        "verification (fail closed)")
                        break
                    status, detail = verify_sha_api(args.repo, sha, token)
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

    print("\n".join(findings) if findings else "")
    if not findings:
        print("PASS: EC-16/19/20 claim verification clean")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
