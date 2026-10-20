#!/usr/bin/env python3
"""EC-16/19/20 claim verification (inert until adopted; called by proposals/ci/verify-claimed-shas.yml).

Modes (combinable except where noted):
  default            : extract `Claimed-SHAs:` from the PR body; verify each full 40-hex SHA exists
                       on the remote. Remote verification uses the GitHub REST API
                       (GET /repos/{repo}/commits/<sha>) -- NOT `git ls-remote`, which resolves only
                       ref tips and would falsely reject mid-branch commits (recorded correction).
  --local            : offline mode; verify SHAs against the local git object store (tests only).
  --require-directive-id-only : EC-19 check only: body must carry a non-empty `Directive-ID:` field.
  --check-tests-in-diff FILE  : EC-20: every identifier in `Acceptance-Tests:` must appear in the
                       diff file list FILE (one path per line) or the PR fails.

Exit codes: 0 = pass; 1 = fail (findings printed to stdout, one per line, `FAIL:` prefixed).
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import urllib.error
import urllib.request

FULL_SHA = re.compile(r"\b[0-9a-f]{40}\b")
FIELDS = ("Directive-ID", "Closes-items", "Base-SHA", "Stack-parent",
          "Claimed-SHAs", "Acceptance-Tests", "Status-summary")


def field_values(body: str, field: str) -> list[str]:
    """Return the values of a template field: the remainder of its line plus following indented lines."""
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


def sha_exists_local(sha: str) -> bool:
    r = subprocess.run(["git", "cat-file", "-e", f"{sha}^{{commit}}"],
                       capture_output=True)
    return r.returncode == 0


def sha_exists_api(repo: str, sha: str, token: str | None) -> bool:
    url = f"https://api.github.com/repos/{repo}/commits/{sha}"
    req = urllib.request.Request(url, method="GET")
    if token:
        req.add_header("Authorization", f"token {token}")
    try:
        with urllib.request.urlopen(req) as resp:
            return resp.status == 200
    except urllib.error.HTTPError as e:
        if e.code in (404, 422):
            return False
        raise


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--body", required=True, help="path to the PR body text")
    ap.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY", ""),
                    help="owner/name for API verification")
    ap.add_argument("--local", action="store_true", help="verify against local git objects")
    ap.add_argument("--require-directive-id-only", action="store_true")
    ap.add_argument("--check-tests-in-diff", metavar="FILE")
    args = ap.parse_args()

    body = open(args.body, encoding="utf-8").read()
    findings: list[str] = []

    if args.require_directive_id_only:
        vals = field_values(body, "Directive-ID")
        if not any(v for v in vals):
            findings.append("FAIL: EC-19 Directive-ID: field missing or empty")
        for f in ("Closes-items", "Status-summary"):
            if not any(v for v in field_values(body, f)):
                findings.append(f"FAIL: EC-19 {f}: field missing or empty")
        print("\n".join(findings) if findings else "PASS: directive fields present")
        return 1 if findings else 0

    directive = any(v for v in field_values(body, "Directive-ID"))
    if not directive:
        findings.append("FAIL: EC-19 Directive-ID: field missing or empty")

    if args.check_tests_in_diff:
        tests = [t for v in field_values(body, "Acceptance-Tests") for t in v.split()]
        diff_paths = set(open(args.check_tests_in_diff, encoding="utf-8").read().split())
        for t in tests:
            if not any(t in p for p in diff_paths):
                findings.append(f"FAIL: EC-20 acceptance test named but not in diff: {t}")
        if not tests:
            findings.append("FAIL: EC-20 Acceptance-Tests: field missing or empty")
        print("\n".join(findings) if findings else "PASS: all named acceptance tests in diff")
        return 1 if findings else 0

    shas = []
    for v in field_values(body, "Claimed-SHAs"):
        shas.extend(FULL_SHA.findall(v))
    if not shas:
        findings.append("FAIL: EC-16 Claimed-SHAs: contains no full 40-hex SHA")
    else:
        token = os.environ.get("GH_TOKEN")
        for sha in sorted(set(shas)):
            ok = sha_exists_local(sha) if args.local else sha_exists_api(args.repo, sha, token)
            if ok:
                print(f"OK  : {sha}")
            else:
                findings.append(f"FAIL: EC-16 claimed SHA not on remote: {sha}")

    print("\n".join(findings))
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
