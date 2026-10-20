#!/usr/bin/env python3
"""EC-21 declared base (inert until adopted). Auditor revision 2026-10-20 (item 3.3).

Previous defect: this checker only proved that the base object exists in the LOCAL object
store and is an ancestor of HEAD. That is not proof the base is on the intended remote, and
`Stack-parent` was never validated. The revised flow verifies, in order:

  1. The declared base is a full 40-hex SHA (malformed -> FAIL).
  2. The declared Stack-parent is a known stack class (main-line | science-line | policy-line).
  3. The repository history is COMPLETE (a shallow checkout fails closed; CI must fetch
     with fetch-depth: 0 / --unshallow).
  4. LOCAL object availability (distinct, informational label: LOCAL_OBJECT_ABSENT).
  5. LOCAL ancestry: HEAD descends from the base (unrelated base -> FAIL).
  6. REMOTE verification: GET /repos/{slug}/commits/{base}; the returned commit SHA must
     EQUAL the request. 404/422 -> FAIL NOT_ON_REMOTE (a real verification result).
     401/403/429/5xx/timeout/network or a mismatched returned SHA -> FAIL
     REMOTE_VERIFICATION_FAILURE (fail closed; never reported as "not on remote").
  7. STACK validation (line-unique history). Neither plain ancestry nor first-parent
     chains can discriminate stack lines in a repo whose lines merge each other: after
     science merged main, the main tip is BOTH an ancestor of the science head AND on
     its first-parent chain, yet work stacked there is NOT on the science line (the
     exact H-1 failure this rule exists for). The implemented semantics:
       main-line    : base must be reachable from origin/main.
       science-line : base must be reachable from origin/science/p0-correctness-2026-10-19
                      AND NOT reachable from origin/main (science-unique history).
       policy-line  : base must be reachable from origin/policy/engineering-constitution
                      AND NOT reachable from origin/main (policy-unique history).
     A base grabbed from the wrong line's tip FAILs (STACK_MISMATCH).

Offline test mode: --api-json FILE supplies a canned API response
  {"status": 200, "body": {"sha": "<40hex>"}} | {"status": 404} | {"status": 403} ...
and --stack-refs JSON overrides the stack-class -> local ref mapping. Production mode
(no fixtures) uses the live GitHub API and origin/<branch> refs.

Exit codes: 0 = pass; 1 = fail (distinct FAIL labels above); 2 = usage error.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.request

FULL_SHA = re.compile(r"^[0-9a-f]{40}$")

DEFAULT_STACK_REFS = {
    "main-line": "origin/main",
    "science-line": "origin/science/p0-correctness-2026-10-19",
    "policy-line": "origin/policy/engineering-constitution",
}


def git_ok(*argv: str) -> bool:
    return subprocess.run(["git", *argv], capture_output=True).returncode == 0


def verify_remote_base(slug: str, sha: str, token: str | None, fixture: dict | None):
    """Return (status, detail): 'present' / 'absent' / 'error'."""
    if fixture is not None:
        status = fixture.get("status")
        body = fixture.get("body") or {}
        if status == 200:
            returned = body.get("sha", "")
            if returned != sha:
                return ("error", f"returned commit SHA {returned!r} does not match requested {sha!r}")
            return ("present", returned)
        if status in (404, 422):
            return ("absent", f"HTTP {status}: commit not found on {slug}")
        return ("error", f"HTTP {status}: remote verification failure (fail closed)")
    url = f"https://api.github.com/repos/{slug}/commits/{sha}"
    req = urllib.request.Request(url, method="GET",
                                 headers={"Accept": "application/vnd.github+json",
                                          "X-GitHub-Api-Version": "2022-11-28"})
    if token:
        req.add_header("Authorization", f"token {token}")
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            returned = json.loads(resp.read().decode("utf-8")).get("sha", "")
            if returned != sha:
                return ("error", f"returned commit SHA {returned!r} does not match requested {sha!r}")
            return ("present", returned)
    except urllib.error.HTTPError as e:
        if e.code in (404, 422):
            return ("absent", f"HTTP {e.code}: commit not found on {slug}")
        return ("error", f"HTTP {e.code}: remote verification failure (fail closed)")
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        return ("error", f"network failure during remote verification (fail closed): {type(e).__name__}")
    except ValueError:
        return ("error", "unparseable API response (fail closed)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", help="declared base SHA (or via --declaration)")
    ap.add_argument("--head", default="HEAD")
    ap.add_argument("--stack-parent", help="main-line | science-line | policy-line (or via --declaration)")
    ap.add_argument("--declaration", help="PR-body-style file with Base-SHA:/Stack-parent: fields")
    ap.add_argument("--repo-slug", default=os.environ.get("GITHUB_REPOSITORY", ""),
                    help="owner/name for the remote commit lookup")
    ap.add_argument("--api-json", help="offline fixture: canned GitHub API response (tests)")
    ap.add_argument("--stack-refs", help="offline JSON mapping stack class -> local ref (tests)")
    a = ap.parse_args()

    base, stack = a.base, a.stack_parent
    if a.declaration:
        text = open(a.declaration, encoding="utf-8").read()
        m = re.search(r"^Base-SHA:\s*(\S+)", text, re.M)
        if m and not base:
            base = m.group(1)
        m = re.search(r"^Stack-parent:\s*(\S+)", text, re.M)
        if m and not stack:
            stack = m.group(1)

    if not base:
        print("FAIL: Base-SHA missing (fail closed)")
        return 1
    if not FULL_SHA.match(base):
        print(f"FAIL: MALFORMED_BASE {base!r} is not a full 40-character lowercase hex SHA")
        return 1
    if not stack:
        print("FAIL: Stack-parent missing (fail closed)")
        return 1

    stack_refs = DEFAULT_STACK_REFS
    if a.stack_refs:
        stack_refs = json.load(open(a.stack_refs, encoding="utf-8"))
    if stack not in stack_refs:
        print(f"FAIL: UNKNOWN_STACK_CLASS {stack!r} (allowed: {sorted(stack_refs)})")
        return 1

    shallow = subprocess.run(["git", "rev-parse", "--is-shallow-repository"],
                             capture_output=True, text=True).stdout.strip()
    if shallow == "true":
        print("FAIL: INCOMPLETE_HISTORY shallow checkout; fetch complete history "
              "(fetch-depth: 0 / git fetch --unshallow) before declaring a base")
        return 1

    if not git_ok("cat-file", "-e", f"{base}^{{commit}}"):
        print(f"FAIL: LOCAL_OBJECT_ABSENT base {base} is not in the local object store "
              f"(fetch it; note: local absence is NOT remote verification)")
        return 1
    if not git_ok("merge-base", "--is-ancestor", base, a.head):
        print(f"FAIL: ANCESTRY HEAD ({a.head}) does not descend from declared base {base} "
              f"(unrelated base)")
        return 1

    fixture = json.load(open(a.api_json, encoding="utf-8")) if a.api_json else None
    status, detail = verify_remote_base(a.repo_slug, base,
                                        os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN"),
                                        fixture)
    if status == "absent":
        print(f"FAIL: NOT_ON_REMOTE base {base} is not present on {a.repo_slug} ({detail})")
        return 1
    if status == "error":
        print(f"FAIL: REMOTE_VERIFICATION_FAILURE base {base} ({detail}); fail closed "
              f"-- not treated as present or absent")
        return 1

    stack_ref = stack_refs[stack]
    on_line = git_ok("merge-base", "--is-ancestor", base, stack_ref)
    if stack == "main-line":
        stack_ok = on_line
        why = f"reachable from {stack_ref}"
    else:
        on_main = git_ok("merge-base", "--is-ancestor", base, "origin/main")
        stack_ok = on_line and not on_main
        why = (f"on {stack_ref} and NOT reachable from origin/main "
               f"(line-unique history; merged-in main-side commits are not this line)")
    if not stack_ok:
        print(f"FAIL: STACK_MISMATCH base {base} does not satisfy the declared stack "
              f"parent {stack} (expected: {why})")
        return 1

    print(f"PASS: EC-21 base {base} is local, on remote ({detail}), an ancestor of {a.head}, "
          f"and valid for the declared stack parent {stack}: {why}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
