#!/usr/bin/env python3
"""Package-scoped forbidden-reference check (Round-19 direction 7).

Encodes the auditor's Round-18 documentation-exemption decision as an
automated, fail-closed check: references to the removed external provider
(gemini / generativeai / the AIzaSy key prefix) are permitted ONLY in
documentation prose under docs/ and in the enumerated baseline exceptions
below. Any match anywhere else in the tree FAILS.

Pattern (identical to the auditor's acceptance grep):
    gemini | generativeai | AIzaSy   (case-insensitive)

Enumerated baseline exceptions (auditor-accepted, each documented):
  - docs/**                              documentation exemption (the pattern
                                          text itself is quoted in the exempt
                                          RECEIPT R10 evidence entry);
  - evidence/incident-push-2/**          pre-existing incident prose on the
                                          accepted baseline (8 instances,
                                          out of package scope, untouched);
  - scripts/ci/forbidden_reference_check.py
                                         this checker: the pattern tokens
                                          appear here as data. The file is
                                          digest-locked by the package
                                          manifest, so it cannot drift
                                          silently.

Implementation uses `git grep -nIiE` (binary-skipping, tracked files) with
pathspec excludes, exactly reproducing the acceptance grep that a human ran
in Rounds 18-19, so machine and human evidence cannot diverge.

Exit: 0 = clean; 1 = findings; 2 = usage/internal error.
"""
import argparse
import subprocess
import sys

PATTERN = "gemini|generativeai|AIzaSy"
EXCLUDES = [
    ":(exclude)docs",
    ":(exclude)evidence/incident-push-2",
    ":(exclude)scripts/ci/forbidden_reference_check.py",
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", default=".")
    a = ap.parse_args()
    argv = ["git", "grep", "-nIiE", PATTERN, "--", "."] + EXCLUDES
    r = subprocess.run(argv, capture_output=True, text=True, cwd=a.repo_root)
    if r.returncode not in (0, 1):
        print(f"FAIL: git grep unusable (exit {r.returncode}): {r.stderr.strip()}")
        return 2
    hits = [l for l in r.stdout.splitlines() if l.strip()]
    if hits:
        print("FAIL: non-documentation provider reference(s) found "
              "(documentation-exemption decision violated):")
        for h in hits:
            print(f"  {h}")
        return 1
    print("PASS: no provider references outside the documented exemptions "
          "(docs/ prose, incident-push-2 baseline prose, this checker's own pattern data)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
