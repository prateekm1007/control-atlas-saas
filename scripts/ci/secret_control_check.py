#!/usr/bin/env python3
"""EC-23 secret-scan control verification (Round-20 hardening).

Auditor finding (Round-20 direction 3): the previous positive control treated
ANY nonzero gitleaks exit as proof that the seeded finding was detected, so a
scanner crash or invocation error could masquerade as a successful control.
The previous all-refs classification also tolerated every generic-api-key
finding under docs/ or evidence/ (fixed by secret_baseline_check.py).

This module verifies the CONTROLS of a scan run:

  positive mode  -- the seeded synthetic key MUST be genuinely detected:
      - exit code must equal the reserved FINDINGS code (101; the scan is
        invoked with --exit-code 101 so "findings found" is distinct from
        execution errors, which exit 1);
      - the JSON report must exist and parse;
      - the report must contain the expected finding for the seeded file with
        the EXPECTED RULE (aws-access-token) actually fired;
      - the report must contain no unexpected additional findings;
      - the assembled key value must be REDACTED: absent from the report
        bytes and from the captured scanner log (the key is passed via an
        environment variable name, never via argv or file content).

  negative mode  -- a clean repository MUST pass:
      - exit code must be 0; a nonzero code is a scanner failure or a false
        positive, both of which fail the control;
      - the JSON report must exist, parse, and contain zero findings.

Every failure path fails closed (exit 1); only a fully verified control exits 0.
"""
import argparse
import json
import os
import sys

FINDINGS_EXIT = 101  # reserved: --exit-code 101 means "leaks found"; errors exit 1


def _load_report(path: str):
    if not path or not os.path.exists(path):
        return None, "report file absent (scanner produced no report)"
    try:
        data = json.load(open(path, encoding="utf-8"))
    except (OSError, ValueError) as e:
        return None, f"report is not parseable JSON ({e})"
    if not isinstance(data, list):
        return None, "report is not a JSON array of findings"
    return data, None


def verify_positive(report_path: str, exit_code: int, expect_rule: str,
                    expect_file: str, key_env: str, log_path: str | None) -> int:
    fails: list[str] = []
    if exit_code == 0:
        fails.append(f"scanner exited 0: the seeded key was NOT detected "
                     f"(exit_code={exit_code}; expected {FINDINGS_EXIT})")
    elif exit_code != FINDINGS_EXIT:
        fails.append(f"scanner exit code {exit_code} is neither 0 nor the reserved "
                     f"findings code {FINDINGS_EXIT} -- execution error masquerading "
                     f"as detection is rejected (fail closed)")
    data, err = _load_report(report_path)
    if err:
        fails.append(err)
        data = []
    expected = [f for f in data
                if f.get("RuleID") == expect_rule and f.get("File") == expect_file]
    others = [f for f in data if f not in expected]
    if not expected:
        fails.append(f"no finding for rule {expect_rule!r} in file {expect_file!r} "
                     f"-- the expected rule did not fire")
    else:
        print("expected finding verified: rule=%s file=%s count=%d line(s)=%s"
              % (expect_rule, expect_file, len(expected),
                 ",".join(str(f.get("StartLine")) for f in expected)))
    if others:
        fails.append("unexpected additional finding(s): %s"
                     % "; ".join("%s@%s:%s" % (f.get("RuleID"), f.get("File"),
                                               f.get("StartLine")) for f in others))
    raw = open(report_path, "rb").read() if report_path and os.path.exists(report_path) else b""
    key = os.environ.get(key_env, "")
    if not key:
        fails.append(f"the seeded key is not available via ${key_env} "
                     f"(redaction cannot be proven -- fail closed)")
    else:
        if key.encode() in raw:
            fails.append("the raw seeded key appears in the scan report "
                         "(redaction FAILED)")
        for f in expected:
            if f.get("Secret") == key:
                fails.append("the finding's Secret field carries the raw key "
                             "(redaction FAILED)")
        print("redaction verified: report Secret/Match fields carry %r; raw key "
              "absent from %d report bytes" % (expected[0].get("Secret", "<none>")
                                               if expected else "<none>", len(raw)))
    if log_path and os.path.exists(log_path) and key:
        if key.encode() in open(log_path, "rb").read():
            fails.append("the raw seeded key appears in the captured scanner log "
                         "(redaction FAILED)")
    if fails:
        for f in fails:
            print("FAIL: " + f)
        return 1
    print("PASS: positive control -- seeded synthetic key genuinely detected by "
          "rule %s in %s, report well-formed, key redacted" % (expect_rule, expect_file))
    return 0


def verify_negative(report_path: str, exit_code: int) -> int:
    fails: list[str] = []
    if exit_code != 0:
        fails.append(f"clean repository scan exited {exit_code} (expected 0): "
                     f"scanner failure or false positive -- control fails")
    data, err = _load_report(report_path)
    if err:
        fails.append(err)
    elif data:
        fails.append("clean repository produced %d finding(s): %s"
                     % (len(data), "; ".join("%s@%s" % (f.get("RuleID"),
                                                        f.get("File")) for f in data)))
    if fails:
        for f in fails:
            print("FAIL: " + f)
        return 1
    print("PASS: negative control -- clean repository scanned with zero findings")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=["positive", "negative"])
    ap.add_argument("--report", required=True, help="gitleaks JSON report path")
    ap.add_argument("--exit-code", type=int, required=True,
                    help="the scanner's observed exit code")
    ap.add_argument("--expect-rule", default="aws-access-token",
                    help="rule that must fire for the seeded key (positive mode)")
    ap.add_argument("--expect-file", default="fixture.txt",
                    help="file that must carry the seeded finding (positive mode)")
    ap.add_argument("--key-env", default="SEEDED_KEY",
                    help="environment variable holding the assembled synthetic key "
                         "(value never printed, never passed via argv)")
    ap.add_argument("--log-file", default=None,
                    help="captured scanner stdout/stderr to additionally check for "
                         "absence of the raw key")
    a = ap.parse_args()
    if a.mode == "positive":
        return verify_positive(a.report, a.exit_code, a.expect_rule,
                               a.expect_file, a.key_env, a.log_file)
    return verify_negative(a.report, a.exit_code)


if __name__ == "__main__":
    sys.exit(main())
