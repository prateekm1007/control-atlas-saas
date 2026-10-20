#!/usr/bin/env python3
"""EC-23 secret-scan control verification (Round-20 hardening; Round-21 disclosure fix).

Auditor finding (Round-20 direction 3): the previous positive control treated
ANY nonzero gitleaks exit as proof that the seeded finding was detected, so a
scanner crash or invocation error could masquerade as a successful control.
The previous all-refs classification also tolerated every generic-api-key
finding under docs/ or evidence/ (fixed by secret_baseline_check.py).

Auditor finding (Round-21 direction 3 -- CRITICAL): the Round-20 success path
printed the report's `Secret` field value, so a FAILED redaction could
disclose the very secret the control detected. This version NEVER prints
`Secret` or `Match` field values or any raw report/log bytes: diagnostics
carry safe metadata only (rule, file, line, counts, booleans). Additional
fail-closed rules:

  * if --log-file is requested, the file MUST exist and be readable --
    a missing or unreadable capture log fails the control (Round-21 3.3);
  * a clean HEAD-tree (or clean-fixture) negative result requires the report
    to exist, parse as a JSON array, and contain ZERO findings -- scanner
    exit status alone never proves cleanliness (Round-21 3.4).

Modes:

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

  negative mode  -- a clean repository/tree MUST pass:
      - exit code must be 0; a nonzero code is a scanner failure or a false
        positive, both of which fail the control;
      - the JSON report must exist, parse, and contain zero findings;
      - if --log-file is requested it must exist and be readable (fail
        closed), and if --key-env is provided the key must be absent from it.

Every failure path fails closed (exit 1); only a fully verified control exits 0.
No path in this module writes secret material to stdout or stderr.
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
        return None, f"report is not parseable JSON ({type(e).__name__})"
    if not isinstance(data, list):
        return None, "report is not a JSON array of findings"
    return data, None


def _read_log_bytes(log_path: str, fails: list[str]) -> bytes | None:
    """Round-21 3.3: a requested capture log must exist and be readable (fail closed)."""
    if not log_path:
        return None
    if not os.path.exists(log_path):
        fails.append(f"captured scanner log {log_path!r} is MISSING -- a requested "
                     f"log must exist and be checked (fail closed)")
        return None
    try:
        return open(log_path, "rb").read()
    except OSError as e:
        fails.append(f"captured scanner log {log_path!r} is unreadable "
                     f"({type(e).__name__}) -- fail closed")
        return None


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
        # Safe metadata ONLY: rule/file/line/count. Secret and Match values are
        # never printed (Round-21 critical fix).
        print("expected finding verified: rule=%s file=%s count=%d line(s)=%s"
              % (expect_rule, expect_file, len(expected),
                 ",".join(str(f.get("StartLine")) for f in expected)))
    if others:
        fails.append("unexpected additional finding(s): %s"
                     % "; ".join("%s@%s:%s" % (f.get("RuleID"), f.get("File"),
                                               f.get("StartLine")) for f in others))
    key = os.environ.get(key_env, "")
    if not key:
        fails.append(f"the seeded key is not available via ${key_env} "
                     f"(redaction cannot be proven -- fail closed)")
    else:
        raw = open(report_path, "rb").read() if report_path and os.path.exists(report_path) else b""
        key_bytes = key.encode()
        report_redacted = key_bytes not in raw
        field_redacted = all(f.get("Secret") != key and f.get("Match") != key
                             for f in expected)
        if not report_redacted:
            fails.append("the raw seeded key appears in the scan report bytes "
                         "(redaction FAILED)")
        if not field_redacted:
            fails.append("a finding's Secret or Match field carries the raw key "
                         "(redaction FAILED)")
        if report_redacted and field_redacted:
            print("redaction verified: the raw key is absent from all %d report "
                  "bytes and from every finding's Secret/Match field" % len(raw))
        else:
            # Round-21 3.1/3.2: disclose ONLY the boolean outcome and the field
            # NAMES that failed; never the values.
            print("redaction FAILED: report_bytes_clean=%s fields_clean=%s "
                  "(values are never printed)" % (report_redacted, field_redacted))
        log = _read_log_bytes(log_path, fails) if log_path else None
        if log is not None and key_bytes in log:
            fails.append("the raw seeded key appears in the captured scanner log "
                         "(redaction FAILED)")
        elif log is not None:
            print("capture log verified: the raw key is absent from all %d log "
                  "bytes" % len(log))
    if fails:
        for f in fails:
            print("FAIL: " + f)
        return 1
    print("PASS: positive control -- seeded synthetic key genuinely detected by "
          "rule %s in %s, report well-formed, key redacted everywhere it was "
          "checked (values never printed)" % (expect_rule, expect_file))
    return 0


def verify_negative(report_path: str, exit_code: int,
                    log_path: str | None = None, key_env: str = "SEEDED_KEY") -> int:
    fails: list[str] = []
    if exit_code != 0:
        fails.append(f"clean scan exited {exit_code} (expected 0): "
                     f"scanner failure or false positive -- control fails")
    # Round-21 3.4: exit 0 is NOT sufficient; the report must exist, parse, and
    # conform to the empty-finding schema (a JSON array with zero entries).
    data, err = _load_report(report_path)
    if err:
        fails.append(err)
    elif data:
        fails.append("clean scan produced %d finding(s): %s"
                     % (len(data), "; ".join("%s@%s" % (f.get("RuleID"),
                                                        f.get("File")) for f in data)))
    else:
        print("negative report verified: exists, parses, zero findings "
              "(exit status alone was not trusted)")
    if log_path:
        log = _read_log_bytes(log_path, fails)
        if log is not None:
            key = os.environ.get(key_env, "")
            if key:
                if key.encode() in log:
                    fails.append("the raw key appears in the captured scanner log "
                                 "(redaction FAILED)")
                else:
                    print("capture log verified: the raw key is absent from all "
                          "%d log bytes" % len(log))
            else:
                print("capture log verified: exists and is readable (%d bytes)" % len(log))
    if fails:
        for f in fails:
            print("FAIL: " + f)
        return 1
    print("PASS: negative control -- clean target scanned with a verified "
          "empty-finding report")
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
                    help="captured scanner stdout/stderr; when given it MUST exist "
                         "and be readable (fail closed) and is checked for absence "
                         "of the raw key")
    a = ap.parse_args()
    if a.mode == "positive":
        return verify_positive(a.report, a.exit_code, a.expect_rule,
                               a.expect_file, a.key_env, a.log_file)
    return verify_negative(a.report, a.exit_code, a.log_file, a.key_env)


if __name__ == "__main__":
    sys.exit(main())
