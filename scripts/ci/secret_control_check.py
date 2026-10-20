#!/usr/bin/env python3
"""EC-23 secret-scan control verification (Round-20 hardening; R21/R22 disclosure fixes).

Auditor finding (Round-20 direction 3): the previous positive control treated
ANY nonzero gitleaks exit as proof that the seeded finding was detected, so a
scanner crash or invocation error could masquerade as a successful control.

Auditor finding (Round-21 direction 3 -- CRITICAL): the Round-20 success path
printed the report's `Secret` field value, so a FAILED redaction could
disclose the very secret the control detected.

Auditor finding (Round-22 direction 2 -- CRITICAL): the Round-21 fix removed
the direct `Secret`/`Match` prints but STILL interpolated report-derived
metadata fields (`RuleID`, `File`, `StartLine`) into failure diagnostics and
negative-control output. A malicious or malformed report can place the
synthetic key in ANY report field, so NO report-derived value may ever reach
stdout/stderr. This version enforces:

  * every report element is TYPE-VALIDATED before it is examined: each row
    must be a JSON object whose `RuleID`/`File` are strings, `StartLine` an
    int (not bool), and `Secret`/`Match` strings when present. A malformed
    JSON array or malformed row value fails closed WITHOUT printing row
    contents (fixed messages: row index, field name, expected type);
  * diagnostics carry ONLY fixed messages, finding counts, argv-supplied
    expectations, and validated integers -- never report-derived strings;
  * the raw-key checks are NOT weakened: the assembled key must be absent
    from the report bytes and the captured scanner log even when the report
    fails schema validation (a malformed report that carries the key is a
    redaction failure, disclosed without printing the value);
  * requested-but-missing/unreadable capture logs fail closed;
  * a clean HEAD-tree (or clean-fixture) negative result requires the report
    to exist, parse as a JSON array, and contain ZERO findings -- scanner
    exit status alone never proves cleanliness;
  * unexpected internal errors are caught and reported with a FIXED message
    (no traceback, no interpolated values) and fail closed.

Modes:

  positive mode  -- the seeded synthetic key MUST be genuinely detected:
      - exit code must equal the reserved FINDINGS code (101; the scan is
        invoked with --exit-code 101 so "findings found" is distinct from
        execution errors, which exit 1);
      - the JSON report must exist and parse to a type-valid findings array;
      - the report must contain the expected finding for the seeded file with
        the EXPECTED RULE (aws-access-token) actually fired;
      - the report must contain no unexpected additional findings;
      - the assembled key value must be REDACTED: absent from the report
        bytes and from the captured scanner log (the key is passed via an
        environment variable name, never via argv or file content).

  negative mode  -- a clean repository/tree MUST pass:
      - exit code must be 0; a nonzero code is a scanner failure or a false
        positive, both of which fail the control;
      - the JSON report must exist, parse, be type-valid, and contain zero
        findings;
      - if --log-file is requested it must exist and be readable (fail
        closed), and if --key-env is provided the key must be absent from it.

Every failure path fails closed (exit 1); only a fully verified control exits 0.
No path in this module writes secret material or report-derived strings to
stdout or stderr.
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


def _validate_findings(data) -> tuple[list | None, str | None]:
    """Round-22: validate every element's type/schema BEFORE examination.

    Returns (data, None) when every row is a dict with the required field
    types, else (None, fixed-message error). Row VALUES are never included in
    the error: a malformed row must fail closed without printing its content.
    """
    if not isinstance(data, list):
        return None, "report is not a JSON array of findings"
    for i, row in enumerate(data):
        if not isinstance(row, dict):
            return None, (f"report row {i} is not a JSON object "
                          f"(malformed row; contents withheld)")
        rid = row.get("RuleID")
        if not isinstance(rid, str) or not rid:
            return None, (f"report row {i} field 'RuleID' is missing or not a "
                          f"non-empty string (contents withheld)")
        f = row.get("File")
        if not isinstance(f, str) or not f:
            return None, (f"report row {i} field 'File' is missing or not a "
                          f"non-empty string (contents withheld)")
        sl = row.get("StartLine")
        if isinstance(sl, bool) or not isinstance(sl, int) or sl < 0:
            return None, (f"report row {i} field 'StartLine' is not a "
                          f"non-negative integer (contents withheld)")
        for k in ("Secret", "Match"):
            v = row.get(k)
            if v is not None and not isinstance(v, str):
                return None, (f"report row {i} field {k!r} is neither absent "
                              f"nor a string (contents withheld)")
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


def _check_redaction(report_path: str, key: str, log_path: str | None,
                     fails: list[str], expected_rows: list[dict]) -> None:
    """The assembled key must be absent from the report bytes, the validated
    finding fields, and the captured log. Only booleans and byte counts are
    printed -- never values (Round-21/22 critical rules)."""
    raw = b""
    if report_path and os.path.exists(report_path):
        try:
            raw = open(report_path, "rb").read()
        except OSError:
            fails.append("the scan report is unreadable as bytes (fail closed)")
    key_bytes = key.encode()
    report_redacted = key_bytes not in raw
    field_redacted = all(r.get("Secret") != key and r.get("Match") != key
                         for r in expected_rows)
    if not report_redacted:
        fails.append("the raw seeded key appears in the scan report bytes "
                     "(redaction FAILED)")
    if not field_redacted:
        fails.append("a finding's Secret or Match field carries the raw key "
                     "(redaction FAILED)")
    if not raw:
        if report_path and os.path.exists(report_path):
            fails.append("the scan report is unreadable as bytes (fail closed)")
        # a missing report already fails elsewhere; do not claim verification
    elif report_redacted and field_redacted:
        print("redaction verified: the raw key is absent from all %d report "
              "bytes and from every finding's Secret/Match field" % len(raw))
    else:
        # Round-21 3.1/3.2 + Round-22: disclose ONLY the boolean outcomes;
        # never the values, and never any other report-derived field.
        print("redaction FAILED: report_bytes_clean=%s fields_clean=%s "
              "(values are never printed)" % (report_redacted, field_redacted))
    log = _read_log_bytes(log_path, fails) if log_path else None
    if log is not None:
        if key_bytes in log:
            fails.append("the raw seeded key appears in the captured scanner log "
                         "(redaction FAILED)")
        else:
            print("capture log verified: the raw key is absent from all %d log "
                  "bytes" % len(log))


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
        data = None
    else:
        data, verr = _validate_findings(data)
        if verr:
            fails.append(verr)
            data = None
    validated = data or []
    expected = [r for r in validated
                if r.get("RuleID") == expect_rule and r.get("File") == expect_file]
    others = [r for r in validated if r not in expected]
    if not expected:
        fails.append(f"no finding for rule {expect_rule!r} in file {expect_file!r} "
                     f"-- the expected rule did not fire")
    else:
        # Safe metadata ONLY: argv-supplied expectations, the count, and
        # schema-validated integer line numbers (Round-22: StartLine values
        # are type-validated ints before printing; no report-derived strings).
        print("expected finding verified: rule=%s file=%s count=%d line(s)=%s"
              % (expect_rule, expect_file, len(expected),
                 ",".join(str(r["StartLine"]) for r in expected)))
    if others:
        # Round-22: report ONLY the count of unexpected findings. RuleID/File/
        # StartLine values are report-derived strings and must never be
        # interpolated into diagnostics (the Round-22 audit found exactly this
        # leak path).
        fails.append("unexpected additional finding(s): %d (rule/file/line values "
                     "are never printed)" % len(others))
    key = os.environ.get(key_env, "")
    if not key:
        fails.append(f"the seeded key is not available via ${key_env} "
                     f"(redaction cannot be proven -- fail closed)")
    else:
        # redaction is checked against the report bytes even when the schema
        # validation failed: a malformed report carrying the key is still a
        # disclosure, and must be reported without printing the value.
        _check_redaction(report_path, key, log_path, fails, expected)
    if fails:
        for f in fails:
            print("FAIL: " + f)
        return 1
    print("PASS: positive control -- seeded synthetic key genuinely detected by "
          "rule %s in %s, report well-formed and type-valid, key redacted "
          "everywhere it was checked (values never printed)"
          % (expect_rule, expect_file))
    return 0


def verify_negative(report_path: str, exit_code: int,
                    log_path: str | None = None, key_env: str = "SEEDED_KEY") -> int:
    fails: list[str] = []
    if exit_code != 0:
        fails.append(f"clean scan exited {exit_code} (expected 0): "
                     f"scanner failure or false positive -- control fails")
    # Round-21 3.4 + Round-22: exit 0 is NOT sufficient; the report must
    # exist, parse, be type-valid, and contain zero findings.
    data, err = _load_report(report_path)
    if err:
        fails.append(err)
    else:
        data, verr = _validate_findings(data)
        if verr:
            fails.append(verr)
        elif data:
            # Round-22: report ONLY the count. RuleID/File values are
            # report-derived and must never be printed.
            fails.append("clean scan produced %d finding(s) (rule/file values "
                         "are never printed)" % len(data))
        else:
            print("negative report verified: exists, parses, type-valid, zero "
                  "findings (exit status alone was not trusted)")
    key = os.environ.get(key_env, "")
    if log_path:
        log = _read_log_bytes(log_path, fails)
        if log is not None:
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
    # Round-22: unexpected internal errors must fail closed with a FIXED
    # message -- a traceback could interpolate unexpected values, and the
    # key-absence guarantee must hold on every path including exceptions.
    try:
        if a.mode == "positive":
            return verify_positive(a.report, a.exit_code, a.expect_rule,
                                   a.expect_file, a.key_env, a.log_file)
        return verify_negative(a.report, a.exit_code, a.log_file, a.key_env)
    except Exception as e:  # noqa: BLE001 -- deliberately broad, fail closed
        print("FAIL: internal error (%s) -- the control fails closed; no "
              "report-derived or secret values are printed"
              % type(e).__name__)
        return 1


if __name__ == "__main__":
    sys.exit(main())
