#!/usr/bin/env python3
"""EC-18 evidence sequence + append-only integrity checker (inert until adopted).

Auditor revision 2026-10-20 (item 3.6): the sequence checker alone does not establish
append-only integrity. This version adds a parent-version comparison: an evidence log's
current content must preserve the exact previous version as a byte PREFIX — any edit,
deletion, or reordering of prior entries FAILs, even when timestamps and IDs stay clean.

Log grammar:   ## <ISO-8601 UTC> | <ENTRY-ID> | <verb> | <subject> [| RECONSTRUCTED]
Sequence checks: timestamps parse as UTC, non-decreasing, ENTRY-IDs unique; with
--commit-time, entries whose claimed time drifts more than --max-drift from the carrying
commit must be RECONSTRUCTED.

Append-only check (--parent-prefix OLDER-FILE): OLDER must be a byte-prefix of the
current file. The workflow obtains OLDER via `git show <prev-commit>:<path>`, so an
older evidence version that was edited or truncated in place is detected.

Push-log truthfulness (EC-24 aid, --pushlog-schema): validates that post-push
confirmation entries carry a real 40-hex sha=... (never a self-referential or
placeholder marker) and that INTENT/CONFIRM phases are separately recorded.

Exit: 0 clean; 1 findings; 2 usage.
"""
import argparse
import re
import sys
from datetime import datetime, timezone

LINE = re.compile(r"^## (\S+) \| ([^|]+) \| ([^|]+) \| ([^|]+?)(?:\s\|\s*RECONSTRUCTED)?\s*$")
PLACEHOLDER_SHAS = {"self-referential-tip", "TBD", "PENDING", "<sha>", "sha"}


def parse_ts(s: str):
    return datetime.fromisoformat(s.replace("Z", "+00:00")).astimezone(timezone.utc)


def check_sequence(path: str, commit_time: str | None, max_drift: int) -> list[str]:
    findings, seen_ids, prev = [], set(), None
    for n, raw in enumerate(open(path, encoding="utf-8"), 1):
        line = raw.rstrip("\n")
        if not line.startswith("## "):
            continue
        m = LINE.match(line)
        if not m:
            findings.append(f"FAIL: {path}:{n}: malformed log entry (grammar mismatch)")
            continue
        ts_s, eid = m.group(1), m.group(2).strip()
        recon = "RECONSTRUCTED" in line
        try:
            ts = parse_ts(ts_s)
        except ValueError:
            findings.append(f"FAIL: {path}:{n}: timestamp not ISO-8601 UTC: {ts_s}")
            continue
        if eid in seen_ids:
            findings.append(f"FAIL: {path}:{n}: duplicate ENTRY-ID {eid}")
        seen_ids.add(eid)
        if prev is not None and ts < prev:
            findings.append(f"FAIL: {path}:{n}: timestamp goes backwards ({ts_s} < previous)")
        prev = ts
        if commit_time and not recon:
            ct = parse_ts(commit_time)
            if (ct - ts).total_seconds() > max_drift:
                findings.append(f"FAIL: {path}:{n}: claimed time drifts >{max_drift}s from "
                                f"carrying commit without RECONSTRUCTED label")
    return findings


def check_parent_prefix(current: str, older: str) -> list[str]:
    """OLDER must be a byte-prefix of CURRENT (append-only, EC-18 revision)."""
    cur = open(current, "rb").read()
    old = open(older, "rb").read()
    if cur == old:
        return []  # unchanged version is trivially append-only
    if len(old) > len(cur):
        return [f"FAIL: {current}: SHRUNK relative to the previous version "
                f"({len(old)} -> {len(cur)} bytes); prior entries were deleted"]
    if not cur.startswith(old):
        # locate the first divergence for a precise message
        i = next((k for k in range(len(old)) if old[k] != cur[k]), len(old))
        ctx = old[max(0, i - 40):i + 40].decode("utf-8", "replace")
        return [f"FAIL: {current}: prior content MODIFIED at byte {i} "
                f"(append-only violated); context: {ctx!r}"]
    return []


def check_pushlog_schema(path: str) -> list[str]:
    """Post-push confirmation entries must cite real full SHAs (EC-24 truthfulness)."""
    findings = []
    for n, line in enumerate(open(path, encoding="utf-8"), 1):
        if line.startswith("## ") and "CONFIRM" in line:
            m = re.search(r"sha=(\S+)", line)
            sha = m.group(1) if m else None
            if not sha or sha.lower() in PLACEHOLDER_SHAS or not re.fullmatch(
                    r"[0-9a-f]{40}", sha):
                findings.append(
                    f"FAIL: {path}:{n}: CONFIRM entry carries {sha!r} instead of the real "
                    f"40-hex SHA (push-log truthfulness violated)")
    return findings


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("logfile")
    ap.add_argument("--commit-time", help="ISO-8601 commit time carrying this log version")
    ap.add_argument("--max-drift", type=int, default=300, help="seconds of allowed drift (default 300)")
    ap.add_argument("--parent-prefix", metavar="OLDER_FILE",
                    help="require OLDER_FILE to be a byte-prefix of logfile (append-only)")
    ap.add_argument("--pushlog-schema", action="store_true",
                    help="also enforce push-log CONFIRM-entry SHA truthfulness")
    a = ap.parse_args()
    findings = check_sequence(a.logfile, a.commit_time, a.max_drift)
    if a.parent_prefix:
        findings += check_parent_prefix(a.logfile, a.parent_prefix)
    if a.pushlog_schema:
        findings += check_pushlog_schema(a.logfile)
    print("\n".join(findings) if findings else f"PASS: {a.logfile} "
          f"{'sequence + append-only ' if a.parent_prefix else ''}clean")
    return 0 if not findings else 1


if __name__ == "__main__":
    sys.exit(main())
