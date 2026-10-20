#!/usr/bin/env python3
"""EC-18 evidence sequence checker (inert until adopted).

Log grammar:  ## <ISO-8601 UTC> | <ENTRY-ID> | <verb> | <subject> [| RECONSTRUCTED]
Checks: timestamps parse as UTC, non-decreasing, ENTRY-IDs unique; with --commit-time, entries
whose claimed time drifts more than --max-drift from the carrying commit must be RECONSTRUCTED.
"""
import argparse, re, sys
from datetime import datetime, timezone

LINE = re.compile(r"^## (\S+) \| ([^|]+) \| ([^|]+) \| ([^|]+?)(?:\s\|\s*RECONSTRUCTED)?\s*$")

def parse_ts(s: str):
    return datetime.fromisoformat(s.replace("Z", "+00:00")).astimezone(timezone.utc)

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("logfile")
    ap.add_argument("--commit-time", help="ISO-8601 commit time carrying this log version")
    ap.add_argument("--max-drift", type=int, default=300, help="seconds of allowed drift (default 300)")
    a = ap.parse_args()
    findings, seen_ids, prev = [], set(), None
    for n, raw in enumerate(open(a.logfile, encoding="utf-8"), 1):
        line = raw.rstrip("\n")
        if not line.startswith("## "):
            continue
        m = LINE.match(line)
        if not m:
            findings.append(f"FAIL: {a.logfile}:{n}: malformed log entry (grammar mismatch)")
            continue
        ts_s, eid, verb, subj = m.group(1), m.group(2).strip(), m.group(3).strip(), m.group(4).strip()
        recon = "RECONSTRUCTED" in line
        try:
            ts = parse_ts(ts_s)
        except ValueError:
            findings.append(f"FAIL: {a.logfile}:{n}: timestamp not ISO-8601 UTC: {ts_s}")
            continue
        if eid in seen_ids:
            findings.append(f"FAIL: {a.logfile}:{n}: duplicate ENTRY-ID {eid}")
        seen_ids.add(eid)
        if prev is not None and ts < prev:
            findings.append(f"FAIL: {a.logfile}:{n}: timestamp goes backwards ({ts_s} < previous)")
        prev = ts
        if a.commit_time and not recon:
            ct = parse_ts(a.commit_time)
            if (ct - ts).total_seconds() > a.max_drift:
                findings.append(f"FAIL: {a.logfile}:{n}: claimed time drifts >{a.max_drift}s from "
                                f"carrying commit without RECONSTRUCTED label")
    print("\n".join(findings) if findings else f"PASS: {a.logfile} sequence clean")
    return 0 if not findings else 1

if __name__ == "__main__":
    sys.exit(main())
