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


def check_pushlog_schema(path: str) -> tuple[list[str], list[str]]:
    """Round-19 hardening (auditor direction 4). Returns (findings, pending).

    Enforced on `## <ts> | <ENTRY-ID> | <verb> | <subject>` entries:
      - ENTRY-IDs must be unique (duplicate IDs across competing histories are
        exactly the P-015-CONFIRM collision class) and match P-<slug>.
      - push-confirm (CONFIRM) entries must carry branch=, remote_sha=<full
        40-hex> (never a placeholder), credential=, and evidence= (malformed
        confirmations / missing required fields fail).
      - push-intent entries must carry branch= and credential=.
      - push (legacy single-phase) entries must carry branch= and credential=.
      - Pairing: every CONFIRM must reconcile a preceding unconfirmed INTENT
        (same branch, FIFO) or a legacy push entry whose sha was a disclosed
        placeholder (the P-011 correction case). A CONFIRM with nothing to
        confirm fails.
      - PENDING (reported, not a failure): an INTENT whose push has no CONFIRM
        yet, and a legacy placeholder-sha push with no CONFIRM — the two-phase
        protocol is incomplete and must not be implied complete (direction 4).
    """
    findings: list[str] = []
    pending: list[str] = []
    entries: list[tuple[int, str, str, str]] = []
    for n, raw in enumerate(open(path, encoding="utf-8"), 1):
        line = raw.rstrip("\n")
        if not line.startswith("## "):
            continue
        parts = [p.strip() for p in line[3:].split("|")]
        if len(parts) < 4:
            findings.append(f"FAIL: {path}:{n}: malformed log entry (grammar mismatch)")
            continue
        ts, eid, verb, subject = parts[0], parts[1], parts[2], parts[3]
        entries.append((n, eid, verb, subject))
        if not re.fullmatch(r"P-[A-Za-z0-9][A-Za-z0-9._-]*", eid):
            findings.append(f"FAIL: {path}:{n}: ENTRY-ID {eid!r} does not match the "
                            f"push-log grammar P-<slug>")
        if verb == "push-confirm":
            for key in ("branch=", "remote_sha=", "credential=", "evidence="):
                if key not in subject:
                    findings.append(f"FAIL: {path}:{n}: CONFIRM entry is missing the "
                                    f"required field {key[:-1]!r} (malformed confirmation)")
            m = re.search(r"remote_sha=([0-9a-f]{40})(?:\s|$)", subject)
            if not m:
                sm = re.search(r"remote_sha=(\S+)", subject)
                got = sm.group(1).rstrip('",') if sm else None
                if got is None or got.lower() in PLACEHOLDER_SHAS or not re.fullmatch(
                        r"[0-9a-f]{40}", got):
                    findings.append(
                        f"FAIL: {path}:{n}: CONFIRM entry carries {got!r} instead of the real "
                        f"40-hex remote SHA (push-log truthfulness violated)")
        elif verb in ("push-intent", "push"):
            for key in ("branch=", "credential="):
                if key not in subject:
                    findings.append(f"FAIL: {path}:{n}: {verb} entry is missing the "
                                    f"required field {key[:-1]!r}")
    # duplicate ENTRY-IDs (first occurrence wins; also enforced by check_sequence,
    # re-checked here so --pushlog-schema-only runs catch it)
    seen: dict[str, int] = {}
    for n, eid, verb, subject in entries:
        if eid in seen:
            findings.append(f"FAIL: {path}:{n}: duplicate ENTRY-ID {eid} "
                            f"(first at line {seen[eid]}; one ID, one meaning)")
        else:
            seen[eid] = n
    # INTENT/CONFIRM pairing, FIFO per branch; CONFIRM may also correct a legacy
    # placeholder-sha push entry (the documented P-011 case)
    open_intents: dict[str, list[tuple[int, str]]] = {}
    legacy_placeholder: dict[str, list[tuple[int, str]]] = {}
    for n, eid, verb, subject in entries:
        bm = re.search(r"branch=(\S+)", subject)
        branch = bm.group(1) if bm else None
        if verb == "push-intent":
            if branch:
                open_intents.setdefault(branch, []).append((n, eid))
        elif verb == "push":
            sm = re.search(r"sha=(\S+)", subject)
            if branch and sm and sm.group(1).lower() in PLACEHOLDER_SHAS:
                legacy_placeholder.setdefault(branch, []).append((n, eid))
        elif verb == "push-confirm":
            if not branch:
                continue
            if open_intents.get(branch):
                open_intents[branch].pop(0)
            elif legacy_placeholder.get(branch):
                legacy_placeholder[branch].pop(0)
            else:
                findings.append(f"FAIL: {path}:{n}: CONFIRM {eid} has no preceding "
                                f"unconfirmed INTENT (or placeholder push) for branch {branch!r}")
    for branch, queue in open_intents.items():
        for n, eid in queue:
            pending.append(f"PENDING: {path}:{n}: {eid} for branch {branch} has no CONFIRM "
                           f"entry (two-phase protocol incomplete; do not imply completion)")
    for branch, queue in legacy_placeholder.items():
        for n, eid in queue:
            pending.append(f"PENDING: {path}:{n}: legacy push {eid} for branch {branch} still "
                           f"carries a placeholder sha with no CONFIRM correction")
    return findings, pending


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("logfile")
    ap.add_argument("--commit-time", help="ISO-8601 commit time carrying this log version")
    ap.add_argument("--max-drift", type=int, default=300, help="seconds of allowed drift (default 300)")
    ap.add_argument("--parent-prefix", metavar="OLDER_FILE",
                    help="require OLDER_FILE to be a byte-prefix of logfile (append-only)")
    ap.add_argument("--pushlog-schema", action="store_true",
                    help="also enforce push-log schema: unique IDs, well-formed INTENT/"
                         "CONFIRM entries with required fields, and INTENT/CONFIRM pairing "
                        "(unpaired INTENTs are reported as PENDING, not failures)")
    a = ap.parse_args()
    findings = check_sequence(a.logfile, a.commit_time, a.max_drift)
    if a.parent_prefix:
        findings += check_parent_prefix(a.logfile, a.parent_prefix)
    pending: list[str] = []
    if a.pushlog_schema:
        schema_findings, pending = check_pushlog_schema(a.logfile)
        findings += schema_findings
    for p in pending:
        print(p)
    print("\n".join(findings) if findings else f"PASS: {a.logfile} "
          f"{'sequence + append-only ' if a.parent_prefix else ''}clean")
    return 0 if not findings else 1


if __name__ == "__main__":
    sys.exit(main())
