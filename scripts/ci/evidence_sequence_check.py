#!/usr/bin/env python3
"""EC-18/EC-24 evidence sequence + append-only + remote-ref verification checker.

Auditor revision 2026-10-20 (item 3.6): the sequence checker alone does not
establish append-only integrity. This version adds a parent-version comparison:
an evidence log's current content must preserve the exact previous version as a
byte PREFIX — any edit, deletion, or reordering of prior entries FAILs, even
when timestamps and IDs stay clean.

Round-20 hardening (auditor direction 4): the checker validated confirmation
field SYNTAX (40-hex remote_sha, required keys) but never proved that the
recorded remote_sha equals the live remote ref for the stated branch, nor that
the evidence= field actually supports the recorded SHA. New --remote-verify
layer:

  - every CONFIRM/RECONCILE entry is checked against the LIVE remote ref
    (git ls-remote; fetch failures and absent refs fail closed);
  - live == remote_sha  => classification CURRENT;
    live != remote_sha  => classification HISTORICAL (valid for the push it
    describes; it does NOT confirm the newer tip — that state is disclosed,
    never smoothed over);
  - the entry's evidence must actually support the recorded SHA: the SHA and
    the claimed branch must appear in the evidence field (or the entry body),
    and any cited Actions run ID whose manifest record shows a DIFFERENT
    head_sha fails (fabricated support); run IDs outside the manifest window
    are disclosed as unverified, not trusted;
  - the branch under test (--branch/--expect-tip) must have its CURRENT tip
    declared: an unpaired INTENT (PENDING, disclosed) or a CURRENT
    CONFIRM/RECONCILE. A tip whose push was never declared FAILS (the
    Round-19 follow-up-push class: pushes that carried confirmations for
    earlier tips without declaring themselves);
  - RECONCILE is a new post-hoc verb for reconciling pushes that already
    happened without a preceding INTENT: it records the push truthfully with
    dated live-ref evidence, never backdates, and cannot pair with or consume
    INTENT records.

Log grammar:   ## <ISO-8601 UTC> | <ENTRY-ID> | <verb> | <subject> [| RECONSTRUCTED]
Sequence checks: timestamps parse as UTC, non-decreasing, ENTRY-IDs unique; with
--commit-time, entries whose claimed time drifts more than --max-drift from the
carrying commit must be RECONSTRUCTED.

Exit: 0 clean; 1 findings; 2 usage.
"""
import argparse
import re
import subprocess
import sys
from datetime import datetime, timezone

LINE = re.compile(r"^## (\S+) \| ([^|]+) \| ([^|]+) \| ([^|]+?)(?:\s\|\s*RECONSTRUCTED)?\s*$")
PLACEHOLDER_SHAS = {"self-referential-tip", "TBD", "PENDING", "<sha>", "sha"}
SHA40 = re.compile(r"\b([0-9a-f]{40})\b")
RUN_ID = re.compile(r"\brun(?:_id)?\s*[=:]?\s*(\d{6,})\b", re.I)


def parse_ts(s: str):
    return datetime.fromisoformat(s.replace("Z", "+00:00")).astimezone(timezone.utc)


def parse_entries(path: str):
    """Return entry dicts {n, ts, eid, verb, subject, body, raw} for `## ` headers."""
    entries, cur = [], None
    for n, raw in enumerate(open(path, encoding="utf-8"), 1):
        line = raw.rstrip("\n")
        if line.startswith("## "):
            if cur:
                entries.append(cur)
            m = LINE.match(line)
            cur = {"n": n, "raw": line, "m": m, "body": []}
            if m:
                cur["ts"], cur["eid"], cur["verb"], cur["subject"] = (
                    m.group(1), m.group(2).strip(), m.group(3).strip(), m.group(4).strip())
        elif cur:
            cur["body"].append(line)
    if cur:
        entries.append(cur)
    return entries


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


def check_pushlog_schema(path: str) -> tuple[list[str], list[str], dict]:
    """Round-19/20 hardening. Returns (findings, pending, state).

    state carries the parsed pairing state for the remote-verification layer:
    open_intents (branch -> [(line, eid)]) after FIFO pairing, and the entry
    list.
    """
    findings: list[str] = []
    pending: list[str] = []
    entries = parse_entries(path)
    open_intents: dict[str, list[tuple[int, str]]] = {}
    legacy_placeholder: dict[str, list[tuple[int, str]]] = {}
    for e in entries:
        if not e["m"]:
            findings.append(f"FAIL: {path}:{e['n']}: malformed log entry (grammar mismatch)")
            continue
        eid, verb, subject = e["eid"], e["verb"], e["subject"]
        if not re.fullmatch(r"P-[A-Za-z0-9][A-Za-z0-9._-]*", eid):
            findings.append(f"FAIL: {path}:{e['n']}: ENTRY-ID {eid!r} does not match the "
                            f"push-log grammar P-<slug>")
        if verb == "push-confirm":
            for key in ("branch=", "remote_sha=", "credential=", "evidence="):
                if key not in subject:
                    findings.append(f"FAIL: {path}:{e['n']}: CONFIRM entry is missing the "
                                    f"required field {key[:-1]!r} (malformed confirmation)")
            m = re.search(r"remote_sha=([0-9a-f]{40})(?:\s|$)", subject)
            if not m:
                sm = re.search(r"remote_sha=(\S+)", subject)
                got = sm.group(1).rstrip('",') if sm else None
                if got is None or got.lower() in PLACEHOLDER_SHAS or not re.fullmatch(
                        r"[0-9a-f]{40}", got):
                    findings.append(
                        f"FAIL: {path}:{e['n']}: CONFIRM entry carries {got!r} instead of the real "
                        f"40-hex remote SHA (push-log truthfulness violated)")
        elif verb == "push-reconcile":
            # Round-20 verb: post-hoc reconciliation of an undeclared past push.
            for key in ("branch=", "remote_sha=", "credential=", "evidence="):
                if key not in subject:
                    findings.append(f"FAIL: {path}:{e['n']}: RECONCILE entry is missing the "
                                    f"required field {key[:-1]!r}")
            m = re.search(r"remote_sha=([0-9a-f]{40})(?:\s|$)", subject)
            if not m:
                findings.append(f"FAIL: {path}:{e['n']}: RECONCILE entry remote_sha must be a "
                                f"real 40-hex SHA (never a placeholder)")
        elif verb in ("push-intent", "push"):
            for key in ("branch=", "credential="):
                if key not in subject:
                    findings.append(f"FAIL: {path}:{e['n']}: {verb} entry is missing the "
                                    f"required field {key[:-1]!r}")
    # duplicate ENTRY-IDs (first occurrence wins; also enforced by check_sequence,
    # re-checked here so schema-only runs catch it)
    seen: dict[str, int] = {}
    for e in entries:
        if not e["m"]:
            continue
        eid = e["eid"]
        if eid in seen:
            findings.append(f"FAIL: {path}:{e['n']}: duplicate ENTRY-ID {eid} "
                            f"(first at line {seen[eid]}; one ID, one meaning)")
        else:
            seen[eid] = e["n"]
    # INTENT/CONFIRM pairing, FIFO per branch; CONFIRM may also correct a legacy
    # placeholder-sha push entry (the documented P-011 case). RECONCILE records
    # are post-hoc by design and NEVER consume INTENT records.
    for e in entries:
        if not e["m"]:
            continue
        verb, subject = e["verb"], e["subject"]
        bm = re.search(r"branch=(\S+)", subject)
        branch = bm.group(1) if bm else None
        if verb == "push-intent":
            if branch:
                open_intents.setdefault(branch, []).append((e["n"], e["eid"]))
        elif verb == "push":
            sm = re.search(r"sha=(\S+)", subject)
            if branch and sm and sm.group(1).lower() in PLACEHOLDER_SHAS:
                legacy_placeholder.setdefault(branch, []).append((e["n"], e["eid"]))
        elif verb == "push-confirm":
            if not branch:
                continue
            if open_intents.get(branch):
                open_intents[branch].pop(0)
            elif legacy_placeholder.get(branch):
                legacy_placeholder[branch].pop(0)
            else:
                findings.append(f"FAIL: {path}:{e['n']}: CONFIRM {e['eid']} has no preceding "
                                f"unconfirmed INTENT (or placeholder push) for branch {branch!r}")
    for branch, queue in open_intents.items():
        for n, eid in queue:
            pending.append(f"PENDING: {path}:{n}: {eid} for branch {branch} has no CONFIRM "
                           f"entry (two-phase protocol incomplete; do not imply completion)")
    for branch, queue in legacy_placeholder.items():
        for n, eid in queue:
            pending.append(f"PENDING: {path}:{n}: legacy push {eid} for branch {branch} still "
                           f"carries a placeholder sha with no CONFIRM correction")
    return findings, pending, {"entries": entries, "open_intents": open_intents}


def ls_remote_tip(remote: str, branch: str) -> tuple[str | None, str | None]:
    """Return (sha, error): the live tip of refs/heads/<branch> on <remote>."""
    ref = f"refs/heads/{branch}"
    try:
        r = subprocess.run(["git", "ls-remote", remote, ref],
                           capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired) as e:
        return None, f"ls-remote failed for {branch!r} ({e})"
    if r.returncode != 0:
        return None, (f"ls-remote exited {r.returncode} for {branch!r}: "
                      f"{r.stderr.strip()[:120]}")
    for line in r.stdout.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[1] == ref and re.fullmatch(r"[0-9a-f]{40}", parts[0]):
            return parts[0], None
    return None, f"remote ref {ref} not found on {remote} (missing remote verification)"


def check_remote_verification(path: str, state: dict, remote: str,
                              branch_under_test: str | None, expect_tip: str | None,
                              manifest: dict | None) -> tuple[list[str], list[str]]:
    """Round-20 layer: live-ref comparison + evidence support + tip declaration."""
    findings: list[str] = []
    info: list[str] = []
    entries = state["entries"]
    open_intents: dict[str, list[tuple[int, str]]] = state["open_intents"]

    live_cache: dict[str, tuple[str | None, str | None]] = {}

    def live(branch: str):
        if branch not in live_cache:
            live_cache[branch] = ls_remote_tip(remote, branch)
        return live_cache[branch]

    # branch -> sha of its latest CURRENT-classified confirm/reconcile
    current_for: dict[str, str] = {}
    for e in entries:
        if not e["m"] or e["verb"] not in ("push-confirm", "push-reconcile"):
            continue
        subject, eid, verb = e["subject"], e["eid"], e["verb"]
        bm = re.search(r"branch=(\S+)", subject)
        sm = re.search(r"remote_sha=([0-9a-f]{40})", subject)
        em = re.search(r'evidence="([^"]*)"', subject)
        if not (bm and sm):
            continue  # field grammar already failed in schema check
        branch, sha = bm.group(1), sm.group(1)
        evidence = em.group(1) if em else ""
        body = "\n".join(e["body"])
        # 1. evidence must support the recorded SHA (field first, body second)
        if sha not in evidence:
            if sha in body:
                info.append(f"NOTE: {eid}: evidence field does not inline the SHA "
                            f"(body citation only; legacy form)")
            else:
                findings.append(f"FAIL: {path}:{e['n']}: {eid} records remote_sha={sha[:12]}… "
                                f"but neither the evidence field nor the entry body contains "
                                f"that SHA (fabricated or unsupported evidence)")
        # 2. claimed branch must match the evidence
        if branch not in evidence and branch not in body:
            findings.append(f"FAIL: {path}:{e['n']}: {eid} claims branch {branch!r} but the "
                            f"evidence does not mention that branch (mismatched branch)")
        # 3. cited run IDs: positive contradiction fails; absent-from-manifest disclosed
        for rid in RUN_ID.findall(evidence):
            if manifest is None:
                continue
            rec = manifest.get("runs", {}).get(rid)
            if rec is None:
                info.append(f"NOTE: {eid}: cited run {rid} is outside the runs-manifest "
                            f"window (unverified, not trusted)")
            elif rec.get("head_sha") != sha:
                findings.append(f"FAIL: {path}:{e['n']}: {eid} cites run {rid} as evidence for "
                                f"{sha[:12]}… but that run's head_sha is "
                                f"{str(rec.get('head_sha'))[:12]}… (evidence does not support "
                                f"the recorded SHA)")
        # 4. live-ref comparison
        tip, err = live(branch)
        if err:
            findings.append(f"FAIL: {path}:{e['n']}: {eid} could not be verified against the "
                            f"live remote ref: {err} (fail closed)")
            continue
        if tip == sha:
            current_for[branch] = sha
            info.append(f"CURRENT: {eid} confirms the live tip of {branch} ({sha[:12]}…)")
        else:
            info.append(f"HISTORICAL: {eid} confirmed {branch} @ {sha[:12]}…; the live tip is "
                        f"{tip[:12]}… — the confirmation is valid for the push it describes "
                        f"and does NOT confirm the newer tip")

    # 5. branch under test: candidate must be the live tip and the tip must be declared
    if branch_under_test:
        tip, err = live(branch_under_test)
        if err:
            findings.append(f"FAIL: branch under test {branch_under_test!r} could not be "
                            f"verified: {err} (fail closed)")
        elif expect_tip and tip != expect_tip:
            findings.append(f"FAIL: live tip of {branch_under_test} is {tip[:12]}…, not the "
                            f"expected candidate {expect_tip[:12]}… — the checked-out "
                            f"candidate is not the live ref (stale or raced)")
        else:
            declared = current_for.get(branch_under_test)
            unpaired = open_intents.get(branch_under_test, [])
            if declared:
                info.append(f"DECLARED: {branch_under_test} current tip is CONFIRMED-CURRENT")
            elif unpaired:
                info.append(f"PENDING: {branch_under_test} current tip {tip[:12]}… is declared "
                            f"by unpaired INTENT {unpaired[0][1]} (line {unpaired[0][0]}); its "
                            f"CONFIRM is delivered by the next legitimate push (disclosed "
                            f"limitation, never backdated)")
            else:
                findings.append(f"FAIL: the live tip of {branch_under_test} "
                                f"({tip[:12]}…) is UNDECLARED — no unpaired INTENT declares "
                                f"its push and no CONFIRM/RECONCILE matches it (the "
                                f"Round-19 follow-up-push class: pushes that carried "
                                f"confirmations for earlier tips without declaring "
                                f"themselves; reconcile truthfully, never backdate)")

    # 6. cross-branch disclosure: branches whose tips moved past their last confirmation
    for e in entries:
        if not e["m"] or e["verb"] != "push":
            continue
        bm = re.search(r"branch=(\S+)", e["subject"])
        sm = re.search(r"sha=([0-9a-f]{40})", e["subject"])
        if not (bm and sm):
            continue
        branch, sha = bm.group(1), sm.group(1)
        tip, err = live(branch)
        if err:
            continue  # legacy entries: disclosed in the confirm layer if applicable
        if tip != sha and current_for.get(branch) is None and not open_intents.get(branch):
            info.append(f"NOTE: legacy push of {branch} recorded {sha[:12]}…; live tip "
                        f"{tip[:12]}… is beyond the log's last confirmation for that branch "
                        f"(disclosed; reconcile by a reviewed record if this branch re-enters "
                        f"the protocol)")
    return findings, info


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("logfile")
    ap.add_argument("--commit-time", help="ISO-8601 commit time carrying this log version")
    ap.add_argument("--max-drift", type=int, default=300, help="seconds of allowed drift (default 300)")
    ap.add_argument("--parent-prefix", metavar="OLDER_FILE",
                    help="require OLDER_FILE to be a byte-prefix of logfile (append-only)")
    ap.add_argument("--pushlog-schema", action="store_true",
                    help="also enforce push-log schema: unique IDs, well-formed INTENT/"
                         "CONFIRM/RECONCILE entries with required fields, and INTENT/CONFIRM "
                         "pairing (unpaired INTENTs are reported as PENDING, not failures)")
    ap.add_argument("--remote-verify", action="store_true",
                    help="Round-20 layer: verify every CONFIRM/RECONCILE against the live "
                         "remote ref, prove the evidence supports the recorded SHA, and "
                         "require the branch under test to declare its current tip")
    ap.add_argument("--remote", default="origin",
                    help="remote (URL or path) for ls-remote verification (default: origin)")
    ap.add_argument("--branch", default=None,
                    help="the branch under test (required with --expect-tip)")
    ap.add_argument("--expect-tip", default=None, metavar="SHA",
                    help="the expected live tip of --branch (in CI: the exact pushed SHA)")
    ap.add_argument("--runs-manifest", default=None, metavar="FILE",
                    help="JSON {\"runs\": {id: {head_sha, conclusion}}} for evidence run-ID "
                         "verification (cited runs with a contradicting head_sha FAIL)")
    a = ap.parse_args()

    findings = check_sequence(a.logfile, a.commit_time, a.max_drift)
    if a.parent_prefix:
        findings += check_parent_prefix(a.logfile, a.parent_prefix)
    pending: list[str] = []
    info: list[str] = []
    if a.pushlog_schema or a.remote_verify:
        schema_findings, pending, state = check_pushlog_schema(a.logfile)
        findings += schema_findings
        if a.remote_verify:
            manifest = None
            if a.runs_manifest:
                try:
                    import json
                    manifest = json.load(open(a.runs_manifest, encoding="utf-8"))
                except (OSError, ValueError) as e:
                    print(f"FAIL: runs manifest unusable ({e}); fail closed")
                    return 1
            remote = a.remote
            if remote == "origin":
                r = subprocess.run(["git", "remote", "get-url", "origin"],
                                   capture_output=True, text=True)
                if r.returncode == 0 and r.stdout.strip():
                    remote = r.stdout.strip()
            rv_findings, rv_info = check_remote_verification(
                a.logfile, state, remote, a.branch, a.expect_tip, manifest)
            findings += rv_findings
            info += rv_info
    for p in pending:
        print(p)
    for i in info:
        print(i)
    print("\n".join(findings) if findings else f"PASS: {a.logfile} "
          f"{'sequence + append-only ' if a.parent_prefix else ''}"
          f"{'schema ' if (a.pushlog_schema or a.remote_verify) else ''}"
          f"{'remote-verified ' if a.remote_verify else ''}clean")
    return 0 if not findings else 1


if __name__ == "__main__":
    sys.exit(main())
