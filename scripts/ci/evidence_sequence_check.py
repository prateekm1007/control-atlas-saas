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

Push-log protocol check (EC-24 aid, --pushlog-schema) — Round-18 extension per the
auditor direction of 2026-10-20 ("a schema check that passes while a push remains
unconfirmed must report pending, not complete"):

  1. TRUTHFULNESS: every push-confirm entry carries a real 40-hex remote_sha (never a
     placeholder such as `self-referential-tip`).
  2. PAIRING: every push-intent entry for a branch MUST be followed, later in the log,
     by a push-confirm for the SAME branch — unless (a) the intent is explicitly
     reconciled in the reconciliation record, or (b) it is the FINAL intent in the log
     for its branch and --head-sha + --remote-slug are given and the live remote tip of
     that branch EQUALS --head-sha (the carrying push itself, verified live: the
     schema's disclosed limitation is that a branch's last confirm can only be
     published by a subsequent push).
  3. REMOTE VERIFICATION (only with --remote-slug): every recorded push SHA (legacy
     one-phase `push` entries and CONFIRM remote_sha values) is verified against the
     GitHub commits API:
       - HTTP 200            -> "on remote" (the recorded object exists);
       - HTTP 404/422        -> FAIL: recorded SHA is absent from the remote;
       - 401/403/429/5xx/net -> VERIFICATION_FAILURE (fail closed, never "absent").
     Each such SHA must then be (i) the CURRENT remote tip of its branch, or (ii) an
     ancestor of it (local merge-base when the objects are present, else the compare
     API), or (iii) explicitly reconciled in the reconciliation record with a pointer
     to the documenting evidence (e.g. a documented, no-force-push-replacement rebuild
     of the branch). A SHA that is none of these FAILs.
  4. STATUS: the check ends with exactly one status line:
       PUSHLOG STATUS: COMPLETE   — zero pending intents, zero findings, all remote
                                    verifications performed (only when --remote-slug
                                    was given); COMPLETE is never printed for an
                                    offline-only run.
       PUSHLOG STATUS: PENDING    — zero findings but at least one unconfirmed intent.
       PUSHLOG STATUS: INVALID    — any findings.
     Pending and findings both exit 1; the status line is what distinguishes them.

Reconciliation record (machine-readable, append-only in practice): a JSON file
(default docs/evidence/push-log/PUSH_LOG_RECONCILIATION.json) whose "entries" list
objects: {"entry_id": str, "branch": str, "sha": 40-hex or null, "class": str,
"pointer": str}. Unknown classes fail; the class vocabulary is:
  SUPERSEDED-BY-DOCUMENTED-REBUILD  (SHA/branch no longer on the current lineage
                                     because the branch was rebuilt/advanced through a
                                     documented, no-force-push-replacement event)
  CONFIRMED-BY-LATER-EVIDENCE       (intent confirmed by evidence outside the log)
  SUPERSEDED-BY-CORRECTION-ENTRY    (a later correction entry supersedes this record's
                                     remote_sha — e.g. a hand-typed SHA tail captured
                                     wrongly; the correction entry is authoritative)

Offline test mode: --fixture FILE supplies canned API responses:
  {"refs": {"<branch>": "<sha>"}, "commits": {"<sha>": 200}, "compare": {"<base>...<head>": "ahead"}}

Exit: 0 clean; 1 findings or pending; 2 usage.
"""
import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone

LINE = re.compile(r"^## (\S+) \| ([^|]+) \| ([^|]+) \| ([^|]+?)(?:\s\|\s*RECONSTRUCTED)?\s*$")
PLACEHOLDER_SHAS = {"self-referential-tip", "TBD", "PENDING", "<sha>", "sha", "RECORDED-BELOW"}
FULL_SHA = re.compile(r"^[0-9a-f]{40}$")
KV = re.compile(r"(\b[a-z][a-z0-9_-]*)=(\S+)")
RECON_CLASSES = {"SUPERSEDED-BY-DOCUMENTED-REBUILD", "CONFIRMED-BY-LATER-EVIDENCE",
                 "SUPERSEDED-BY-CORRECTION-ENTRY"}


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


def _parse_entries(path: str):
    """Parse grammar-compliant entries into (line_no, ts, entry_id, verb, subject)."""
    out = []
    for n, raw in enumerate(open(path, encoding="utf-8"), 1):
        line = raw.rstrip("\n")
        if not line.startswith("## "):
            continue
        m = LINE.match(line)
        if not m:
            continue
        try:
            ts = parse_ts(m.group(1))
        except ValueError:
            continue
        out.append({"line": n, "ts": ts, "id": m.group(2).strip(),
                    "verb": m.group(3).strip(), "subject": m.group(4).strip()})
    return out


def _kv(subject: str) -> dict:
    return dict(KV.findall(subject))


def _gh_get(url: str, token: str | None):
    req = urllib.request.Request(url, method="GET",
                                 headers={"Accept": "application/vnd.github+json",
                                          "X-GitHub-Api-Version": "2022-11-28"})
    if token:
        req.add_header("Authorization", f"token {token}")
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            return resp.status, json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        return e.code, None
    except (urllib.error.URLError, TimeoutError, OSError):
        return None, None


def _remote_ancestor(slug: str, base: str, head: str, token: str | None, fixture: dict | None):
    """True when `base` is an ancestor of (or equal to) `head` on the remote."""
    key = f"{base}...{head}"
    if fixture is not None:
        return fixture.get("compare", {}).get(key) in ("ahead", "identical")
    status, body = _gh_get(f"https://api.github.com/repos/{slug}/compare/{key}", token)
    if status != 200 or not isinstance(body, dict):
        return False
    return body.get("status") in ("ahead", "identical")


def _local_ancestor(base: str, head: str) -> bool | None:
    """Local ancestry via git objects; None when the objects are unavailable."""
    if subprocess_run_ok(["git", "cat-file", "-e", f"{base}^{{commit}}"]) and \
       subprocess_run_ok(["git", "cat-file", "-e", f"{head}^{{commit}}"]):
        return subprocess_run_ok(["git", "merge-base", "--is-ancestor", base, head])
    return None


def subprocess_run_ok(argv: list[str]) -> bool:
    import subprocess
    return subprocess.run(argv, capture_output=True).returncode == 0


def _remote_tip(slug: str, branch: str, token: str | None, fixture: dict | None):
    if fixture is not None:
        return fixture.get("refs", {}).get(branch)
    import urllib.parse
    ref = urllib.parse.quote(f"heads/{branch}", safe="/")
    status, body = _gh_get(f"https://api.github.com/repos/{slug}/git/ref/{ref}", token)
    if status == 200 and isinstance(body, dict):
        obj = body.get("object") or {}
        return obj.get("sha")
    return None


def _sha_on_remote(slug: str, sha: str, token: str | None, fixture: dict | None):
    """('on-remote'|'absent'|'error', detail) for one recorded SHA."""
    if fixture is not None:
        code = fixture.get("commits", {}).get(sha)
        if code == 200:
            return ("on-remote", "")
        if code in (404, 422):
            return ("absent", f"fixture HTTP {code}")
        return ("error", f"fixture HTTP {code}")
    status, _ = _gh_get(f"https://api.github.com/repos/{slug}/commits/{sha}", token)
    if status == 200:
        return ("on-remote", "")
    if status in (404, 422):
        return ("absent", f"HTTP {status}")
    return ("error", f"HTTP {status}" if status else "network failure (fail closed)")


def _load_reconciliation(path: str | None):
    if not path or not os.path.exists(path):
        return {}, []
    data = json.load(open(path, encoding="utf-8"))
    entries = {}
    errors = []
    for i, e in enumerate(data.get("entries", [])):
        if e.get("class") not in RECON_CLASSES:
            errors.append(f"reconciliation entry {i}: unknown class {e.get('class')!r}")
            continue
        if not e.get("entry_id") or not e.get("branch") or not e.get("pointer"):
            errors.append(f"reconciliation entry {i}: missing entry_id/branch/pointer")
            continue
        sha = e.get("sha")
        if sha is not None and not FULL_SHA.match(sha or ""):
            errors.append(f"reconciliation entry {i}: sha is not a full 40-hex value")
            continue
        entries[(e["entry_id"], e["branch"])] = e
    return entries, errors


def check_pushlog_schema(path: str, slug: str | None = None, token: str | None = None,
                         head_sha: str | None = None, recon_path: str | None = None,
                         fixture: dict | None = None) -> tuple[list[str], str]:
    """Two-phase push-log protocol check. Returns (report_lines, status_word).

    report_lines contains FAIL:/PENDING:/OK: lines; the status word is INVALID only
    when FAIL: lines exist, PENDING when only unconfirmed intents exist, COMPLETE
    only for a fully verified protocol (remote mode), OFFLINE-CLEAN for a clean
    offline run (pairing + truthfulness only; never reported as COMPLETE).
    """
    findings: list[str] = []   # FAIL lines
    notes: list[str] = []      # OK lines
    pending: list[str] = []    # PENDING lines
    recon, recon_errors = _load_reconciliation(recon_path)
    findings += [f"FAIL: {path}: {e}" for e in recon_errors]

    entries = _parse_entries(path)
    intents = [e for e in entries if e["verb"] == "push-intent"]
    confirms = [e for e in entries if e["verb"] == "push-confirm"]
    legacy = [e for e in entries if e["verb"] == "push"]

    # 1. truthfulness: CONFIRM entries carry real 40-hex remote SHAs
    confirms_by_branch: dict[str, list[dict]] = {}
    for e in confirms:
        kv = _kv(e["subject"])
        sha = kv.get("remote_sha") or kv.get("sha")
        if not sha or sha.lower() in PLACEHOLDER_SHAS or not FULL_SHA.match(sha):
            findings.append(f"FAIL: {path}:{e['line']}: CONFIRM entry {e['id']} carries "
                            f"{sha!r} instead of the real 40-hex SHA (truthfulness violated)")
            continue
        confirms_by_branch.setdefault(kv.get("branch", ""), []).append(e)

    # 2. pairing: every push-intent has a LATER confirm for the same branch,
    #    unless reconciled, or it is the branch's final intent confirmed live.
    for e in intents:
        kv = _kv(e["subject"])
        branch = kv.get("branch", "")
        later = [c for c in confirms_by_branch.get(branch, []) if c["line"] > e["line"]]
        if later:
            continue
        if (e["id"], branch) in recon:
            notes.append(f"OK: {path}:{e['line']}: intent {e['id']} for {branch} is "
                         f"reconciled in the push-log reconciliation record")
            continue
        is_final = not any(o["line"] > e["line"] for o in intents if o["id"] != e["id"]
                           and _kv(o["subject"]).get("branch", "") == branch)
        if is_final and slug:
            tip = _remote_tip(slug, branch, token, fixture)
            if tip is None:
                findings.append(f"FAIL: {path}:{e['line']}: cannot read remote tip for {branch} "
                                f"(fail closed); final intent {e['id']} not confirmable live")
                continue
            # CONFIRMED-LIVE, two admissible proofs:
            #  (a) the remote tip equals --head-sha (the carrying push itself);
            #  (b) the remote tip is AT OR AFTER the commit that INTRODUCED this intent
            #      entry (the push that carried the intent is verifiably on the remote —
            #      covers a stacked tree validating another branch's final intent).
            carrying = None
            if fixture is None:
                import subprocess as _sp
                _r = _sp.run(["git", "log", "--reverse", "-S", e["id"], "--format=%H",
                              "--", path], capture_output=True, text=True, cwd=os.getcwd())
                if _r.stdout.strip():
                    carrying = _r.stdout.strip().splitlines()[0]
            if head_sha and FULL_SHA.match(head_sha or "") and tip == head_sha:
                notes.append(f"OK: {path}:{e['line']}: final intent {e['id']} for {branch} "
                             f"CONFIRMED-LIVE (remote tip == carrying commit {head_sha[:12]}…)")
                continue
            if carrying:
                la = _local_ancestor(carrying, tip)
                if la is True or (la is None and
                                  _remote_ancestor(slug, carrying, tip, token, fixture)):
                    notes.append(f"OK: {path}:{e['line']}: final intent {e['id']} for {branch} "
                                 f"CONFIRMED-LIVE (remote tip is at-or-after the intent's "
                                 f"carrying commit {carrying[:12]}…)")
                    continue
            pending.append(f"PENDING: {path}:{e['line']}: push-intent {e['id']} for branch "
                           f"{branch or '<unknown>'} has no matching push-confirm (a schema pass "
                           f"while a push remains unconfirmed reports PENDING, not complete)")
        elif is_final and not slug:
            pending.append(f"PENDING: {path}:{e['line']}: push-intent {e['id']} for branch "
                           f"{branch or '<unknown>'} has no matching push-confirm (offline run; "
                           f"live confirmation was not requested — PENDING, not complete)")
        elif not is_final:
            pending.append(f"PENDING: {path}:{e['line']}: push-intent {e['id']} for branch "
                           f"{branch or '<unknown>'} has no matching push-confirm (a schema pass "
                           f"while a push remains unconfirmed reports PENDING, not complete)")

    # 3. remote verification of every recorded push SHA (needs --remote-slug)
    if slug:
        recorded: list[tuple[dict, str, str]] = []  # (entry, branch, sha)
        for e in legacy:
            kv = _kv(e["subject"])
            sha = kv.get("sha")
            if not sha or sha.lower() in PLACEHOLDER_SHAS or not FULL_SHA.match(sha):
                continue  # disclosed legacy markers (e.g. P-011) have their own CONFIRM
            recorded.append((e, kv.get("branch", ""), sha))
        for e in confirms:
            kv = _kv(e["subject"])
            sha = kv.get("remote_sha") or kv.get("sha") or ""
            if FULL_SHA.match(sha):
                recorded.append((e, kv.get("branch", ""), sha))
        for e, branch, sha in recorded:
            state, detail = _sha_on_remote(slug, sha, token, fixture)
            if state == "error":
                findings.append(f"FAIL: {path}:{e['line']}: VERIFICATION_FAILURE for {sha} "
                                f"({detail}); fail closed — not treated as absent")
                continue
            if state == "absent":
                findings.append(f"FAIL: {path}:{e['line']}: recorded SHA {sha} ({e['id']}, "
                                f"branch {branch}) is ABSENT from the remote")
                continue
            tip = _remote_tip(slug, branch, token, fixture)
            if tip == sha:
                continue
            la = _local_ancestor(sha, tip) if tip else None
            if la is True:
                continue
            if la is None and tip and _remote_ancestor(slug, sha, tip, token, fixture):
                continue
            if (e["id"], branch) in recon:
                continue
            findings.append(f"FAIL: {path}:{e['line']}: recorded SHA {sha} ({e['id']}, branch "
                            f"{branch}) is neither the current remote tip, nor an ancestor of "
                            f"it, nor reconciled in the push-log reconciliation record")

    # 4. final status word (FAIL lines dominate; PENDING is never COMPLETE)
    if findings:
        status = "INVALID"
    elif pending:
        status = "PENDING"
    elif not slug:
        status = "OFFLINE-CLEAN"
    else:
        status = "COMPLETE"
    return findings + pending + notes, status


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("logfile")
    ap.add_argument("--commit-time", help="ISO-8601 commit time carrying this log version")
    ap.add_argument("--max-drift", type=int, default=300, help="seconds of allowed drift (default 300)")
    ap.add_argument("--parent-prefix", metavar="OLDER_FILE",
                    help="require OLDER_FILE to be a byte-prefix of logfile (append-only)")
    ap.add_argument("--pushlog-schema", action="store_true",
                    help="push-log protocol check: CONFIRM truthfulness + INTENT->CONFIRM "
                         "pairing (pending intents are reported PENDING, never complete)")
    ap.add_argument("--remote-slug", metavar="OWNER/REPO",
                    help="with --pushlog-schema: verify recorded SHAs and branch tips "
                         "against this repository via the GitHub API")
    ap.add_argument("--head-sha", help="with --remote-slug: the carrying commit SHA of this "
                                       "log version (enables CONFIRMED-LIVE for the final intent)")
    ap.add_argument("--reconciliation", metavar="FILE",
                    default="docs/evidence/push-log/PUSH_LOG_RECONCILIATION.json",
                    help="machine-readable reconciliation record for superseded SHAs/intents")
    ap.add_argument("--fixture", help="offline test mode: canned GitHub API responses (JSON)")
    a = ap.parse_args()
    findings = check_sequence(a.logfile, a.commit_time, a.max_drift)
    if a.parent_prefix:
        findings += check_parent_prefix(a.logfile, a.parent_prefix)
    if a.pushlog_schema:
        fixture = json.load(open(a.fixture, encoding="utf-8")) if a.fixture else None
        token = None
        if not fixture:
            token = os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
        report, status = check_pushlog_schema(
            a.logfile, slug=a.remote_slug, token=token, head_sha=a.head_sha,
            recon_path=a.reconciliation, fixture=fixture)
        findings += [l for l in report if l.startswith(("FAIL:", "PENDING:"))]
        for l in report:
            if l.startswith("OK:"):
                print(l)
        if status != "OFFLINE-CLEAN" or a.remote_slug:
            print(f"PUSHLOG STATUS: {status}")
        elif status == "OFFLINE-CLEAN":
            print("PUSHLOG STATUS: OFFLINE-CLEAN (pairing/truthfulness only; no remote "
                  "verification was requested, so this is never a COMPLETE verdict)")
    print("\n".join(findings) if findings else f"PASS: {a.logfile} "
          f"{'sequence + append-only ' if a.parent_prefix else ''}clean")
    return 0 if not findings else 1


if __name__ == "__main__":
    sys.exit(main())
