#!/usr/bin/env python3
"""EC-18/EC-24 evidence sequence + append-only + remote-ref verification checker.

Revision history:
  Round-19/20: append-only parent-prefix comparison; two-phase INTENT/CONFIRM
    schema; push-reconcile verb; live remote-ref verification (CURRENT vs
    HISTORICAL); evidence-support and runs-manifest cross-checks.

  Round-21 hardening (auditor directions 4-6):
    * IDENTIFIER PAIRING (direction 4): a CONFIRM pairs the INTENT (or legacy
      placeholder push) with the SAME identifier base for the SAME branch --
      P-058-CONFIRM consumes P-058-INTENT, never another intent. FIFO/branch
      pairing alone is no longer accepted.
    * TERMINAL-INTENT EXACT-CANDIDATE BINDING (direction 4): for the branch
      under test, at most ONE unmatched intent -- the terminal one -- may
      represent the current candidate's push, and only if (a) the commit that
      introduced the intent IS the exact candidate under test (an intent can
      never declare later commits merely by remaining unmatched), (b) its
      expected_remote_tip_before_push equals the newest CONFIRM/RECONCILE
      state recorded for the branch before it (no intervening ref update) and
      is an ancestor of the candidate, or is ABSENT with no prior branch
      record, and (c) the live remote tip is the exact expected tip. That
      state is labeled PENDING, never COMPLETE. Excess unmatched intents FAIL
      unless separately voided by an intent-reconcile record.
    * HISTORICAL PROOF (direction 5): a confirmation whose SHA differs from
      the live tip is no longer classified HISTORICAL on that fact alone. The
      SHA must be proven a real commit and supported: by ancestry from the
      branch's live tip, or -- when the SHA left the branch lineage -- by a
      reviewed evidence-back record whose SHA is a real remote commit
      (object store or commits API via --repo+token) with verified run
      citations.
    * MERGE TRANSITION (direction 6): a push to a protected branch is
      declared ONLY by a VERIFIED merge-intent: the candidate must be a merge
      commit whose first parent is the declared pre-merge base and whose
      second parent is the declared (approved) PR head, the merge-intent must
      ride the PR branch history, and the PR must be merged with an approving
      review from a user other than the PR author (GitHub API, or fixture
      files in tests). On the source branch, the active merge-intent must
      bind to the exact candidate and the declared base must still be the
      live tip of the target. Raced or unrelated updates to the protected
      branch FAIL.

  Round-22 hardening (auditor directions 3-5):
    * CLAIM SEPARATION (direction 3): commit existence, branch-push
      provenance, and CI success are now verified and reported as SEPARATE
      claims. The runs manifest and the Actions-API verifier retain and
      validate each cited run's EVENT TYPE and HEAD_BRANCH alongside
      head_sha and conclusion; a push-event run on a different branch
      contradicts the entry and FAILs; a failed/cancelled run is disclosed
      and can never be represented as successful gate evidence.
    * INDEPENDENT PUSH PROVENANCE (direction 3): an off-lineage HISTORICAL
      confirmation must prove the BRANCH PUSH, not merely commit existence:
      (a) a push-event Actions run (event=push, head_branch, head_sha), or
      (b) a live Events-API PushEvent (ref + before/after), or
      (c) an evidence-back carrying a recorded push_event= citation
      (id/ref/before/after) validated by a machine-verifiable review_ref=
      into the canonical auditor-adjudicated record, or
      (d) an explicit push-not-proven reclassification (valid review_ref=):
      the push is NOT PROVEN -- disclosed, never counted as evidence.
      Otherwise FAIL. A commits-API 200 or a free-text rebuild narrative is
      NOT push provenance.
    * INTENT RECONCILIATION (direction 4): every intent-reconcile record
      must target a PRECEDING eligible intent, match its branch semantics
      (push-intent: branch; merge-intent: source_branch), void it exactly
      once, and never void a consumed/confirmed intent. The unmatched-intent
      invariant is enforced GLOBALLY: no branch may hold more than one
      unvoided unmatched intent.
    * PRE-PUSH STATE BINDING (direction 4): the terminal intent's
      expected_remote_tip_before_push is bound to the TRUSTED push-event
      before/after state (GitHub PushEvent): ref must match the branch,
      after must be the candidate, before must equal the declared pre-push
      tip (or the all-zero SHA for a newly created branch with an ABSENT
      declaration). PR-event runs never substitute for push-event
      provenance.
    * MERGE-TRANSITION PROOF (direction 5): the PR API's merge_commit_sha
      must equal the protected branch's live candidate; source head.ref,
      exact head SHA, protected base.ref, merged state, merged_at, and a
      FRESH non-author approval (latest review state per reviewer, on the
      exact head, submitted before the merge) are all verified; a
      pre-merge verification mode (--pre-merge-verify) checks the current
      PR/head/review state before any merge executes.

Log grammar:   ## <ISO-8601 UTC> | <ENTRY-ID> | <verb> | <subject> [| RECONSTRUCTED]
Verbs: push, push-intent, push-confirm, push-reconcile (schema v3),
  merge-intent, evidence-back, intent-reconcile (schema v4),
  push-not-proven (schema v5).
Sequence checks: timestamps parse as UTC, non-decreasing, ENTRY-IDs unique; with
--commit-time, entries whose claimed time drifts more than --max-drift from the
carrying commit must be RECONSTRUCTED.

Exit: 0 clean; 1 findings; 2 usage.
"""
import argparse
import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.request

LINE = re.compile(r"^## (\S+) \| ([^|]+) \| ([^|]+) \| ([^|]+?)(?:\s\|\s*RECONSTRUCTED)?\s*$")
PLACEHOLDER_SHAS = {"self-referential-tip", "TBD", "PENDING", "<sha>", "sha"}
SHA40 = re.compile(r"\b([0-9a-f]{40})\b")
ZERO_SHA = "0" * 40
RUN_ID = re.compile(r"\brun(?:_id)?\s*[=:]?\s*(\d{6,})\b", re.I)
PUSH_EVENT_CITE = re.compile(
    r"push_event=(\d+)\s+ref=(\S+)\s+before=([0-9a-f]{40}|0{40})\s+after=([0-9a-f]{40})")
# Round-22 prior-tip citation: a LATER push on the same ref whose `before`
# field records the backed SHA -- GitHub's own record that the ref pointed
# there (refs move only by pushes), used when the push's own event has
# rolled off the retained window and the API elides its `after`.
PRIOR_TIP_CITE = re.compile(
    r"push_event=(\d+)\s+ref=(\S+)\s+prior_tip=([0-9a-f]{40})")
CONFIRM_VERBS = ("push-confirm", "push-reconcile")
API = "https://api.github.com"


def parse_ts(s: str):
    from datetime import datetime, timezone
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


def id_base(eid: str) -> str:
    """Identifier base for pairing: P-058-CONFIRM / P-058-INTENT / P-058 -> P-058."""
    return re.sub(r"-(INTENT|CONFIRM)$", "", eid)


def field(subject: str, name: str):
    m = re.search(rf"(?:^|\s){name}=(\"[^\"]*\"|\S+)", subject)
    if not m:
        return None
    v = m.group(1)
    return v[1:-1] if v.startswith('"') else v


# ------------------------------------------------------------- review references

def validate_review_ref(logrepo, review_ref: str, needles: list[str]) -> str | None:
    """Round-22 direction 3.5: a review/decision reference must resolve to a
    tracked canonical record whose adjudicating section mentions the entry.

    Format: <repo-relative-path>#<anchor>, e.g.
      docs/decisions/OWNER_DECISIONS.md#Part-14
    The anchor names a `## ` heading (dashes read as spaces); the section runs
    from that heading to the next `## ` heading. At least one of NEEDLES (the
    backed/reclassified entry IDs or SHAs) must appear in that section, so
    "reviewed" is never merely a string written by the same record's author:
    it must point at the canonical auditor-adjudicated record. Returns an
    error string, or None when the reference validates."""
    if not review_ref:
        return "review_ref is missing"
    path, _, anchor = review_ref.rpartition("#")
    if not path or not anchor:
        return (f"review_ref {review_ref!r} is malformed (need "
                f"<repo-relative-path>#<heading-anchor>)")
    if not logrepo.ok:
        return ("review_ref cannot be validated -- the push log is not inside "
                "a git repository (fail closed)")
    full = os.path.join(logrepo.toplevel, path)
    if not os.path.exists(full):
        return f"review_ref names a file absent from the candidate tree: {path}"
    try:
        text = open(full, encoding="utf-8").read()
    except OSError:
        return f"review_ref names an unreadable file: {path}"
    want = anchor.replace("-", " ").strip().lower()
    lines = text.splitlines()
    start = None
    for i, ln in enumerate(lines):
        if ln.startswith("## "):
            heading = ln[3:].strip().lower()
            if want and want in heading:
                start = i
                break
    if start is None:
        return (f"review_ref anchor {anchor!r} does not name any '## ' heading "
                f"in {path}")
    section: list[str] = []
    for ln in lines[start + 1:]:
        if ln.startswith("## "):
            break
        section.append(ln)
    body = "\n".join(section)
    if not any(nd in body for nd in needles if nd):
        return (f"review_ref {review_ref!r} resolves, but its adjudicating "
                f"section does not mention the backed entry "
                f"({', '.join(n for n in needles if n)}) -- a review reference "
                f"must bind to the record that actually adjudicates this entry")
    return None


# ---------------------------------------------------------------- sequence / prefix

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
        i = next((k for k in range(len(old)) if old[k] != cur[k]), len(old))
        ctx = old[max(0, i - 40):i + 40].decode("utf-8", "replace")
        return [f"FAIL: {current}: prior content MODIFIED at byte {i} "
                f"(append-only violated); context: {ctx!r}"]
    return []


# ---------------------------------------------------------------- schema (pairing by ID)

def check_pushlog_schema(path: str, logrepo=None) -> tuple[list[str], list[str], dict]:
    """Schema + identifier-based INTENT/CONFIRM pairing + Round-22
    reconciliation/reclassification validation. Returns (findings, pending, state)."""
    findings: list[str] = []
    pending: list[str] = []
    entries = parse_entries(path)

    def subj_field(e, key):
        return field(e["subject"], key)

    # 1. per-verb required fields + value grammar
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
            for key in ("branch=", "remote_sha=", "credential=", "evidence="):
                if key not in subject:
                    findings.append(f"FAIL: {path}:{e['n']}: RECONCILE entry is missing the "
                                    f"required field {key[:-1]!r}")
            if not re.search(r"remote_sha=([0-9a-f]{40})(?:\s|$)", subject):
                findings.append(f"FAIL: {path}:{e['n']}: RECONCILE entry remote_sha must be a "
                                f"real 40-hex SHA (never a placeholder)")
        elif verb == "push-intent":
            for key in ("branch=", "credential="):
                if key not in subject:
                    findings.append(f"FAIL: {path}:{e['n']}: push-intent entry is missing the "
                                    f"required field {key[:-1]!r}")
            ert = subj_field(e, "expected_remote_tip_before_push")
            if ert is not None:
                tok = ert.split()[0] if ert.split() else ""
                if not (re.fullmatch(r"[0-9a-f]{40}", tok) or tok == "ABSENT"):
                    findings.append(f"FAIL: {path}:{e['n']}: expected_remote_tip_before_push "
                                    f"must be a real 40-hex SHA or ABSENT, not {tok!r}")
        elif verb == "push":
            for key in ("branch=", "credential="):
                if key not in subject:
                    findings.append(f"FAIL: {path}:{e['n']}: push entry is missing the "
                                    f"required field {key[:-1]!r}")
        elif verb == "merge-intent":
            for key in ("target=", "pr=", "source_branch=", "expected_pre_merge_base=",
                        "expected_head=", "credential=", "authorized_by="):
                if key not in subject:
                    findings.append(f"FAIL: {path}:{e['n']}: merge-intent entry is missing the "
                                    f"required field {key[:-1]!r}")
            for f_name in ("expected_pre_merge_base",):
                v = subj_field(e, f_name)
                if v is not None and not re.fullmatch(r"[0-9a-f]{40}", v):
                    findings.append(f"FAIL: {path}:{e['n']}: merge-intent {f_name} must be a "
                                    f"real 40-hex SHA")
            eh = subj_field(e, "expected_head")
            if eh is not None and not (re.fullmatch(r"[0-9a-f]{40}", eh)
                                       or eh == "DEFERRED"):
                findings.append(f"FAIL: {path}:{e['n']}: merge-intent expected_head must be "
                                f"a real 40-hex SHA or DEFERRED (the head is the commit "
                                f"that introduces this record; a record cannot embed its "
                                f"own commit's SHA)")
            pr = subj_field(e, "pr")
            if pr is not None and not re.fullmatch(r"\d+", pr):
                findings.append(f"FAIL: {path}:{e['n']}: merge-intent pr= must be a number")
        elif verb == "evidence-back":
            for key in ("backs=", "branch=", "remote_sha=", "credential=", "evidence="):
                if key not in subject:
                    findings.append(f"FAIL: {path}:{e['n']}: evidence-back entry is missing the "
                                    f"required field {key[:-1]!r}")
            if not re.search(r"remote_sha=([0-9a-f]{40})(?:\s|$)", subject):
                findings.append(f"FAIL: {path}:{e['n']}: evidence-back remote_sha must be a "
                                f"real 40-hex SHA")
        elif verb == "intent-reconcile":
            for key in ("voids=", "branch=", "credential=", "evidence="):
                if key not in subject:
                    findings.append(f"FAIL: {path}:{e['n']}: intent-reconcile entry is missing "
                                    f"the required field {key[:-1]!r}")
        elif verb == "push-not-proven":
            for key in ("reclassifies=", "branch=", "remote_sha=", "credential=",
                        "review_ref=", "evidence="):
                if key not in subject:
                    findings.append(f"FAIL: {path}:{e['n']}: push-not-proven entry is missing "
                                    f"the required field {key[:-1]!r}")
            if not re.search(r"remote_sha=([0-9a-f]{40})(?:\s|$)", subject):
                findings.append(f"FAIL: {path}:{e['n']}: push-not-proven remote_sha must be a "
                                f"real 40-hex SHA")
        else:
            findings.append(f"FAIL: {path}:{e['n']}: unknown push-log verb {verb!r}")

    # 2. duplicate ENTRY-IDs
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

    by_eid = {e["eid"]: e for e in entries if e["m"]}

    # 2b. INTENT-RECONCILE validation (Round-22 direction 4): a void must
    # target a PRECEDING eligible intent, match its branch semantics, and
    # void it exactly once. Cross-branch voids, duplicate/conflicting voids,
    # voids of consumed (confirmed) intents, and invalid target verbs FAIL.
    voided_intents: set[str] = set()
    void_targets: dict[str, dict] = {}
    for e in entries:
        if not e["m"] or e["verb"] != "intent-reconcile":
            continue
        voids = field(e["subject"], "voids")
        target = by_eid.get(voids or "")
        if target is None or target["verb"] not in ("push-intent", "merge-intent"):
            findings.append(f"FAIL: {path}:{e['n']}: intent-reconcile {e['eid']} names "
                            f"voids={voids!r} which is not an intent entry (dangling void)")
            continue
        if target["n"] >= e["n"]:
            findings.append(f"FAIL: {path}:{e['n']}: intent-reconcile {e['eid']} voids "
                            f"{voids} which does not PRECEDE it (a reconciliation "
                            f"must target a preceding eligible intent)")
            continue
        # branch scope: push-intent -> its own branch; merge-intent -> its
        # source_branch (unambiguous source/target semantics, not a free-form
        # branch value)
        if target["verb"] == "push-intent":
            scope = field(target["subject"], "branch")
        else:
            scope = field(target["subject"], "source_branch")
        rbranch = field(e["subject"], "branch")
        if scope != rbranch:
            findings.append(f"FAIL: {path}:{e['n']}: intent-reconcile {e['eid']} declares "
                            f"branch={rbranch!r} but the voided intent {voids} belongs to "
                            f"branch scope {scope!r} (cross-branch voiding rejected)")
            continue
        if voids in void_targets:
            findings.append(f"FAIL: {path}:{e['n']}: intent-reconcile {e['eid']} voids "
                            f"{voids} a second time (first at line {void_targets[voids]['n']}); "
                            f"an intent is voided exactly once -- duplicate or conflicting "
                            f"voids are rejected)")
            continue
        # a reconciliation must not substitute for or silently consume a
        # confirmation: voiding an intent that a PRECEDING CONFIRM already
        # completed is rejected (a void followed by its CONFIRM is rejected
        # by the pairing rule instead -- both directions fail closed)
        consumed = any(x["m"] and x["verb"] == "push-confirm"
                       and id_base(x["eid"]) == id_base(voids)
                       and field(x["subject"], "branch") == rbranch
                       and x["n"] < e["n"]
                       for x in entries)
        if consumed:
            findings.append(f"FAIL: {path}:{e['n']}: intent-reconcile {e['eid']} voids "
                            f"{voids}, which is already completed by its CONFIRM "
                            f"(a reconciliation must not consume or substitute for a "
                            f"confirmation; void a live intent)")
            continue
        void_targets[voids] = e
        voided_intents.add(voids)

    # 2c. PUSH-NOT-PROVEN reclassification validation (Round-22 direction 3):
    # reclassifies a CONFIRM/RECONCILE with matching branch/sha, exactly once,
    # with a review_ref that resolves into the canonical adjudicated record.
    reclassified: set[str] = set()
    reclass_targets: dict[str, dict] = {}
    for e in entries:
        if not e["m"] or e["verb"] != "push-not-proven":
            continue
        rc = field(e["subject"], "reclassifies")
        target = by_eid.get(rc or "")
        if target is None or target["verb"] not in CONFIRM_VERBS:
            findings.append(f"FAIL: {path}:{e['n']}: push-not-proven {e['eid']} names "
                            f"reclassifies={rc!r} which is not a CONFIRM/RECONCILE entry "
                            f"(dangling reclassification)")
            continue
        if target["n"] >= e["n"]:
            findings.append(f"FAIL: {path}:{e['n']}: push-not-proven {e['eid']} reclassifies "
                            f"{rc} which does not PRECEDE it")
            continue
        tb = field(target["subject"], "branch")
        ts_ = field(target["subject"], "remote_sha")
        if tb != field(e["subject"], "branch") or ts_ != field(e["subject"], "remote_sha"):
            findings.append(f"FAIL: {path}:{e['n']}: push-not-proven {e['eid']} branch/sha "
                            f"do not match the reclassified entry {rc} "
                            f"(a reclassification must describe the same push)")
            continue
        if rc in reclass_targets:
            findings.append(f"FAIL: {path}:{e['n']}: push-not-proven {e['eid']} reclassifies "
                            f"{rc} a second time (first at line {reclass_targets[rc]['n']}); "
                            f"exactly one reclassification per confirmation")
            continue
        rr_err = validate_review_ref(logrepo, field(e["subject"], "review_ref"),
                                      [e["eid"], rc, (ts_ or "")[:12]])
        if rr_err:
            findings.append(f"FAIL: {path}:{e['n']}: push-not-proven {e['eid']}: {rr_err}")
            continue
        reclass_targets[rc] = e
        reclassified.add(rc)

    # 3. IDENTIFIER pairing (Round-21 direction 4 + Round-22 void rules): a
    # CONFIRM consumes the unpaired, UNVOIDED intent (or legacy placeholder
    # push) with the SAME id base for the SAME branch. A voided intent is not
    # eligible for pairing: a reconciliation must not silently consume or
    # substitute for a confirmation.
    open_by_base: dict[str, list[dict]] = {}
    for e in entries:
        if not e["m"]:
            continue
        verb, subject = e["verb"], e["subject"]
        branch = field(subject, "branch")
        if verb == "push-intent":
            open_by_base.setdefault(id_base(e["eid"]), []).append(
                {"e": e, "branch": branch, "placeholder": False})
        elif verb == "push":
            sm = re.search(r"sha=(\S+)", subject)
            if branch and sm and sm.group(1).rstrip('",').lower() in PLACEHOLDER_SHAS:
                open_by_base.setdefault(id_base(e["eid"]), []).append(
                    {"e": e, "branch": branch, "placeholder": True})
    unpaired: list[dict] = []
    for e in entries:
        if not e["m"] or e["verb"] != "push-confirm":
            continue
        branch = field(e["subject"], "branch")
        base = id_base(e["eid"])
        queue = open_by_base.get(base, [])
        match = next((q for q in queue
                      if not q.get("consumed") and q["branch"] == branch
                      and q["e"]["eid"] not in voided_intents), None)
        if match is not None:
            match["consumed"] = True
        else:
            voided_here = [q for q in queue
                           if not q.get("consumed") and q["branch"] == branch
                           and q["e"]["eid"] in voided_intents]
            if voided_here:
                findings.append(
                    f"FAIL: {path}:{e['n']}: CONFIRM {e['eid']} pairs the intent "
                    f"{voided_here[0]['e']['eid']}, which was VOIDED by an "
                    f"intent-reconcile record (a voided intent is no longer "
                    f"eligible; confirm a live intent)")
            else:
                findings.append(
                    f"FAIL: {path}:{e['n']}: CONFIRM {e['eid']} has no preceding unconfirmed "
                    f"INTENT (or placeholder push) with the identifier {base} for branch "
                    f"{branch!r} (Round-21 identifier pairing: a confirmation must match its "
                    f"own intent and must not consume another)")
    for q in [q for qs in open_by_base.values() for q in qs
              if not q.get("consumed") and q["e"]["eid"] not in voided_intents]:
        unpaired.append(q)
        kind = "legacy push" if q["placeholder"] else "push-intent"
        pending.append(f"PENDING: {path}:{q['e']['n']}: {q['e']['eid']} ({kind}) for branch "
                       f"{q['branch']} has no CONFIRM entry (two-phase protocol incomplete; "
                       f"do not imply completion)")

    # 3b. GLOBAL unmatched-intent invariant (Round-22 direction 4): no branch
    # may hold more than one unvoided, unconsumed intent. Only one valid
    # terminal intent per branch may remain pending; excess intents must be
    # validly reconciled. This applies to EVERY branch, not only the branch
    # under test.
    by_branch: dict[str, list[dict]] = {}
    for q in unpaired:
        if q["branch"]:
            by_branch.setdefault(q["branch"], []).append(q)
    for branch, qs in sorted(by_branch.items()):
        if len(qs) > 1:
            ids = ", ".join(q["e"]["eid"] for q in qs)
            findings.append(f"FAIL: {path}: multiple unmatched intents for branch "
                            f"{branch}: {ids} -- only one valid terminal intent per "
                            f"branch may remain pending; reconcile the excess "
                            f"intents (intent-reconcile) or complete them")

    # 4. evidence-back target integrity (+ Round-22 review_ref validation
    # when the reference is present; the provenance REQUIREMENT is enforced
    # in the remote-verification layer, which knows the API/manifest sources)
    for e in entries:
        if not e["m"] or e["verb"] != "evidence-back":
            continue
        backs = field(e["subject"], "backs")
        target = by_eid.get(backs or "")
        if target is None or target["verb"] not in CONFIRM_VERBS:
            findings.append(f"FAIL: {path}:{e['n']}: evidence-back {e['eid']} names "
                            f"backs={backs!r} which is not a CONFIRM/RECONCILE entry "
                            f"(dangling backing)")
        else:
            tb = field(target["subject"], "branch")
            ts_ = field(target["subject"], "remote_sha")
            if tb != field(e["subject"], "branch") or ts_ != field(e["subject"], "remote_sha"):
                findings.append(f"FAIL: {path}:{e['n']}: evidence-back {e['eid']} branch/sha "
                                f"do not match the backed entry {backs} "
                                f"(backing must describe the same push)")
        rr = field(e["subject"], "review_ref")
        if rr is not None:
            rr_err = validate_review_ref(logrepo, rr, [e["eid"], backs or "",
                                                        (field(e["subject"], "remote_sha") or "")[:12]])
            if rr_err:
                findings.append(f"FAIL: {path}:{e['n']}: evidence-back {e['eid']}: {rr_err}")
        # a recorded push_event= citation must be syntactically valid
        ev_text = (field(e["subject"], "evidence") or "") + "\n" + "\n".join(e["body"])
        for m in PUSH_EVENT_CITE.finditer(ev_text):
            if m.group(3) == m.group(4):
                findings.append(f"FAIL: {path}:{e['n']}: evidence-back {e['eid']} records a "
                                f"push_event citation whose before equals its after "
                                f"(malformed ref-transition evidence)")

    state = {"entries": entries, "unpaired": unpaired, "by_eid": by_eid,
             "voided_intents": voided_intents, "reclassified": reclassified}
    return findings, pending, state


# ---------------------------------------------------------------- git + remote helpers

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


def ls_remote_archive(remote: str, branch: str) -> tuple[str | None, str | None, str | None]:
    """Return (commit_sha, tag_ref, error): the commit an archive tag preserves
    for a PRUNED branch.

    D-V3/D-V4 (auditor, Round-26): a branch that is pruned is first preserved
    as an ANNOTATED tag at its recorded SHA, and the push-log checker must
    resolve that tag as remote evidence. Two spellings are accepted:
    refs/tags/archive/<branch> (the D-V4 policy spelling) and
    refs/tags/archive/<branch-with-slashes-as-hyphens> (the spelling the
    auditor named for repair/r23-governance-recovery). The peeled commit of
    the annotated tag is authoritative; a lightweight tag resolves to its
    commit directly.
    """
    candidates = [f"refs/tags/archive/{branch}",
                  f"refs/tags/archive/{branch.replace('/', '-')}"]
    try:
        # A glob over the archive namespace is required: the peeled
        # ("^{}") line of an annotated tag only appears when the pattern
        # matches it too, and the peeled commit -- not the tag object -- is
        # the SHA the branch pointed at.
        r = subprocess.run(["git", "ls-remote", "--tags", remote, "refs/tags/archive/*"],
                           capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired) as e:
        return None, None, f"ls-remote failed for the archive tag of {branch!r} ({e})"
    if r.returncode != 0:
        return None, None, (f"ls-remote exited {r.returncode} for the archive tag of "
                            f"{branch!r}: {r.stderr.strip()[:120]}")
    direct: dict[str, str] = {}
    peeled: dict[str, str] = {}
    for line in r.stdout.splitlines():
        parts = line.split()
        if len(parts) != 2 or not re.fullmatch(r"[0-9a-f]{40}", parts[0]):
            continue
        sha, ref = parts
        if ref.endswith("^{}"):
            peeled[ref[:-3]] = sha
        else:
            direct[ref] = sha
    for ref in candidates:
        if ref in peeled:
            return peeled[ref], ref, None
        if ref in direct:
            return None, ref, (f"the archive tag {ref} carries no peeled commit "
                                f"(an archive tag must be ANNOTATED, so the preserved "
                                f"commit is provable)")
    return None, None, f"no archive tag for {branch!r} on {remote}"


class LogRepo:
    """The git repository containing the push log (for ancestry, object and
    introducing-commit queries). A log outside any repository is legal for
    sequence/schema checks but fails closed wherever git proof is required."""

    def __init__(self, log_path: str):
        self.ok = False
        self.toplevel = None
        self.rel = None
        self.log_dir = os.path.dirname(os.path.abspath(log_path))
        self.log_name = os.path.basename(log_path)
        try:
            r = subprocess.run(["git", "-C", self.log_dir,
                                "rev-parse", "--show-toplevel"],
                               capture_output=True, text=True, timeout=30)
            if r.returncode == 0 and r.stdout.strip():
                self.toplevel = r.stdout.strip()
                try:
                    self.rel = os.path.relpath(
                        os.path.realpath(log_path),
                        os.path.realpath(self.toplevel)).replace("\\", "/")
                except (OSError, ValueError):
                    self.rel = self.log_name
                self.ok = True
        except (OSError, subprocess.TimeoutExpired):
            pass

    def _git(self, *args):
        return subprocess.run(["git", "-C", self.toplevel, *args],
                              capture_output=True, text=True, timeout=60)

    def _git_in_log_dir(self, *args):
        """Run git FROM the log's own directory.

        The pickaxe query needs a pathspec. A path computed from the OS path
        string breaks whenever the two spellings of the same directory differ
        -- on Windows the 8.3 short form of a temp directory (C:\\Users\\ADMINI~1\\...)
        versus the long form git reports (C:/Users/Administrator/...) -- which
        produced a pathspec OUTSIDE the repository, made every
        introducing-commit lookup fail closed on a perfectly valid log, and
        turned the whole push-log gate red on that platform. Running git from
        the log's directory with the file name as the pathspec is
        spelling-independent.
        """
        return subprocess.run(["git", *args], cwd=self.log_dir,
                              capture_output=True, text=True, timeout=60)

    def object_exists(self, sha: str) -> bool:
        return self.ok and self._git("cat-file", "-e", f"{sha}^{{commit}}").returncode == 0

    def is_ancestor(self, sha: str, tip: str) -> bool | None:
        if not self.ok:
            return None
        if not (self.object_exists(sha) and self.object_exists(tip)):
            return None
        return self._git("merge-base", "--is-ancestor", sha, tip).returncode == 0

    def parents(self, sha: str) -> list[str] | None:
        if not self.ok or not self.object_exists(sha):
            return None
        r = self._git("rev-list", "--parents", "-n", "1", sha)
        if r.returncode != 0:
            return None
        return r.stdout.split()

    def introducing_commit(self, eid: str, verb: str) -> str | None:
        """The commit that introduced the entry header line for <eid> | <verb>."""
        if not self.ok:
            return None
        needle = f"{eid} | {verb}"
        r = self._git_in_log_dir("log", "-S", needle, "--format=%H", "--", self.log_name)
        if r.returncode != 0:
            return None
        lines = [l for l in r.stdout.splitlines() if l.strip()]
        return lines[0] if lines else None


def _normalize_push_event(ev: dict) -> dict | None:
    """Accept either the raw Events-API shape or the compact fixture shape.
    The `after` field is OPTIONAL: the events feed elides it for some push
    events (observed live on 2026-10-20), and a prior-tip citation relies
    only on ref + before."""
    if not isinstance(ev, dict):
        return None
    payload = ev.get("payload") or {}
    ref = ev.get("ref") or payload.get("ref")
    before = ev.get("before") or payload.get("before")
    after = ev.get("after") or payload.get("after")
    if not (isinstance(ref, str) and isinstance(before, str)):
        return None
    if after is not None and not isinstance(after, str):
        after = None
    return {"id": str(ev.get("id", "")), "ref": ref, "before": before,
            "after": after, "created_at": ev.get("created_at", "")}


class RunVerifier:
    """Verify cited Actions runs and push events against the manifest, the
    GitHub Actions API, and the repository Events API.

    Fail-closed policy (Round-21 direction 5.6 + Round-22 direction 3): with a
    manifest configured (the CI gate always configures one), a citation
    outside the manifest window MUST verify through the API or the entry
    FAILS -- a manifest omission is never an approval. With no manifest and
    no API configured (local tool mode), citations are disclosed as
    unverified rather than trusted; the gate never runs in that mode.

    Round-22: run records retain and validate `event` and `head_branch`
    alongside head_sha and conclusion, so push-event provenance can be
    distinguished from CI success; the Events API supplies genuine
    PushEvent before/after state for the pre-push binding.
    """

    def __init__(self, manifest: dict | None, repo: str | None, token_env: str,
                 push_events_fixture: list | None = None):
        self.manifest = manifest or {}
        self.repo = repo
        self.token = os.environ.get(token_env, "") if token_env else ""
        self._cache: dict[str, dict | None] = {}
        self.has_api = bool(repo and self.token)
        self.events_fixture = ([_normalize_push_event(ev)
                                for ev in (push_events_fixture or [])]
                               if push_events_fixture is not None else None)
        self.events_fixture = [ev for ev in self.events_fixture if ev] \
            if self.events_fixture else self.events_fixture
        self._events_cache: list[dict] | None = None
        self._events_status: str | None = None

    def _api_run(self, rid: str) -> dict | None:
        req = urllib.request.Request(
            f"{API}/repos/{self.repo}/actions/runs/{rid}",
            headers={"Accept": "application/vnd.github+json",
                     **({"Authorization": "Bearer " + self.token} if self.token else {})})
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                d = json.loads(r.read() or b"{}")
                return {"head_sha": d.get("head_sha"),
                        "head_branch": d.get("head_branch"),
                        "event": d.get("event"),
                        "conclusion": d.get("conclusion")}
        except (urllib.error.HTTPError, urllib.error.URLError, OSError, ValueError):
            return None

    def verify(self, rid: str) -> tuple[dict | None, str]:
        """Returns (record, classification). classification in
        {manifest, api, absent-from-manifest-no-api, api-unavailable, unverified}."""
        if rid in self._cache:
            rec = self._cache[rid]
        else:
            rec = None
            src = None
            if self.manifest:
                rec = self.manifest.get("runs", {}).get(rid)
                if rec is not None:
                    src = "manifest"
                elif self.has_api:
                    rec = self._api_run(rid)
                    src = "api" if rec else "api-unavailable"
                else:
                    src = "absent-from-manifest-no-api"
            elif self.has_api:
                rec = self._api_run(rid)
                src = "api" if rec else "api-unavailable"
            else:
                src = "no-source-configured"
            rec = dict(rec) if rec else None
            if rec is not None:
                rec["_source"] = src
            self._cache[rid] = rec
        if rec is not None:
            return rec, rec["_source"]
        if self.manifest and not self.has_api:
            return None, "absent-from-manifest-no-api"
        if self.manifest or self.has_api:
            return None, "unverified"
        return None, "no-source-configured"

    def api_commit_exists(self, sha: str) -> bool | None:
        """commits API existence check; None when unavailable."""
        if not (self.repo and self.token):
            return None
        req = urllib.request.Request(
            f"{API}/repos/{self.repo}/commits/{sha}",
            headers={"Accept": "application/vnd.github+json",
                     "Authorization": "Bearer " + self.token})
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                return r.status == 200
        except urllib.error.HTTPError as e:
            if e.code in (404, 422):
                return False
            return None
        except (urllib.error.URLError, OSError):
            return None

    # ------------------------- Events API (Round-22 direction 4)


    def api_is_ancestor(self, base: str, head: str) -> bool | None:
        """Is <base> an ancestor of <head> per the GitHub compare API?

        Used for a PRUNED branch whose archive tag holds the branch's final
        tip: an earlier confirmation of that branch is proven by ancestry
        inside the archived history. Fail-closed: any doubt is None.
        """
        if not self.token or not self.repo:
            return None
        try:
            req = urllib.request.Request(
                f"{API}/repos/{self.repo}/compare/{base}...{head}",
                headers={"Accept": "application/vnd.github+json",
                         **({"Authorization": "Bearer " + self.token} if self.token else {})})
            with urllib.request.urlopen(req, timeout=30) as r:
                d = json.loads(r.read() or b"{}")
                return str(d.get("status")) in ("ahead", "identical")
        except (urllib.error.URLError, urllib.error.HTTPError, OSError, ValueError):
            return None

    def events_available(self) -> bool:
        return self.events_fixture is not None or self.has_api

    def _fetch_push_events(self) -> None:
        if self._events_cache is not None:
            return
        if self.events_fixture is not None:
            self._events_cache = list(self.events_fixture)
            self._events_status = "fixture"
            return
        if not self.has_api:
            self._events_cache = []
            self._events_status = "no-source"
            return
        evs: list[dict] = []
        import time
        # the Events API serves at most 3 pages / 300 events per repo; a 422
        # carrying "pagination is limited" AFTER a successful page is the
        # documented END OF DATA for this resource, not an outage (observed
        # live 2026-10-20: page 4 returns exactly that). Genuine failures on
        # any page still fail closed once the bounded retries are exhausted.
        page = 1
        while page <= 3:
            req = urllib.request.Request(
                f"{API}/repos/{self.repo}/events?per_page=100&page={page}",
                headers={"Accept": "application/vnd.github+json",
                         "Authorization": "Bearer " + self.token})
            data = None
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(req, timeout=30) as r:
                        data = json.loads(r.read() or b"[]")
                    break
                except urllib.error.HTTPError as e:
                    body = b""
                    try:
                        body = e.read()
                    except Exception:  # noqa: BLE001
                        pass
                    if e.code == 422 and b"pagination is limited" in body \
                            and page > 1:
                        data = []
                        break
                    if attempt < 2:
                        time.sleep(2 * (attempt + 1))
                except (urllib.error.URLError, OSError, ValueError):
                    if attempt < 2:
                        time.sleep(2 * (attempt + 1))
            if data is None:
                self._events_cache = []
                self._events_status = "unavailable"
                return
            got = data if isinstance(data, list) else []
            for ev in got:
                if isinstance(ev, dict) and ev.get("type") == "PushEvent":
                    nev = _normalize_push_event(ev)
                    if nev:
                        evs.append(nev)
            if len(got) < 100:
                break
            page += 1
        self._events_cache = evs
        self._events_status = "ok"

    def push_event_for(self, sha: str) -> tuple[dict | None, str]:
        """The latest PushEvent whose payload.after == sha.

        Returns (event, status): status in {found, not-found, no-source,
        unavailable}. 'unavailable' is a fetch failure (fail closed when the
        binding is required); 'not-found' means the event is absent from the
        visible window (rolled off, or the SHA was never pushed) -- never
        treated as proof."""
        self._fetch_push_events()
        matches = [e for e in self._events_cache if e.get("after") == sha]
        if matches:
            return matches[-1], "found"
        if self._events_status in ("fixture", "ok"):
            return None, "not-found"
        return None, self._events_status or "no-source"

    def prior_tip_event(self, event_id: str, sha: str) -> tuple[dict | None, str]:
        """Round-22 prior-tip re-verification: the retained PushEvent with
        the recorded id whose `before` == sha (a later push on the ref whose
        before field records the ref's prior tip). Returns (event, status)
        with the same status vocabulary as push_event_for."""
        self._fetch_push_events()
        for e in self._events_cache:
            if e.get("id") == event_id and e.get("before") == sha:
                return e, "found"
        if self._events_status in ("fixture", "ok"):
            return None, "not-found"
        return None, self._events_status or "no-source"

    def push_run_for(self, branch: str, sha: str) -> str | None:
        """A verified Actions run with event=push, head_branch == branch and
        head_sha == sha (the manifest first, else the Actions API queried by
        head_sha). A run of ANY conclusion qualifies: a failed run still
        proves the push event occurred (direction 3) and is never counted as
        successful gate evidence."""
        for rid, rec in (self.manifest.get("runs", {}) or {}).items():
            if (rec.get("event") == "push" and rec.get("head_branch") == branch
                    and rec.get("head_sha") == sha):
                return rid
        if not self.has_api:
            return None
        req = urllib.request.Request(
            f"{API}/repos/{self.repo}/actions/runs?head_sha={sha}&per_page=100",
            headers={"Accept": "application/vnd.github+json",
                     "Authorization": "Bearer " + self.token})
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                data = json.loads(r.read() or b"{}")
        except (urllib.error.HTTPError, urllib.error.URLError, OSError, ValueError):
            return None
        for run in data.get("workflow_runs", []) or []:
            if (run.get("event") == "push"
                    and run.get("head_branch") == branch
                    and run.get("head_sha") == sha):
                return str(run.get("id"))
        return None

    def transition_evidence(self, branch: str, declared: str, candidate: str):
        """Round-22 binding evidence for a ref transition, ADAPTED to this
        repository's live feed realities (verified 2026-10-20): the events
        feed ELIDES PushEvent.after for this repo (only `before` is
        reliable), and ref updates made through the git-data API generate
        NO PushEvent at all (their pushes are proven only by push-event
        Actions runs). Preference order, every step independently verified:

          1. 'push-event': a retained PushEvent on refs/heads/<branch> with
             after == candidate (when the feed provides after); its before
             must equal the declared pre-push tip (all-zero for ABSENT),
             else a mismatch reason is returned (race / false ABSENT);
          2. 'push-event-before': a retained PushEvent on the ref whose
             before == the declared pre-push tip (after elided), PLUS a
             verified push-event run on the exact candidate tying the
             candidate to the branch push;
          3. 'push-run': a verified push-EVENT Actions run on the exact
             candidate (event=push, head_branch, head_sha) -- GitHub's own
             record that the push event occurred on this branch and head;
          4. None: no push-event evidence exists (fail closed).

        Returns (kind, obj, mismatch): kind in {'push-event',
        'push-event-before', 'push-run', None}; obj is the event dict or the
        run id; mismatch is a failure reason for 'push-event' mismatches."""
        ref = f"refs/heads/{branch}"
        self._fetch_push_events()
        on_ref = [e for e in self._events_cache if e.get("ref") == ref]
        # 1. after-match (strongest: proves both ends of the transition)
        for e in reversed(on_ref):
            if e.get("after") == candidate:
                before = e.get("before")
                if declared == "ABSENT":
                    if before != ZERO_SHA:
                        return ("push-event", e,
                                f"false ABSENT: the PushEvent records a prior "
                                f"tip {str(before)[:12]}… (the ref existed "
                                f"before this push)")
                elif isinstance(before, str) and before != declared:
                    return ("push-event", e,
                            f"the PushEvent's before is {before[:12]}…, not the "
                            f"declared pre-push tip {declared[:12]}… (a prior "
                            f"unlogged ref update, a raced update, or an older "
                            f"ancestor was declared; the log does not describe "
                            f"the actual ref transition)")
                return ("push-event", e, None)
        # 2. before-match + push-run on the candidate (after elided)
        want_before = ZERO_SHA if declared == "ABSENT" else declared
        rid = self.push_run_for(branch, candidate)
        for e in reversed(on_ref):
            if e.get("before") == want_before and rid is not None:
                return ("push-event-before", e, None)
        # 3. push-event run on the exact candidate alone
        if rid is not None:
            return ("push-run", rid, None)
        return (None, None, None)

    # ------------------------- PR / reviews (Round-21/22 merge verification)

    def api_pull(self, number: str) -> dict | None:
        if not (self.repo and self.token):
            return None
        req = urllib.request.Request(
            f"{API}/repos/{self.repo}/pulls/{number}",
            headers={"Accept": "application/vnd.github+json",
                     "Authorization": "Bearer " + self.token})
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read() or b"{}")
        except (urllib.error.HTTPError, urllib.error.URLError, OSError, ValueError):
            return None

    def api_pull_reviews(self, number: str) -> list | None:
        """All reviews, paginated (Round-22: a single 100-entry page must not
        silently truncate the review history)."""
        if not (self.repo and self.token):
            return None
        out: list = []
        page = 1
        while page <= 5:
            req = urllib.request.Request(
                f"{API}/repos/{self.repo}/pulls/{number}/reviews"
                f"?per_page=100&page={page}",
                headers={"Accept": "application/vnd.github+json",
                         "Authorization": "Bearer " + self.token})
            try:
                with urllib.request.urlopen(req, timeout=30) as r:
                    data = json.loads(r.read() or b"[]")
            except (urllib.error.HTTPError, urllib.error.URLError, OSError, ValueError):
                return None
            got = data if isinstance(data, list) else []
            out.extend(got)
            if len(got) < 100:
                break
            page += 1
        return out


def _valid_approvals(reviews: list, author: str | None, head: str | None,
                     before_iso: str | None) -> tuple[list[dict], list[str]]:
    """Round-22 direction 5: an approval counts only when it is from a
    non-author reviewer, applies to the EXACT declared head (commit_id), was
    submitted BEFORE the merge (when post-merge), and is the reviewer's
    LATEST review state (later changes/dismissals supersede earlier
    approvals). Returns (valid_approvals, rejection_reasons_sample)."""
    latest: dict[str, dict] = {}
    for rv in reviews:
        if not isinstance(rv, dict):
            continue
        user = (rv.get("user") or {}).get("login")
        if not user:
            continue
        prev = latest.get(user)
        try:
            ts = parse_ts(rv.get("submitted_at") or "1970-01-01T00:00:00Z")
        except ValueError:
            continue
        if prev is None or ts >= prev["_ts"]:
            rv2 = dict(rv)
            rv2["_ts"] = ts
            latest[user] = rv2
    valid: list[dict] = []
    rejected: list[str] = []
    for user, rv in sorted(latest.items()):
        if user == author:
            rejected.append(f"{user}: author")
            continue
        if rv.get("state") != "APPROVED":
            rejected.append(f"{user}: latest state {rv.get('state')!r}")
            continue
        cid = rv.get("commit_id")
        if head and cid != head:
            rejected.append(f"{user}: approval bound to commit "
                            f"{(cid or 'none')[:12]}…, not the exact head")
            continue
        if before_iso:
            try:
                if rv["_ts"] > parse_ts(before_iso):
                    rejected.append(f"{user}: submitted after the merge")
                    continue
            except ValueError:
                rejected.append(f"{user}: unparseable merge timestamp")
                continue
        valid.append(rv)
    return valid, rejected


# ---------------------------------------------------------------- remote verification

def check_remote_verification(path: str, state: dict, remote: str,
                              branch_under_test: str | None, expect_tip: str | None,
                              runs: RunVerifier | None, logrepo: LogRepo,
                              pr_fixture: dict | None, reviews_fixture: list | None,
                              ) -> tuple[list[str], list[str]]:
    """Round-20/21/22 layers: live-ref comparison, evidence support, claim-
    separated historical proof (commit existence vs push provenance vs CI
    success), terminal-intent binding with push-event before/after state, and
    the protected-branch merge transition."""
    findings: list[str] = []
    info: list[str] = []
    entries = state["entries"]
    live_cache: dict[str, tuple[str | None, str | None]] = {}

    def live(branch: str):
        if branch not in live_cache:
            live_cache[branch] = ls_remote_tip(remote, branch)
        return live_cache[branch]

    archive_cache: dict[str, tuple[str | None, str | None, str | None]] = {}

    def archive(branch: str):
        if branch not in archive_cache:
            archive_cache[branch] = ls_remote_archive(remote, branch)
        return archive_cache[branch]

    def entry_runs(e) -> tuple[str, list[str]]:
        evm = re.search(r'evidence="([^"]*)"', e["subject"])
        text = (evm.group(1) if evm else "") + "\n" + "\n".join(e["body"])
        return text, sorted(set(RUN_ID.findall(text)))

    def verify_citations(e, sha: str, branch: str) -> tuple[list[str], list[str]]:
        """Verify every run citation in the entry (field + body). Returns
        (ci_success_rids, push_provenance_rids). Round-22: each cited run's
        EVENT TYPE and HEAD_BRANCH are validated against the claim -- a
        push-event run on a different branch contradicts the entry and FAILs;
        a failed/cancelled run is disclosed and never counts as success."""
        _, rids = entry_runs(e)
        supporting: list[str] = []
        provenance: list[str] = []
        for rid in rids:
            rec, cls = runs.verify(rid) if runs else (None, "no-source-configured")
            if rec is None:
                if cls == "no-source-configured":
                    info.append(f"NOTE: {e['eid']}: cited run {rid} is unverified "
                                f"(no runs-verification source configured; the CI gate "
                                f"always configures one -- local tool mode)")
                else:
                    findings.append(
                        f"FAIL: {path}:{e['n']}: {e['eid']} cites run {rid} which could not "
                        f"be verified ({cls}) -- an unverifiable citation fails the gate "
                        f"(a manifest omission is never an approval)")
                continue
            if rec.get("head_sha") != sha:
                findings.append(f"FAIL: {path}:{e['n']}: {e['eid']} cites run {rid} as "
                                f"evidence for {sha[:12]}… but that run's head_sha is "
                                f"{str(rec.get('head_sha'))[:12]}… (evidence does not "
                                f"support the recorded SHA)")
                continue
            ev_type = rec.get("event")
            hb = rec.get("head_branch")
            if ev_type == "push" and hb is not None and hb != branch:
                findings.append(f"FAIL: {path}:{e['n']}: {e['eid']} cites run {rid} as "
                                f"push evidence for branch {branch!r} but that run is a "
                                f"push event on branch {hb!r} -- a push event on another "
                                f"branch does not support this entry (Round-22 event/"
                                f"head_branch validation)")
                continue
            if ev_type == "push" and hb == branch:
                provenance.append(rid)
            if rec.get("conclusion") == "success":
                supporting.append(rid)
            else:
                info.append(f"NOTE: {e['eid']}: cited run {rid} ran on the exact head "
                            f"(event={ev_type!r}) but concluded {rec.get('conclusion')!r} -- "
                            f"CI-success is NOT claimed from it; recorded honestly; a failed "
                            f"or cancelled run can never be represented as successful gate "
                            f"evidence (it may still evidence the push event itself)")
        return supporting, provenance

    # ---- pass 1: CONFIRM/RECONCILE entries
    current_for: dict[str, str] = {}
    for e in entries:
        if not e["m"] or e["verb"] not in CONFIRM_VERBS:
            continue
        subject, eid = e["subject"], e["eid"]
        branch = field(subject, "branch")
        sha = field(subject, "remote_sha")
        if not (branch and sha and re.fullmatch(r"[0-9a-f]{40}", sha)):
            continue  # field grammar already failed in schema check
        evidence = field(subject, "evidence") or ""
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
        # 2. claimed branch must be named IN THE EVIDENCE FIELD itself
        if branch not in evidence:
            findings.append(f"FAIL: {path}:{e['n']}: {eid} claims branch {branch!r} but the "
                            f"evidence does not mention that branch (mismatched branch)")
        # 3. run citations (evidence field AND body) with event/head_branch
        #    validation; CI-success and push-provenance are tracked separately
        supporting, provenance = verify_citations(e, sha, branch) if runs else ([], [])
        # 4. live-ref comparison
        tip, err = live(branch)
        if err:
            # A reviewed `push-not-proven` reclassification is disclosed, never
            # counted as evidence: the branch push is NOT PROVEN (the commit
            # may still exist). This is the same treatment the off-lineage
            # HISTORICAL path gives, applied when the branch ref itself is
            # gone and the archive tag cannot preserve the record.
            if eid in state.get("reclassified", set()):
                rc_eid = next((rc["eid"] for rc in entries
                               if rc["m"] and rc["verb"] == "push-not-proven"
                               and field(rc["subject"], "reclassifies") == eid), eid)
                info.append(f"NOT-PROVEN: {eid} confirmed {branch} @ {sha[:12]}… is "
                            f"RECLASSIFIED by {rc_eid}: the branch push is NOT PROVEN "
                            f"(commit existence may be proven; push provenance is not) -- "
                            f"the branch ref is gone and its archive tag cannot preserve this "
                            f"record; disclosed, never counted as evidence for any declaration")
                continue
            # D-V3/D-V4: a pruned branch keeps its content in an ANNOTATED
            # archive tag at the recorded SHA; that tag is remote evidence.
            # An archive tag at a DIFFERENT SHA does not preserve the
            # confirmation -- both paths stay fail-closed.
            archived, tag_ref, aerr = archive(branch)
            if archived is None:
                findings.append(f"FAIL: {path}:{e['n']}: {eid} could not be verified against the "
                                f"live remote ref: {err} (fail closed)")
                continue
            if archived == sha:
                info.append(f"ARCHIVED: {eid} confirmed {branch} @ {sha[:12]}…; the branch ref "
                            f"was pruned and its content is preserved by {tag_ref} at the "
                            f"identical SHA (the archive tag is the remote evidence; the ref is "
                            f"gone, the commit is not); CI-success support: {len(supporting)} "
                            f"verified successful run(s)")
                continue
            # The archive tag holds the branch's FINAL tip. An older
            # confirmation describes an earlier push of that branch, so it is
            # proven by ANCESTRY inside the archived history - the same rule
            # the live-ref path uses, applied to the archive tag.
            anc = logrepo.is_ancestor(sha, archived)
            if anc is not True and runs is not None:
                anc = runs.api_is_ancestor(sha, archived)
            if anc is True:
                info.append(f"ARCHIVED-HISTORICAL: {eid} confirmed {branch} @ {sha[:12]}…; the "
                            f"branch ref was pruned and its content is preserved by {tag_ref} at "
                            f"{archived[:12]}…; the recorded SHA is an ancestor of the archived "
                            f"commit (the confirmation stays valid for the push it describes and "
                            f"does NOT confirm the archived tip); CI-success support: "
                            f"{len(supporting)} verified successful run(s)")
                continue
            findings.append(
                f"FAIL: {path}:{e['n']}: {eid} could not be verified against the live remote "
                f"ref: {err}; the archive tag {tag_ref} points at {archived[:12]}…, which is not "
                f"the recorded {sha[:12]}… and does not contain it in its history (an archive tag "
                f"at another SHA does not preserve this confirmation -- fail closed)")
            continue
        if tip == sha:
            current_for[branch] = sha
            info.append(f"CURRENT: {eid} confirms the live tip of {branch} ({sha[:12]}…); "
                        f"CI-success support: {len(supporting)} verified successful run(s)")
            continue
        # 5. HISTORICAL -- Round-22 direction 3: commit existence, push
        #    provenance, and CI success are SEPARATE claims
        anc = logrepo.is_ancestor(sha, tip)
        if anc is True:
            info.append(f"HISTORICAL: {eid} confirmed {branch} @ {sha[:12]}…; the live tip "
                        f"is {tip[:12]}… — SUPPORTED-BY-ANCESTRY (the recorded SHA is on "
                        f"the branch's current lineage; the confirmation is valid for the "
                        f"push it describes and does NOT confirm the newer tip); "
                        f"CI-success support: {len(supporting)} verified successful run(s)")
            continue
        # off-lineage: an explicit reclassification discloses NOT PROVEN
        if eid in state.get("reclassified", set()):
            rc_eid = next((rc["eid"] for rc in entries
                           if rc["m"] and rc["verb"] == "push-not-proven"
                           and field(rc["subject"], "reclassifies") == eid), eid)
            info.append(f"NOT-PROVEN: {eid} confirmed {branch} @ {sha[:12]}… is "
                        f"RECLASSIFIED by {rc_eid}: the branch push is NOT PROVEN "
                        f"(commit existence may be proven; push provenance is not) -- "
                        f"disclosed, never counted as evidence for any declaration")
            continue
        # off-lineage: a reviewed evidence-back record is required
        backs = [b for b in entries if b["m"] and b["verb"] == "evidence-back"
                 and field(b["subject"], "backs") == eid
                 and field(b["subject"], "branch") == branch
                 and field(b["subject"], "remote_sha") == sha]
        if not backs:
            findings.append(
                f"FAIL: {path}:{e['n']}: {eid} is classified HISTORICAL ({branch} @ "
                f"{sha[:12]}…, live tip {tip[:12]}…), the SHA is not on the branch's "
                f"current lineage, and no reviewed evidence-back record supports it -- "
                f"a different SHA alone does not prove the historical push was real, and "
                f"commit existence is not push provenance (Round-22 direction 3)")
            continue
        back = backs[-1]
        # claim (i): commit existence -- object store or commits API
        obj = logrepo.object_exists(sha)
        api_exists = runs.api_commit_exists(sha) if runs else None
        if not obj and api_exists is not True:
            findings.append(
                f"FAIL: {path}:{e['n']}: the historical SHA {sha[:12]}… backed by "
                f"{back['eid']} is not present in the local object store and could not "
                f"be verified as a real remote commit (fail closed)")
            continue
        # the back's own citations verify (fail-closed policy shared above)
        back_support, back_prov = verify_citations(back, sha, branch) if runs else ([], [])
        # claim (ii): push provenance -- independent branch-push proof
        provenance_src = None
        if back_prov:
            provenance_src = (f"push-event run {back_prov[-1]} "
                              f"(event=push, head_branch={branch}, head_sha={sha[:12]}…)")
        else:
            ev, evs = runs.push_event_for(sha) if runs else (None, "no-source")
            rr = field(back["subject"], "review_ref")
            if evs == "unavailable":
                info.append(f"NOTE: {back['eid']}: the Events API is unavailable; "
                            f"live PushEvent verification is skipped (never treated "
                            f"as proof)")
            elif ev is not None and ev.get("ref") == f"refs/heads/{branch}" \
                    and ev.get("after") == sha:
                if rr:
                    provenance_src = (f"live Events-API PushEvent {ev.get('id')} "
                                      f"(ref={ev.get('ref')}, before={ev['before'][:12]}…, "
                                      f"after={ev['after'][:12]}…), adjudicated by "
                                      f"review_ref={rr}")
                else:
                    findings.append(
                        f"FAIL: {path}:{back['n']}: {back['eid']} is backed by a live "
                        f"PushEvent but carries no review_ref -- Events-API evidence "
                        f"is transient (retention window) and must be recorded and "
                        f"adjudicated to count durably (fail closed)")
            elif ev is not None:
                info.append(f"NOTE: {back['eid']}: a PushEvent exists for the SHA but on "
                            f"ref {ev.get('ref')!r}, not refs/heads/{branch} (not proof "
                            f"for this branch)")
            # recorded push_event citation + adjudication
            if provenance_src is None:
                ev_text = ((field(back["subject"], "evidence") or "")
                           + "\n" + "\n".join(back["body"]))
                m = PUSH_EVENT_CITE.search(ev_text)
                if m and m.group(4) == sha:
                    if m.group(2) != f"refs/heads/{branch}":
                        findings.append(
                            f"FAIL: {path}:{back['n']}: {back['eid']} records a push_event "
                            f"citation on ref {m.group(2)} which is not the claimed branch "
                            f"ref (mismatched ref-transition evidence)")
                    elif not rr:
                        findings.append(
                            f"FAIL: {path}:{back['n']}: {back['eid']} records push_event "
                            f"evidence but carries no review_ref -- recorded ref-transition "
                            f"evidence must be adjudicated against the canonical record "
                            f"(fail closed)")
                    else:
                        provenance_src = (f"recorded PushEvent {m.group(1)} "
                                          f"(ref={m.group(2)}, before={m.group(3)[:12]}…, "
                                          f"after={m.group(4)[:12]}…), adjudicated by "
                                          f"review_ref={rr}")
                        lev, _ls = runs.push_event_for(sha) if runs else (None, "no-source")
                        if lev is not None and lev.get("id") == m.group(1):
                            info.append(f"NOTE: {back['eid']}: the recorded PushEvent "
                                        f"{m.group(1)} is still observable live "
                                        f"(re-verified)")
                        elif runs and runs.events_available():
                            info.append(f"NOTE: {back['eid']}: the recorded PushEvent "
                                        f"{m.group(1)} is no longer observable in the "
                                        f"live events window (retention limit); the "
                                        f"adjudicated record carries the proof")
                # Round-22 prior-tip citation: a LATER push on the same ref
                # whose before field records the backed SHA as the ref's
                # prior tip (used when the push's own event rolled off and
                # the API elides its after)
                if provenance_src is None:
                    m2 = PRIOR_TIP_CITE.search(ev_text)
                    if m2 and m2.group(3) == sha:
                        if m2.group(2) != f"refs/heads/{branch}":
                            findings.append(
                                f"FAIL: {path}:{back['n']}: {back['eid']} records a "
                                f"prior-tip push_event citation on ref {m2.group(2)} "
                                f"which is not the claimed branch ref (mismatched "
                                f"ref-transition evidence)")
                        elif not rr:
                            findings.append(
                                f"FAIL: {path}:{back['n']}: {back['eid']} records a "
                                f"prior-tip push_event citation but carries no "
                                f"review_ref -- recorded ref-transition evidence must "
                                f"be adjudicated against the canonical record "
                                f"(fail closed)")
                        else:
                            provenance_src = (f"recorded PushEvent {m2.group(1)} "
                                              f"(ref={m2.group(2)}, "
                                              f"prior_tip={m2.group(3)[:12]}…): a later "
                                              f"push on the same ref records the branch's "
                                              f"prior tip as the backed SHA -- the ref "
                                              f"provably pointed there (refs move only "
                                              f"by pushes), adjudicated by "
                                              f"review_ref={rr}")
                            lev2, _ls2 = (runs.prior_tip_event(m2.group(1), sha)
                                          if runs else (None, "no-source"))
                            if lev2 is not None:
                                info.append(f"NOTE: {back['eid']}: the recorded "
                                            f"prior-tip PushEvent {m2.group(1)} is "
                                            f"still observable live (re-verified)")
                            elif runs and runs.events_available():
                                info.append(f"NOTE: {back['eid']}: the recorded "
                                            f"prior-tip PushEvent {m2.group(1)} is no "
                                            f"longer observable in the live events "
                                            f"window (retention limit); the "
                                            f"adjudicated record carries the proof")
        if provenance_src:
            ci_note = (f"CI-success: {len(back_support)} verified successful run(s)"
                       if back_support else
                       "CI-success: none of the cited runs succeeded (failures disclosed "
                       "separately; a failed run may evidence the push event but is never "
                       "successful gate evidence)")
            info.append(f"HISTORICAL: {eid} confirmed {branch} @ {sha[:12]}… -- "
                        f"SUPPORTED-BY-EVIDENCE-BACK {back['eid']}: PUSH-PROVEN by "
                        f"{provenance_src}; commit-existence proof "
                        f"{'object-store' if obj else 'commits-API'}; {ci_note}")
        else:
            # local tool mode (no runs/events source configured): the provenance
            # gap is DISCLOSED, not silently passed -- the CI gate always
            # configures the manifest (+ API) and fails exactly this case
            # closed (the established Round-21 no-source policy).
            no_source = (runs is None
                         or (not runs.manifest and not runs.has_api
                             and runs.events_fixture is None))
            if no_source:
                info.append(
                    f"NOTE: {eid} ({branch} @ {sha[:12]}…, off-lineage) is backed by "
                    f"{back['eid']}, but push provenance is UNVERIFIED in this mode "
                    f"(no runs-manifest, API, or events source configured; local tool "
                    f"mode discloses -- the CI gate always configures sources and "
                    f"fails this case closed)")
            else:
                findings.append(
                    f"FAIL: {path}:{e['n']}: {eid} ({branch} @ {sha[:12]}…, off-lineage) is "
                    f"backed by {back['eid']}, but the backing proves only COMMIT EXISTENCE: "
                    f"no independent branch-push proof exists (no push-event run with "
                    f"event=push + head_branch={branch}, no Events-API PushEvent, no "
                    f"adjudicated push_event citation) and no push-not-proven "
                    f"reclassification discloses the gap -- commit existence is not push "
                    f"provenance (Round-22 direction 3; fail closed)")

    # ---- pass 2: merge-intents
    merge_intents = [e for e in entries if e["m"] and e["verb"] == "merge-intent"]
    voided = state.get("voided_intents", set())

    def mi_field(e, name):
        return field(e["subject"], name)

    if branch_under_test:
        targeting = [e for e in merge_intents if mi_field(e, "target") == branch_under_test]

        def resolve_head(e):
            """expected_head=DEFERRED resolves to the commit introducing the record
            (a record cannot embed its own commit's SHA); a 40-hex value is used
            verbatim (the intent then rides an ancestor of that head)."""
            v = mi_field(e, "expected_head")
            if v == "DEFERRED":
                intro = logrepo.introducing_commit(e["eid"], "merge-intent")
                return intro, intro is not None
            return v, True

        if targeting:
            active = targeting[-1]
            # Round-22 direction 5: fix merge-intent supersession -- resolve
            # DEFERRED heads to their introducing commits before comparing
            for old in targeting[:-1]:
                if old["eid"] in voided:
                    continue
                old_head, old_ok = resolve_head(old)
                new_head, _new_ok = resolve_head(active)
                if old_ok and old_head and new_head \
                        and logrepo.is_ancestor(old_head, new_head):
                    info.append(f"NOTE: merge-intent {old['eid']} is SUPERSEDED by "
                                f"{active['eid']} (same lineage)")
                else:
                    findings.append(f"FAIL: {path}:{active['n']}: stale merge-intent "
                                    f"{old['eid']} is not on the lineage of the active one "
                                    f"and is not voided (void it with intent-reconcile)")
            # post-merge verification: the candidate must BE the verified merge
            base = mi_field(active, "expected_pre_merge_base")
            head, head_ok = resolve_head(active)
            pr_no = mi_field(active, "pr")
            src_branch = mi_field(active, "source_branch")
            merge_ok = head_ok
            if not head_ok:
                findings.append(f"FAIL: {path}:{active['n']}: the merge-intent's declared "
                                f"head (DEFERRED) could not be resolved to its "
                                f"introducing commit (fail closed)")
            if expect_tip:
                parents = logrepo.parents(expect_tip)
                if not parents or len(parents) != 3:
                    merge_ok = False
                    findings.append(
                        f"FAIL: {path}:{active['n']}: the live tip of {branch_under_test} "
                        f"is not a merge commit (a protected branch may only advance by a "
                        f"reviewed-PR merge; parents="
                        f"{parents if parents else 'unverifiable'})")
                else:
                    _m, p1, p2 = parents
                    if p1 != base:
                        merge_ok = False
                        findings.append(
                            f"FAIL: {path}:{active['n']}: the merge commit's first parent "
                            f"is {p1[:12]}…, not the declared pre-merge base "
                            f"{(base or '?')[:12]}… (raced or unrelated update)")
                    if p2 != head:
                        merge_ok = False
                        findings.append(
                            f"FAIL: {path}:{active['n']}: the merge commit's second parent "
                            f"is {p2[:12]}…, not the declared approved head "
                            f"{(head or '?')[:12]}…")
                    intro = logrepo.introducing_commit(active["eid"], "merge-intent")
                    if not intro or not logrepo.is_ancestor(intro, p2):
                        merge_ok = False
                        findings.append(
                            f"FAIL: {path}:{active['n']}: the merge-intent does not ride "
                            f"the merged PR branch history (introduced by "
                            f"{(intro[:12] + '…') if intro else 'nothing'})")
                # Round-22 direction 4: bind the merge to the trusted
                # push-event evidence on the protected branch (same adapted
                # preference chain as the terminal-intent binding)
                if runs and runs.events_available():
                    kind, obj, mism = runs.transition_evidence(
                        branch_under_test, base or "", expect_tip)
                    if kind is None:
                        merge_ok = False
                        findings.append(
                            f"FAIL: {path}:{active['n']}: no push-event evidence "
                            f"exists for the protected-branch candidate "
                            f"{expect_tip[:12]}… (no PushEvent on "
                            f"refs/heads/{branch_under_test} and no verified "
                            f"push-event Actions run on the exact head) — the "
                            f"merge transition cannot be verified (fail closed)")
                    elif kind == "push-event" and mism:
                        merge_ok = False
                        findings.append(
                            f"FAIL: {path}:{active['n']}: {mism} — the "
                            f"protected-branch transition is not the declared "
                            f"one (raced or unrelated update)")
                    elif kind == "push-event":
                        info.append(f"BOUND: the merge PushEvent verifies the "
                                    f"protected-branch transition (before="
                                    f"{str(obj.get('before'))[:12]}…, after="
                                    f"{str(obj.get('after'))[:12]}…)")
                    elif kind == "push-event-before":
                        info.append(f"BOUND: the merge PushEvent records the "
                                    f"transition from the declared pre-merge "
                                    f"base {str(obj.get('before'))[:12]}… (the "
                                    f"feed elides after; the push-event run on "
                                    f"the merge commit ties the candidate)")
                    else:  # push-run
                        info.append(f"BOUND: push-event run {obj} on the merge "
                                    f"commit (event=push, head_branch="
                                    f"{branch_under_test}) ties the candidate to "
                                    f"a genuine push event; the feed retains no "
                                    f"transition PushEvent for the protected "
                                    f"branch (disclosed)")
            else:
                merge_ok = False
                findings.append(f"FAIL: {path}:{active['n']}: merge verification requires "
                                f"the exact candidate (--expect-tip)")
            # PR + approval verification (GitHub API, or test fixtures)
            if pr_fixture is not None:
                pr = pr_fixture
            else:
                pr = runs.api_pull(pr_no) if (runs and pr_no) else None
            if pr is None:
                merge_ok = False
                findings.append(
                    f"FAIL: {path}:{active['n']}: PR {pr_no} could not be verified "
                    f"(no fixture, API unavailable) -- fail closed")
            else:
                # Round-22 direction 5: the PR API's merge_commit_sha must be
                # EXACTLY the protected-branch candidate being checked
                mcs = pr.get("merge_commit_sha")
                if expect_tip and mcs != expect_tip:
                    merge_ok = False
                    findings.append(
                        f"FAIL: {path}:{active['n']}: PR {pr_no} merge_commit_sha is "
                        f"{(mcs or 'none')[:12] if mcs else 'none'}…, not the live "
                        f"protected-branch candidate {expect_tip[:12]}… (merge-commit "
                        f"parents alone are insufficient)")
                if not pr.get("merged"):
                    merge_ok = False
                    findings.append(f"FAIL: {path}:{active['n']}: PR {pr_no} is not "
                                    f"merged (merged={pr.get('merged')!r})")
                merged_at = pr.get("merged_at")
                if not merged_at:
                    merge_ok = False
                    findings.append(f"FAIL: {path}:{active['n']}: PR {pr_no} carries no "
                                    f"merged_at timestamp (merged state unverifiable)")
                if not ((pr.get("head") or {}).get("sha") == head):
                    merge_ok = False
                    findings.append(
                        f"FAIL: {path}:{active['n']}: PR {pr_no} head "
                        f"{str((pr.get('head') or {}).get('sha'))[:12]}… does not "
                        f"match the declared head {(head or '?')[:12]}…")
                # Round-22: the PR's source branch ref must match the declared
                # source_branch (a wrong-branch merge is not the declared one)
                head_ref = (pr.get("head") or {}).get("ref")
                if src_branch and head_ref != src_branch:
                    merge_ok = False
                    findings.append(
                        f"FAIL: {path}:{active['n']}: PR {pr_no} head ref is "
                        f"{head_ref!r}, not the declared source branch {src_branch!r}")
                if not ((pr.get("base") or {}).get("ref") == branch_under_test):
                    merge_ok = False
                    findings.append(
                        f"FAIL: {path}:{active['n']}: PR {pr_no} base is "
                        f"{(pr.get('base') or {}).get('ref')!r}, not the protected "
                        f"branch {branch_under_test!r}")
            if reviews_fixture is not None:
                reviews = reviews_fixture
            else:
                reviews = runs.api_pull_reviews(pr_no) if (runs and pr_no) else None
            if reviews is None:
                merge_ok = False
                findings.append(
                    f"FAIL: {path}:{active['n']}: the reviews of PR {pr_no} could "
                    f"not be verified (no fixture, API unavailable) -- fail closed")
            else:
                author = (pr or {}).get("user", {}).get("login")
                approved, rejected = _valid_approvals(
                    reviews, author, head, (pr or {}).get("merged_at"))
                if not approved:
                    merge_ok = False
                    detail = "; ".join(rejected[:5]) if rejected else "no reviews"
                    findings.append(
                        f"FAIL: {path}:{active['n']}: PR {pr_no} has no approving "
                        f"review from a user other than the PR author that applies to "
                        f"the exact declared head and predates the merge (required "
                        f"independent approval not satisfied; review states: {detail})")
                else:
                    info.append(f"APPROVAL-VERIFIED: {len(approved)} fresh non-author "
                                f"approval(s) on the exact head "
                                f"({', '.join((r.get('user') or {}).get('login', '?') for r in approved)})")
            if merge_ok:
                info.append(f"MERGE-VERIFIED: {active['eid']} declares the live tip "
                            f"of {branch_under_test} as a reviewed-PR merge "
                            f"(base {base[:12]}… + approved head {head[:12]}…)")
                info.append(f"DECLARED: {branch_under_test} current tip is declared by "
                            f"the verified merge-intent {active['eid']}")

    # ---- pass 3: branch under test -- terminal-intent binding (directions 4)
    if branch_under_test:
        tip, err = live(branch_under_test)
        declared_by_merge = any(mi_field(e, "target") == branch_under_test
                                for e in merge_intents)
        if err:
            findings.append(f"FAIL: branch under test {branch_under_test!r} could not be "
                            f"verified: {err} (fail closed)")
        elif expect_tip and tip != expect_tip:
            findings.append(f"FAIL: live tip of {branch_under_test} is {tip[:12]}…, not the "
                            f"expected candidate {expect_tip[:12]}… — the checked-out "
                            f"candidate is not the live ref (stale or raced)")
        else:
            # pre-merge validation of the active merge-intent FOR THIS SOURCE branch
            src_mis = [e for e in merge_intents
                       if mi_field(e, "source_branch") == branch_under_test]
            if src_mis:
                active = src_mis[-1]
                # Round-22 direction 5: merge-intent supersession is resolved
                # via the introducing commits (DEFERRED heads included) on the
                # source side as well -- a prior intent on the same lineage is
                # SUPERSEDED by the active one; anything else must be voided
                for old in src_mis[:-1]:
                    if old["eid"] in voided:
                        continue
                    old_head, old_ok = resolve_head(old)
                    new_head, _new_ok = resolve_head(active)
                    if old_ok and old_head and new_head \
                            and logrepo.is_ancestor(old_head, new_head):
                        info.append(f"NOTE: merge-intent {old['eid']} is SUPERSEDED "
                                    f"by {active['eid']} (same lineage)")
                    else:
                        findings.append(f"FAIL: {path}:{active['n']}: stale "
                                        f"merge-intent {old['eid']} is not on the "
                                        f"lineage of the active one and is not "
                                        f"voided (void it with intent-reconcile)")
                head, head_ok = resolve_head(active)
                base = mi_field(active, "expected_pre_merge_base")
                tgt = mi_field(active, "target")
                if not head_ok:
                    findings.append(f"FAIL: {path}:{active['n']}: the merge-intent's "
                                    f"declared head (DEFERRED) could not be resolved")
                elif expect_tip and head != expect_tip:
                    findings.append(
                        f"FAIL: {path}:{active['n']}: the active merge-intent "
                        f"{active['eid']} declares head {head[:12]}… but the "
                        f"candidate is {expect_tip[:12]}… — redeclare the merge-intent "
                        f"on the exact final candidate before merging")
                if base and expect_tip:
                    if not logrepo.object_exists(base):
                        findings.append(f"FAIL: {path}:{active['n']}: the declared "
                                        f"pre-merge base {base[:12]}… is not a real commit")
                    elif not logrepo.is_ancestor(base, expect_tip):
                        findings.append(f"FAIL: {path}:{active['n']}: the declared "
                                        f"pre-merge base {base[:12]}… is not an ancestor "
                                        f"of the candidate")
                if tgt:
                    ltip, lerr = live(tgt)
                    if lerr:
                        findings.append(f"FAIL: {path}:{active['n']}: the merge target "
                                        f"{tgt!r} could not be verified: {lerr}")
                    elif base and ltip != base:
                        findings.append(
                            f"FAIL: {path}:{active['n']}: the declared pre-merge base "
                            f"{base[:12]}… is not the live tip of {tgt} (live: "
                            f"{ltip[:12]}…); the target moved -- rebase or redeclare "
                            f"before merging")
                info.append(f"NOTE: merge-intent {active['eid']} binds this candidate to "
                            f"a future reviewed merge into {tgt} (base "
                            f"{(base or '?')[:12]}…); the merge itself is executed only "
                            f"through the protected PR path")
            else:
                info.append(f"NOTE: no merge-intent is declared for "
                            f"{branch_under_test!r} (a merge into a protected branch "
                            f"requires one, added on the exact final candidate, before "
                            f"that branch may merge)")
            # normal declaration path (intents/confirms) -- skipped when the
            # branch is protected and declared by a verified merge-intent
            unpaired_here = [q for q in state["unpaired"] if q["branch"] == branch_under_test]
            declared = current_for.get(branch_under_test)
            if declared and not unpaired_here:
                info.append(f"DECLARED: {branch_under_test} current tip is CONFIRMED-CURRENT")
            elif unpaired_here:
                terminal = unpaired_here[-1]
                te = terminal["e"]
                if terminal["placeholder"]:
                    findings.append(f"FAIL: {path}:{te['n']}: the branch under test's tip "
                                    f"is declared only by a legacy placeholder push "
                                    f"({te['eid']}); a real declaration is required")
                else:
                    if not logrepo.ok:
                        findings.append(
                            f"FAIL: {path}:{te['n']}: cannot verify the terminal-intent "
                            f"binding -- the push log is not inside a git repository "
                            f"(fail closed)")
                    else:
                        intro = logrepo.introducing_commit(te["eid"], "push-intent")
                        if not intro:
                            findings.append(
                                f"FAIL: {path}:{te['n']}: the introducing commit of the "
                                f"terminal intent {te['eid']} could not be found in the "
                                f"candidate's history (fail closed)")
                        elif expect_tip and intro != expect_tip:
                            findings.append(
                                f"FAIL: {path}:{te['n']}: terminal intent {te['eid']} does "
                                f"not bind to the candidate: it was introduced by "
                                f"{intro[:12]}…, not the candidate "
                                f"{expect_tip[:12]}… — commits after the intent are "
                                f"undeclared (an intent cannot declare later commits "
                                f"merely by remaining unmatched)")
                        else:
                            info.append(f"BOUND: terminal intent {te['eid']} was "
                                        f"introduced by the exact candidate under test")
                        ert = field(te["subject"], "expected_remote_tip_before_push")
                        if ert is None:
                            findings.append(
                                f"FAIL: {path}:{te['n']}: terminal intent {te['eid']} "
                                f"lacks expected_remote_tip_before_push (required for "
                                f"pre-push verification)")
                        else:
                            tok = ert.split()[0] if ert.split() else ""
                            prior = [c for c in entries
                                     if c["m"] and c["verb"] in CONFIRM_VERBS
                                     and field(c["subject"], "branch") == branch_under_test
                                     and c["n"] < te["n"]]
                            newest_prior = prior[-1] if prior else None
                            newest_sha = (field(newest_prior["subject"], "remote_sha")
                                          if newest_prior else None)
                            if re.fullmatch(r"[0-9a-f]{40}", tok):
                                bad_decl = False
                                if expect_tip and not logrepo.is_ancestor(tok, expect_tip):
                                    # a DOCUMENTED REBUILD makes a non-ancestral
                                    # pre-push tip legitimate: when the log itself
                                    # carries a reviewed evidence-back OR a
                                    # push-not-proven reclassification for this
                                    # exact branch+SHA (records that document the
                                    # old push and the rebuild), the candidate
                                    # descends from the rebuild, not from the
                                    # declared tip
                                    backs_tok = [b for b in entries
                                                 if b["m"] and b["verb"] == "evidence-back"
                                                 and field(b["subject"], "branch") == branch_under_test
                                                 and field(b["subject"], "remote_sha") == tok]
                                    reclass_tok = [rc for rc in entries
                                                   if rc["m"] and rc["verb"] == "push-not-proven"
                                                   and field(rc["subject"], "branch") == branch_under_test
                                                   and field(rc["subject"], "remote_sha") == tok]
                                    if not backs_tok and not reclass_tok:
                                        bad_decl = True
                                        findings.append(
                                            f"FAIL: {path}:{te['n']}: declared pre-push tip "
                                            f"{tok[:12]}… is not an ancestor of the candidate "
                                            f"(the push was not a fast-forward from the "
                                            f"declared state)")
                                if newest_sha != tok:
                                    bad_decl = True
                                    findings.append(
                                        f"FAIL: {path}:{te['n']}: declared pre-push tip "
                                        f"{tok[:12]}… does not match the newest recorded "
                                        f"branch state {(newest_sha or 'ABSENT')[:12]}… "
                                        f"({newest_prior['eid'] if newest_prior else 'no prior record'}) — "
                                        f"an intervening ref update or a stale declaration; "
                                        f"reconcile before pushing")
                                if not bad_decl:
                                    info.append(f"BOUND: the declared pre-push tip "
                                                f"{tok[:12]}… matches the newest recorded "
                                                f"branch state and is an ancestor of the "
                                                f"candidate (no intervening ref update)")
                            elif tok == "ABSENT":
                                if newest_prior is not None:
                                    findings.append(
                                        f"FAIL: {path}:{te['n']}: terminal intent declares "
                                        f"expected_remote_tip_before_push=ABSENT but the "
                                        f"log already records branch state "
                                        f"{newest_sha[:12]}… ({newest_prior['eid']})")
                                else:
                                    info.append("BOUND: first push of this ref name "
                                                "(no prior branch record contradicts "
                                                "ABSENT)")
                            else:
                                findings.append(
                                    f"FAIL: {path}:{te['n']}: expected_remote_tip_before_push "
                                    f"must be a 40-hex SHA or ABSENT (schema violation)")
                        # Round-22 direction 4: bind the declared pre-push
                        # state to the TRUSTED push-event evidence, adapted
                        # to this repository's feed (after is elided;
                        # API-driven ref updates emit no PushEvent -- the
                        # push-event run on the candidate is the binding
                        # evidence in that class; see transition_evidence)
                        if expect_tip and runs and runs.events_available():
                            decl2 = (field(te["subject"],
                                           "expected_remote_tip_before_push") or "")
                            decl_tok = decl2.split()[0] if decl2.split() else ""
                            kind, obj, mism = runs.transition_evidence(
                                branch_under_test, decl_tok, expect_tip)
                            if kind is None:
                                findings.append(
                                    f"FAIL: {path}:{te['n']}: no push-event "
                                    f"evidence exists for the candidate "
                                    f"{expect_tip[:12]}… (no PushEvent on "
                                    f"refs/heads/{branch_under_test} matches the "
                                    f"declared pre-push state and no verified "
                                    f"push-event Actions run exists on the exact "
                                    f"head) — fail closed")
                            elif kind == "push-event" and mism:
                                findings.append(
                                    f"FAIL: {path}:{te['n']}: {mism}")
                            elif kind == "push-event":
                                info.append(
                                    f"BOUND: the PushEvent verifies the ref "
                                    f"transition (before={str(obj.get('before'))[:12]}…, "
                                    f"after={str(obj.get('after'))[:12]}…)")
                            elif kind == "push-event-before":
                                info.append(
                                    f"BOUND: the PushEvent records the transition "
                                    f"from the declared pre-push tip "
                                    f"{str(obj.get('before'))[:12]}… on "
                                    f"refs/heads/{branch_under_test} (the feed "
                                    f"elides after; the push-event run on the "
                                    f"exact candidate ties the candidate to the "
                                    f"branch push)")
                            else:  # push-run
                                info.append(
                                    f"BOUND: push-event run {obj} (event=push, "
                                    f"head_branch={branch_under_test}, "
                                    f"head_sha={expect_tip[:12]}…) ties the "
                                    f"candidate to a genuine push event; the "
                                    f"events feed retains no transition PushEvent "
                                    f"matching the declared pre-push state for "
                                    f"this ref (API-driven ref updates generate "
                                    f"none; disclosed) — the declared state "
                                    f"remains verified against the log's newest "
                                    f"recorded branch state and the live ref")
                    info.append(f"PENDING: {branch_under_test} current tip "
                                f"{(tip or '?')[:12]}… is declared by the unpaired "
                                f"terminal intent {te['eid']} (line {te['n']}); its "
                                f"CONFIRM is delivered by the next legitimate push "
                                f"(disclosed limitation, never backdated; labeled "
                                f"PENDING, never COMPLETE)")
            else:
                if declared_by_merge:
                    pass  # protected-branch declaration handled in pass 2
                else:
                    findings.append(f"FAIL: the live tip of {branch_under_test} "
                                    f"({(tip or '?')[:12]}…) is UNDECLARED — no unpaired "
                                    f"INTENT declares its push and no CONFIRM/RECONCILE "
                                    f"matches it (reconcile truthfully, never backdate)")

    # ---- pass 4: cross-branch legacy disclosures (unchanged from Round-20)
    for e in entries:
        if not e["m"] or e["verb"] != "push":
            continue
        branch = field(e["subject"], "branch")
        sha_m = re.search(r"sha=([0-9a-f]{40})", e["subject"])
        if not (branch and sha_m):
            continue
        sha = sha_m.group(1)
        tip, err = live(branch)
        if err:
            continue  # legacy entries: disclosed in the confirm layer if applicable
        if tip != sha and current_for.get(branch) is None \
                and not any(q["branch"] == branch for q in state["unpaired"]):
            info.append(f"NOTE: legacy push of {branch} recorded {sha[:12]}…; live tip "
                        f"{tip[:12]}… is beyond the log's last confirmation for that branch "
                        f"(disclosed; reconcile by a reviewed record if this branch re-enters "
                        f"the protocol)")
    return findings, info


# ---------------------------------------------------------- pre-merge verification

def check_pre_merge(path: str, state: dict, remote: str, branch_under_test: str,
                    expect_tip: str, runs: RunVerifier, logrepo: LogRepo,
                    pr_no: str, pr_fixture: dict | None,
                    reviews_fixture: list | None) -> tuple[list[str], list[str]]:
    """Round-22 direction 5: a FUNCTIONING pre-merge gate. On the source
    branch's PR event, verify the CURRENT PR/head/review state before any
    merge executes: the PR's exact head must be the candidate, its base must
    be the merge-intent's declared target, and a fresh non-author approval
    must apply to the exact head. A post-merge verifier is not a substitute
    for this check -- the merge must be blocked until the approval exists."""
    findings: list[str] = []
    info: list[str] = []
    entries = state["entries"]

    def pm_field(e, name):
        return field(e["subject"], name)

    merge_intents = [e for e in entries if e["m"] and e["verb"] == "merge-intent"
                     and pm_field(e, "source_branch") == branch_under_test]
    if not merge_intents:
        findings.append(f"FAIL: {path}: no merge-intent is declared for source branch "
                        f"{branch_under_test!r}; a merge into a protected branch requires "
                        f"one, bound to the exact final candidate, BEFORE merging")
        return findings, info
    active = merge_intents[-1]
    tgt = pm_field(active, "target")
    eh = pm_field(active, "expected_head")
    if eh == "DEFERRED":
        intro = logrepo.introducing_commit(active["eid"], "merge-intent")
        if not intro:
            findings.append(f"FAIL: {path}:{active['n']}: the merge-intent's declared "
                            f"head (DEFERRED) could not be resolved (fail closed)")
        elif expect_tip and intro != expect_tip:
            findings.append(f"FAIL: {path}:{active['n']}: the active merge-intent "
                            f"{active['eid']} binds head {intro[:12]}… but the candidate "
                            f"is {expect_tip[:12]}… — redeclare the merge-intent on the "
                            f"exact final candidate before merging")
    elif expect_tip and eh != expect_tip:
        findings.append(f"FAIL: {path}:{active['n']}: the active merge-intent "
                        f"{active['eid']} declares head {eh[:12]}… but the candidate is "
                        f"{expect_tip[:12]}… — redeclare before merging")

    if pr_fixture is not None:
        pr = pr_fixture
    else:
        pr = runs.api_pull(pr_no)
    if pr is None:
        findings.append(f"FAIL: {path}: PR {pr_no} could not be verified pre-merge "
                        f"(API unavailable or unauthorized) -- fail closed; the "
                        f"authenticated PR and review endpoints must be reachable "
                        f"(workflow permission pull-requests: read)")
        return findings, info
    head_sha = (pr.get("head") or {}).get("sha")
    head_ref = (pr.get("head") or {}).get("ref")
    base_ref = (pr.get("base") or {}).get("ref")
    if expect_tip and head_sha != expect_tip:
        findings.append(f"FAIL: {path}: PR {pr_no} head is "
                        f"{(head_sha or 'none')[:12] if head_sha else 'none'}…, not the "
                        f"candidate {expect_tip[:12]}… (the PR does not carry the "
                        f"checked-out candidate)")
    if head_ref != branch_under_test:
        findings.append(f"FAIL: {path}: PR {pr_no} head ref is {head_ref!r}, not the "
                        f"declared source branch {branch_under_test!r}")
    if tgt and base_ref != tgt:
        findings.append(f"FAIL: {path}: PR {pr_no} base is {base_ref!r}, not the "
                        f"merge-intent's declared target {tgt!r}")
    if pr.get("merged"):
        findings.append(f"FAIL: {path}: PR {pr_no} is already merged; this is the "
                        f"pre-merge verification (the post-merge path verifies the "
                        f"merge itself)")
        return findings, info
    if reviews_fixture is not None:
        reviews = reviews_fixture
    else:
        reviews = runs.api_pull_reviews(pr_no)
    if reviews is None:
        findings.append(f"FAIL: {path}: the reviews of PR {pr_no} could not be "
                        f"verified pre-merge (API unavailable or unauthorized) -- "
                        f"fail closed")
        return findings, info
    author = (pr.get("user") or {}).get("login")
    approved, rejected = _valid_approvals(reviews, author,
                                          head_sha or expect_tip, None)
    if not approved:
        detail = "; ".join(rejected[:5]) if rejected else "no reviews"
        findings.append(
            f"FAIL: {path}: PRE-MERGE BLOCKED -- PR {pr_no} has no current approving "
            f"review from a non-author reviewer applied to the exact head "
            f"{(head_sha or (expect_tip or '?'))[:12]}… (review states: {detail}). "
            f"The merge must wait for the required independent approval on THIS "
            f"head; this check passing is a precondition of the protected merge, "
            f"not a substitute for it")
    else:
        info.append(f"PRE-MERGE-VERIFIED: PR {pr_no} carries "
                    f"{len(approved)} current non-author approval(s) on the exact "
                    f"head; base={base_ref!r}; the merge may proceed through the "
                    f"protected PR path once every required check is green")
    return findings, info


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("logfile")
    ap.add_argument("--commit-time", help="ISO-8601 commit time carrying this log version")
    ap.add_argument("--max-drift", type=int, default=300, help="seconds of allowed drift (default 300)")
    ap.add_argument("--parent-prefix", metavar="OLDER_FILE",
                    help="require OLDER_FILE to be a byte-prefix of logfile (append-only)")
    ap.add_argument("--pushlog-schema", action="store_true",
                    help="also enforce push-log schema: unique IDs, well-formed entries with "
                         "required fields, identifier-based INTENT/CONFIRM pairing, "
                         "intent-reconcile and push-not-proven validation, and the global "
                         "unmatched-intent invariant (unpaired INTENTs are reported as "
                         "PENDING, not failures)")
    ap.add_argument("--remote-verify", action="store_true",
                    help="Round-20/21/22 layer: verify every CONFIRM/RECONCILE against the "
                         "live remote ref, prove the evidence supports the recorded SHA with "
                         "event/head_branch validation, require historical pushes to prove "
                         "the BRANCH PUSH (commit existence alone is not push provenance), "
                         "bind the terminal intent to the exact candidate and to push-event "
                         "before/after state, and verify the protected-branch merge "
                         "transition (merge_commit_sha, source ref, merged_at, fresh "
                         "non-author approval)")
    ap.add_argument("--remote", default="origin",
                    help="remote (URL or path) for ls-remote verification (default: origin)")
    ap.add_argument("--branch", default=None,
                    help="the branch under test (required with --expect-tip)")
    ap.add_argument("--expect-tip", default=None, metavar="SHA",
                    help="the expected live tip of --branch (in CI: the exact pushed SHA)")
    ap.add_argument("--runs-manifest", default=None, metavar="FILE",
                    help="JSON {\"runs\": {id: {head_sha, head_branch?, event?, conclusion}}} "
                         "for run-citation verification")
    ap.add_argument("--repo", default=None, metavar="SLUG",
                    help="owner/repo for GitHub API verification (runs outside the manifest "
                         "window, commit existence, push events, PR + review state); "
                         "requires the token env var")
    ap.add_argument("--api-token-env", default="GH_TOKEN", metavar="VAR",
                    help="environment variable carrying the API token (default GH_TOKEN; "
                         "the value is never printed or written anywhere)")
    ap.add_argument("--pr-json", default=None, metavar="FILE",
                    help="test fixture replacing the PR API response for merge verification")
    ap.add_argument("--reviews-json", default=None, metavar="FILE",
                    help="test fixture replacing the PR-reviews API response")
    ap.add_argument("--push-events-json", default=None, metavar="FILE",
                    help="test fixture replacing the Events-API PushEvent source "
                         "(a JSON list of PushEvent records)")
    ap.add_argument("--pre-merge-verify", action="store_true",
                    help="Round-22 pre-merge gate: with --branch (source) and --pr-number, "
                         "verify the CURRENT PR/head/review state before any merge "
                         "executes; fails while no fresh non-author approval applies to "
                         "the exact head")
    ap.add_argument("--pr-number", default=None, metavar="N",
                    help="the PR under test (used with --pre-merge-verify)")
    a = ap.parse_args()

    # a missing log file must fail CLOSED with a clear message, not an
    # unhandled traceback (the Round-21-followup audit hit exactly that on the
    # AT-28 direct-push scenario where the branch worktree lacked the log)
    if not os.path.exists(a.logfile):
        print(f"FAIL: evidence log {a.logfile!r} does not exist (fail closed; "
              f"a branch without its evidence log is never a pass)")
        return 1

    findings = check_sequence(a.logfile, a.commit_time, a.max_drift)
    if a.parent_prefix:
        findings += check_parent_prefix(a.logfile, a.parent_prefix)
    pending: list[str] = []
    info: list[str] = []
    logrepo = LogRepo(a.logfile)
    if a.pushlog_schema or a.remote_verify:
        schema_findings, pending, state = check_pushlog_schema(a.logfile, logrepo)
        findings += schema_findings
        if a.remote_verify:
            manifest = None
            if a.runs_manifest:
                try:
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
            push_events_fixture = None
            if a.push_events_json:
                try:
                    loaded = json.load(open(a.push_events_json, encoding="utf-8"))
                    push_events_fixture = loaded if isinstance(loaded, list) else None
                    if push_events_fixture is None:
                        print("FAIL: --push-events-json must be a JSON list of PushEvent "
                              "records; fail closed")
                        return 1
                except (OSError, ValueError) as e:
                    print(f"FAIL: --push-events-json unusable ({e}); fail closed")
                    return 1
            runs = RunVerifier(manifest, a.repo, a.api_token_env,
                               push_events_fixture=push_events_fixture)
            pr_fixture = reviews_fixture = None
            if a.pr_json:
                try:
                    pr_fixture = json.load(open(a.pr_json, encoding="utf-8"))
                except (OSError, ValueError) as e:
                    print(f"FAIL: --pr-json unusable ({e}); fail closed")
                    return 1
            if a.reviews_json:
                try:
                    reviews_fixture = json.load(open(a.reviews_json, encoding="utf-8"))
                except (OSError, ValueError) as e:
                    print(f"FAIL: --reviews-json unusable ({e}); fail closed")
                    return 1
            rv_findings, rv_info = check_remote_verification(
                a.logfile, state, remote, a.branch, a.expect_tip,
                runs, logrepo, pr_fixture, reviews_fixture)
            findings += rv_findings
            info += rv_info
            if a.pre_merge_verify:
                if not (a.branch and a.expect_tip and a.pr_number):
                    print("FAIL: --pre-merge-verify requires --branch, --expect-tip, "
                          "and --pr-number")
                    return 1
                pm_findings, pm_info = check_pre_merge(
                    a.logfile, state, remote, a.branch, a.expect_tip,
                    runs, logrepo, a.pr_number, pr_fixture, reviews_fixture)
                findings += pm_findings
                info += pm_info
    for p in pending:
        print(p)
    for i in info:
        print(i)
    print("\n".join(findings) if findings else f"PASS: {a.logfile} "
          f"{'sequence + append-only ' if a.parent_prefix else ''}"
          f"{'schema ' if (a.pushlog_schema or a.remote_verify) else ''}"
          f"{'remote-verified ' if a.remote_verify else ''}"
          f"{'pre-merge-verified ' if a.pre_merge_verify else ''}"
          f"clean")
    return 0 if not findings else 1


if __name__ == "__main__":
    sys.exit(main())
