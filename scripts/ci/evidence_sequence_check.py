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
      reviewed evidence-back record whose SHA exists (object store or the
      commits API) and whose cited runs verify. Every run citation anywhere
      in an entry (evidence field or body) must verify with the configured
      source (manifest, else the Actions API when --repo is given); failed or
      cancelled runs are disclosed and cannot prove a push; a citation that
      cannot be verified FAILS the gate (a manifest omission is never an
      approval). A textual transcript naming a branch and SHA is support
      evidence, never proof by itself.
    * PROTECTED-MAIN MERGE TRANSITION (direction 6): a push to a protected
      branch (e.g. main) is declared by a VERIFIED merge-intent: the candidate
      must be a merge commit whose first parent is the declared pre-merge
      base and whose second parent is the declared (approved) PR head, the
      merge-intent must ride the PR branch history, and the PR must be merged
      with an approving review from a user other than the PR author (GitHub
      API, or fixture files in tests). On the source branch, the active
      merge-intent must bind to the exact candidate and the declared base
      must still be the live tip of the target. Raced or unrelated updates
      to the protected branch FAIL.

Log grammar:   ## <ISO-8601 UTC> | <ENTRY-ID> | <verb> | <subject> [| RECONSTRUCTED]
Verbs: push, push-intent, push-confirm, push-reconcile (schema v3),
      merge-intent, evidence-back, intent-reconcile (schema v4).
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
RUN_ID = re.compile(r"\brun(?:_id)?\s*[=:]?\s*(\d{6,})\b", re.I)
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

def check_pushlog_schema(path: str) -> tuple[list[str], list[str], dict]:
    """Schema + identifier-based INTENT/CONFIRM pairing. Returns (findings, pending, state)."""
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

    # 3. IDENTIFIER pairing (Round-21 direction 4): a CONFIRM consumes the
    # unpaired intent/legacy-placeholder-push with the SAME id base for the
    # SAME branch. Branch/FIFO pairing alone is not accepted.
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
        match = next((q for q in queue if not q.get("consumed") and q["branch"] == branch), None)
        if match is not None:
            match["consumed"] = True
        else:
            findings.append(
                f"FAIL: {path}:{e['n']}: CONFIRM {e['eid']} has no preceding unconfirmed "
                f"INTENT (or placeholder push) with the identifier {base} for branch "
                f"{branch!r} (Round-21 identifier pairing: a confirmation must match its "
                f"own intent and must not consume another)")
    for q in [q for qs in open_by_base.values() for q in qs if not q.get("consumed")]:
        unpaired.append(q)
        kind = "legacy push" if q["placeholder"] else "push-intent"
        pending.append(f"PENDING: {path}:{q['e']['n']}: {q['e']['eid']} ({kind}) for branch "
                       f"{q['branch']} has no CONFIRM entry (two-phase protocol incomplete; "
                       f"do not imply completion)")

    # 4. evidence-back / intent-reconcile target integrity
    by_eid = {e["eid"]: e for e in entries if e["m"]}
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
    for e in entries:
        if not e["m"] or e["verb"] != "intent-reconcile":
            continue
        voids = field(e["subject"], "voids")
        target = by_eid.get(voids or "")
        if target is None or target["verb"] not in ("push-intent", "merge-intent"):
            findings.append(f"FAIL: {path}:{e['n']}: intent-reconcile {e['eid']} names "
                            f"voids={voids!r} which is not an intent entry (dangling void)")

    state = {"entries": entries, "unpaired": unpaired, "by_eid": by_eid}
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


class LogRepo:
    """The git repository containing the push log (for ancestry, object and
    introducing-commit queries). A log outside any repository is legal for
    sequence/schema checks but fails closed wherever git proof is required."""

    def __init__(self, log_path: str):
        self.ok = False
        self.toplevel = None
        self.rel = None
        try:
            r = subprocess.run(["git", "-C", os.path.dirname(os.path.abspath(log_path)),
                                "rev-parse", "--show-toplevel"],
                               capture_output=True, text=True, timeout=30)
            if r.returncode == 0 and r.stdout.strip():
                self.toplevel = r.stdout.strip()
                self.rel = os.path.relpath(os.path.abspath(log_path), self.toplevel)
                self.ok = True
        except (OSError, subprocess.TimeoutExpired):
            pass

    def _git(self, *args):
        return subprocess.run(["git", "-C", self.toplevel, *args],
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
        r = self._git("log", "-S", needle, "--format=%H", "--", self.rel)
        if r.returncode != 0:
            return None
        lines = [l for l in r.stdout.splitlines() if l.strip()]
        return lines[0] if lines else None


class RunVerifier:
    """Verify cited Actions runs against the manifest, else the GitHub API.

    Fail-closed policy (Round-21 direction 5.6): with a manifest configured
    (the CI gate always configures one), a citation outside the manifest
    window MUST verify through the API or the entry FAILS -- a manifest
    omission is never an approval. With no manifest and no API configured
    (local tool mode), citations are disclosed as unverified rather than
    trusted; the gate never runs in that mode.
    """

    def __init__(self, manifest: dict | None, repo: str | None, token_env: str):
        self.manifest = manifest or {}
        self.repo = repo
        self.token = os.environ.get(token_env, "") if token_env else ""
        self._cache: dict[str, dict | None] = {}
        self.has_api = bool(repo and self.token)

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
        if not (self.repo and self.token):
            return None
        req = urllib.request.Request(
            f"{API}/repos/{self.repo}/pulls/{number}/reviews?per_page=100",
            headers={"Accept": "application/vnd.github+json",
                     "Authorization": "Bearer " + self.token})
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read() or b"[]")
        except (urllib.error.HTTPError, urllib.error.URLError, OSError, ValueError):
            return None


# ---------------------------------------------------------------- remote verification

def check_remote_verification(path: str, state: dict, remote: str,
                              branch_under_test: str | None, expect_tip: str | None,
                              runs: RunVerifier | None, logrepo: LogRepo,
                              pr_fixture: dict | None, reviews_fixture: list | None,
                              ) -> tuple[list[str], list[str]]:
    """Round-20/21 layers: live-ref comparison, evidence support, historical
    proof, terminal-intent binding, and the protected-branch merge transition."""
    findings: list[str] = []
    info: list[str] = []
    entries = state["entries"]
    live_cache: dict[str, tuple[str | None, str | None]] = {}

    def live(branch: str):
        if branch not in live_cache:
            live_cache[branch] = ls_remote_tip(remote, branch)
        return live_cache[branch]

    def entry_runs(e) -> tuple[str, list[str]]:
        evm = re.search(r'evidence="([^"]*)"', e["subject"])
        text = (evm.group(1) if evm else "") + "\n" + "\n".join(e["body"])
        return text, sorted(set(RUN_ID.findall(text)))

    def verify_citations(e, sha: str) -> list[str]:
        """Verify every run citation in the entry (field + body). Returns the
        supporting (successful, exact-head) run ids; fails closed per policy."""
        _, rids = entry_runs(e)
        supporting: list[str] = []
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
            if rec.get("conclusion") == "success":
                supporting.append(rid)
            else:
                info.append(f"NOTE: {e['eid']}: cited run {rid} ran on the exact head "
                            f"but concluded {rec.get('conclusion')!r} -- recorded "
                            f"honestly; a failed or cancelled run cannot prove a "
                            f"successful push")
        return supporting

    # ---- pass 1: CONFIRM/RECONCILE entries
    current_for: dict[str, str] = {}
    confirms_by_eid: dict[str, dict] = {}
    for e in entries:
        if not e["m"] or e["verb"] not in CONFIRM_VERBS:
            continue
        confirms_by_eid[e["eid"]] = e
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
        # 2. claimed branch must match the evidence
        if branch not in evidence and branch not in body:
            findings.append(f"FAIL: {path}:{e['n']}: {eid} claims branch {branch!r} but the "
                            f"evidence does not mention that branch (mismatched branch)")
        # 3. run citations (evidence field AND body -- Round-21 direction 5.3)
        supporting = verify_citations(e, sha) if runs else []
        # 4. live-ref comparison
        tip, err = live(branch)
        if err:
            findings.append(f"FAIL: {path}:{e['n']}: {eid} could not be verified against the "
                            f"live remote ref: {err} (fail closed)")
            continue
        if tip == sha:
            current_for[branch] = sha
            info.append(f"CURRENT: {eid} confirms the live tip of {branch} ({sha[:12]}…)")
            continue
        # 5. HISTORICAL -- Round-21 direction 5: prove the push was real
        anc = logrepo.is_ancestor(sha, tip)
        if anc is True:
            info.append(f"HISTORICAL: {eid} confirmed {branch} @ {sha[:12]}…; the live tip "
                        f"is {tip[:12]}… — SUPPORTED-BY-ANCESTRY (the recorded SHA is on "
                        f"the branch's current lineage; the confirmation is valid for the "
                        f"push it describes and does NOT confirm the newer tip)")
            continue
        # off-lineage (or ancestry not determinable): a reviewed evidence-back
        # record is required (direction 5.5)
        backs = [b for b in entries if b["m"] and b["verb"] == "evidence-back"
                 and field(b["subject"], "backs") == eid
                 and field(b["subject"], "branch") == branch
                 and field(b["subject"], "remote_sha") == sha]
        if not backs:
            findings.append(
                f"FAIL: {path}:{e['n']}: {eid} is classified HISTORICAL ({branch} @ "
                f"{sha[:12]}…, live tip {tip[:12]}…), the SHA is not on the branch's "
                f"current lineage, and no reviewed evidence-back record supports it -- "
                f"a different SHA alone does not prove the historical push was real "
                f"(Round-21 direction 5)")
            continue
        back = backs[-1]
        # the backed SHA must be a real commit: object store or commits API
        obj = logrepo.object_exists(sha)
        api_exists = runs.api_commit_exists(sha) if runs else None
        if not obj and api_exists is not True:
            findings.append(
                f"FAIL: {path}:{e['n']}: the historical SHA {sha[:12]}… backed by "
                f"{back['eid']} is not present in the local object store and could not "
                f"be verified as a real remote commit (fail closed)")
            continue
        # the back's own citations verify (fail-closed policy shared above)
        verify_citations(back, sha)
        # connecting facts (recorded; anything machine-checked strengthens the record)
        connected = False
        _, back_rids = entry_runs(back)
        for rid in back_rids + sorted(set(RUN_ID.findall(
                (field(e["subject"], "evidence") or "") + "\n" + body))):
            rec, _cls = runs.verify(rid) if runs else (None, "no-source-configured")
            if rec and rec.get("head_sha") == sha and rec.get("head_branch") == branch:
                connected = True
                info.append(f"HISTORICAL: {eid} confirmed {branch} @ {sha[:12]}… -- "
                            f"SUPPORTED-BY-EVIDENCE-BACK {back['eid']} "
                            f"(CONNECTED-BY-RUN {rid}, real-commit proof "
                            f"{'object-store' if obj else 'commits-API'})")
                break
        if not connected:
            info.append(f"HISTORICAL: {eid} confirmed {branch} @ {sha[:12]}… -- "
                        f"SUPPORTED-BY-EVIDENCE-BACK {back['eid']} (reviewed record; "
                        f"real-commit proof {'object-store' if obj else 'commits-API'}; "
                        f"no run-connection exists for this head, disclosed)")

    # ---- pass 2: merge-intents
    merge_intents = [e for e in entries if e["m"] and e["verb"] == "merge-intent"]
    voided = {field(e["subject"], "voids") for e in entries
              if e["m"] and e["verb"] == "intent-reconcile"}

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
            for old in targeting[:-1]:
                if old["eid"] in voided:
                    continue
                old_head, new_head = mi_field(old, "expected_head"), mi_field(active, "expected_head")
                if old_head and new_head and logrepo.is_ancestor(old_head, new_head):
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
                if not pr.get("merged"):
                    merge_ok = False
                    findings.append(f"FAIL: {path}:{active['n']}: PR {pr_no} is not "
                                    f"merged (merged={pr.get('merged')!r})")
                if not ((pr.get("head") or {}).get("sha") == head):
                    merge_ok = False
                    findings.append(
                        f"FAIL: {path}:{active['n']}: PR {pr_no} head "
                        f"{str((pr.get('head') or {}).get('sha'))[:12]}… does not "
                        f"match the declared head {(head or '?')[:12]}…")
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
                approved = [r for r in reviews if r.get("state") == "APPROVED"
                            and r.get("user", {}).get("login") != author]
                if not approved:
                    merge_ok = False
                    findings.append(
                        f"FAIL: {path}:{active['n']}: PR {pr_no} has no approving "
                        f"review from a user other than the PR author (required "
                        f"independent approval is not satisfied)")
            if merge_ok:
                info.append(f"MERGE-VERIFIED: {active['eid']} declares the live tip "
                            f"of {branch_under_test} as a reviewed-PR merge "
                            f"(base {base[:12]}… + approved head {head[:12]}…)")
                info.append(f"DECLARED: {branch_under_test} current tip is declared by "
                            f"the verified merge-intent {active['eid']}")

    # ---- pass 3: branch under test -- terminal-intent binding (direction 4)
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
                if declared_by_merge and not unpaired_here:
                    pass  # handled above
                terminal = unpaired_here[-1]
                excess = unpaired_here[:-1]
                for q in excess:
                    if q["e"]["eid"] not in voided:
                        findings.append(
                            f"FAIL: {path}:{q['e']['n']}: multiple unmatched intents for "
                            f"the branch under test: {q['e']['eid']} is not the terminal "
                            f"intent and is not voided by an intent-reconcile record "
                            f"(excess unmatched intents must be separately reconciled)")
                # exact-candidate binding for the terminal intent
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
                                if expect_tip and not logrepo.is_ancestor(tok, expect_tip):
                                    findings.append(
                                        f"FAIL: {path}:{te['n']}: declared pre-push tip "
                                        f"{tok[:12]}… is not an ancestor of the candidate "
                                        f"(the push was not a fast-forward from the "
                                        f"declared state)")
                                elif newest_sha != tok:
                                    findings.append(
                                        f"FAIL: {path}:{te['n']}: declared pre-push tip "
                                        f"{tok[:12]}… does not match the newest recorded "
                                        f"branch state {(newest_sha or 'ABSENT')[:12]}… "
                                        f"({newest_prior['eid'] if newest_prior else 'no prior record'}) — "
                                        f"an intervening ref update or a stale declaration; "
                                        f"reconcile before pushing")
                                else:
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("logfile")
    ap.add_argument("--commit-time", help="ISO-8601 commit time carrying this log version")
    ap.add_argument("--max-drift", type=int, default=300, help="seconds of allowed drift (default 300)")
    ap.add_argument("--parent-prefix", metavar="OLDER_FILE",
                    help="require OLDER_FILE to be a byte-prefix of logfile (append-only)")
    ap.add_argument("--pushlog-schema", action="store_true",
                    help="also enforce push-log schema: unique IDs, well-formed entries with "
                         "required fields, and identifier-based INTENT/CONFIRM pairing "
                         "(unpaired INTENTs are reported as PENDING, not failures)")
    ap.add_argument("--remote-verify", action="store_true",
                    help="Round-20/21 layer: verify every CONFIRM/RECONCILE against the live "
                         "remote ref, prove the evidence supports the recorded SHA, require "
                         "historical pushes to be real (ancestry or a reviewed evidence-back), "
                         "bind the terminal intent to the exact candidate, and verify the "
                         "protected-branch merge transition")
    ap.add_argument("--remote", default="origin",
                    help="remote (URL or path) for ls-remote verification (default: origin)")
    ap.add_argument("--branch", default=None,
                    help="the branch under test (required with --expect-tip)")
    ap.add_argument("--expect-tip", default=None, metavar="SHA",
                    help="the expected live tip of --branch (in CI: the exact pushed SHA)")
    ap.add_argument("--runs-manifest", default=None, metavar="FILE",
                    help="JSON {\"runs\": {id: {head_sha, head_branch?, conclusion}}} for "
                         "run-citation verification")
    ap.add_argument("--repo", default=None, metavar="SLUG",
                    help="owner/repo for GitHub API verification (runs outside the manifest "
                         "window, commit existence, PR + review state); requires the token "
                         "env var")
    ap.add_argument("--api-token-env", default="GH_TOKEN", metavar="VAR",
                    help="environment variable carrying the API token (default GH_TOKEN; "
                         "the value is never printed or written anywhere)")
    ap.add_argument("--pr-json", default=None, metavar="FILE",
                    help="test fixture replacing the PR API response for merge verification")
    ap.add_argument("--reviews-json", default=None, metavar="FILE",
                    help="test fixture replacing the PR-reviews API response")
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
            runs = RunVerifier(manifest, a.repo, a.api_token_env)
            logrepo = LogRepo(a.logfile)
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
    for p in pending:
        print(p)
    for i in info:
        print(i)
    print("\n".join(findings) if findings else f"PASS: {a.logfile} "
          f"{'sequence + append-only ' if a.parent_prefix else ''}"
          f"{'schema ' if (a.pushlog_schema or a.remote_verify) else ''}"
          f"{'remote-verified ' if a.remote_verify else ''}"
          f"clean")
    return 0 if not findings else 1


if __name__ == "__main__":
    sys.exit(main())
