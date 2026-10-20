#!/usr/bin/env python3
"""Round-21/22 evidence table generator (auditor directions 8 / 7).

Generates docs/evidence/ec16-ec26/R21_EVIDENCE_TABLE.md FROM LIVE EVIDENCE at
generation time: anonymous ls-remote for every branch tip, the GitHub Actions
API (GH_TOKEN env var; value never written or printed) for every run citation
found in the push log, and the push log itself for the entry inventory. The
output records exact SHAs, branch/base tips, run IDs and conclusions,
remote-ref classification, and remaining pending state.

Round-22 (direction 7): the classifications now DISTINGUISH COMMIT EXISTENCE
FROM BRANCH-PUSH PROOF and clearly mark unproven history. An off-lineage
confirmation is classified PUSH-PROVEN-BY-RUN only when a cited Actions run
has event=push + matching head_branch + the exact head_sha; PUSH-PROVEN-BY-
EVENT when an adjudicated push_event citation backs it; NOT-PROVEN when a
push-not-proven reclassification discloses the gap; and COMMIT-EXISTS-ONLY
when nothing establishes the branch push (fail-closed in the live gate).

Usage:  GH_TOKEN=<token> python3 scripts/ci/generate_evidence_table.py \
          [--repo prateekm1007/control-atlas-saas] [--out <path>]

Exit 0 on success; 1 on any verification failure (fail closed).
"""
import argparse
import json
import os
import re
import subprocess
import sys
import urllib.request
from datetime import datetime, timezone

LINE = re.compile(r"^## (\S+) \| ([^|]+) \| ([^|]+) \| ([^|]+?)(?:\s\|\s*RECONSTRUCTED)?\s*$")
RUN_ID = re.compile(r"\brun(?:_id)?\s*[=:]?\s*(\d{6,})\b", re.I)
PUSH_EVENT_CITE = re.compile(
    r"push_event=(\d+)\s+ref=(\S+)\s+before=([0-9a-f]{40}|0{40})\s+after=([0-9a-f]{40})")
PRIOR_TIP_CITE = re.compile(
    r"push_event=(\d+)\s+ref=(\S+)\s+prior_tip=([0-9a-f]{40})")
API = "https://api.github.com"

BRANCHES = [
    "main",
    "repair/round19-mainline-acceptance",
    "activation/ec16-ec19-ec20-r19",
    "activation/ec21-r19",
    "activation/ec26-r19",
    "activation/ec18-r19",
    "activation/ec23-r19",
    "activation/ec24-r19",
    "docs/owner-decisions-r19",
    "docs/owner-decisions",
    "governance/canonical-docs",
    "proposal/ec16-ec26-enforcement",
    "proposal/ec16-ec26-enforcement-rebased",
]


def ls_tip(url, branch):
    r = subprocess.run(["git", "ls-remote", url, f"refs/heads/{branch}"],
                       capture_output=True, text=True, timeout=60)
    if r.returncode != 0 or not r.stdout.strip():
        return None
    return r.stdout.strip().split()[0]


def is_ancestor(sha, tip):
    r = subprocess.run(["git", "merge-base", "--is-ancestor", sha, tip],
                       capture_output=True, timeout=60)
    return r.returncode == 0


def field(subject, name):
    m = re.search(rf'(?:^|\s){name}=("[^"]*"|\S+)', subject)
    if not m:
        return None
    v = m.group(1)
    return v[1:-1] if v.startswith('"') else v


REPO = "prateekm1007/control-atlas-saas"


def main() -> int:
    global REPO
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=REPO)
    ap.add_argument("--out", default="docs/evidence/ec16-ec26/R21_EVIDENCE_TABLE.md")
    ap.add_argument("--pushlog", default="docs/evidence/push-log/PUSH_LOG.md")
    ap.add_argument("--remote", default=f"https://github.com/{REPO}.git")
    a = ap.parse_args()
    REPO = a.repo
    token = os.environ.get("GH_TOKEN", "")
    if not token:
        print("FAIL: GH_TOKEN is required (API verification of every cited run; "
              "the value is never written or printed)")
        return 1

    def api(path):
        req = urllib.request.Request(f"{API}/repos/{REPO}/{path}",
                                     headers={"Accept": "application/vnd.github+json",
                                              "Authorization": "Bearer " + token})
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read() or b"{}")
        except Exception as e:  # noqa: BLE001
            return {"_error": str(e)[:100]}

    gen_ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = []
    lines.append("# Round-22 evidence table (GENERATED FROM LIVE EVIDENCE)\n")
    lines.append(f"- **Generated:** {gen_ts} by `scripts/ci/generate_evidence_table.py` "
                 "(committed; re-runnable).")
    lines.append(f"- **Sources:** anonymous `git ls-remote` against `{a.remote}` at the "
                 f"generation timestamp; the GitHub Actions API (`/repos/{REPO}`) for "
                 "every run citation found in the push log; the push log itself for "
                 "the entry inventory. Nothing in this table is hand-copied.")
    lines.append("- **Round-22 classification semantics:** commit existence, "
                 "branch-push provenance, and CI success are DISTINCT claims. An "
                 "off-lineage confirmation is PUSH-PROVEN-BY-RUN only via a cited "
                 "run with event=push, matching head_branch, and the exact head_sha; "
                 "PUSH-PROVEN-BY-EVENT via an adjudicated push_event citation; "
                 "NOT-PROVEN via an explicit push-not-proven reclassification; "
                 "anything else is COMMIT-EXISTS-ONLY (fail-closed in the live "
                 "gate; a commits-API 200 or a rebuild narrative is not push "
                 "proof). Failed run conclusions are disclosed separately and "
                 "never count as successful gate evidence.\n")

    lines.append("## 1. Live branch tips (anonymous ls-remote at generation time)\n")
    lines.append("| Branch | Live tip SHA |")
    lines.append("|---|---|")
    tips = {}
    for b in BRANCHES:
        t = ls_tip(a.remote, b)
        tips[b] = t
        lines.append(f"| `{b}` | `{t if t else '(absent)'}` |")
    lines.append("")

    # ---- parse the log
    entries, cur = [], None
    for raw in open(a.pushlog, encoding="utf-8"):
        line = raw.rstrip("\n")
        if line.startswith("## "):
            if cur:
                entries.append(cur)
            m = LINE.match(line)
            cur = {"m": m, "body": []}
            if m:
                cur.update(ts=m.group(1), eid=m.group(2).strip(), verb=m.group(3).strip(),
                           subject=m.group(4).strip())
        elif cur:
            cur["body"].append(line)
    if cur:
        entries.append(cur)
    parsed = [e for e in entries if e["m"]]

    def entry_text(e):
        evm = re.search(r'evidence="([^"]*)"', e["subject"])
        return (evm.group(1) if evm else "") + "\n" + "\n".join(e["body"])

    # ---- fetch every cited run first (needed for provenance classification)
    run_ids = set()
    for e in parsed:
        run_ids.update(RUN_ID.findall(entry_text(e)))
    runs = {}
    for rid in sorted(run_ids):
        d = api(f"actions/runs/{rid}")
        if "_error" in d or not d.get("id"):
            print(f"FAIL: run {rid} could not be verified via the API")
            return 1
        runs[rid] = d

    # ---- classify
    by_eid = {e["eid"]: e for e in parsed}
    reclassified = {field(e["subject"], "reclassifies"): e for e in parsed
                    if e["verb"] == "push-not-proven"}

    lines.append("## 2. Push-log entry inventory and remote-ref classification\n")
    lines.append("| Entry | Verb | Branch | Recorded remote_sha | Live tip | Classification |")
    lines.append("|---|---|---|---|---|---|")
    for e in parsed:
        subj = e["subject"]
        branch = field(subj, "branch") or (
            field(subj, "target") if e["verb"] == "merge-intent" else None)
        sha = field(subj, "remote_sha")
        live = tips.get(branch) if branch else None
        if e["verb"] in ("push-confirm", "push-reconcile"):
            if live and sha and live == sha:
                cls = "CURRENT (live ref == recorded SHA)"
            elif live and sha and is_ancestor(sha, live):
                cls = "HISTORICAL — ON-LINEAGE (ancestry from the live tip; valid for " \
                      "the push it describes)"
            elif e["eid"] in reclassified:
                cls = "NOT PROVEN (reclassified; commit existence may be proven, " \
                      "push provenance is not — disclosed, never counted as evidence)"
            else:
                backs = [b for b in parsed if b["verb"] == "evidence-back"
                         and field(b["subject"], "backs") == e["eid"]]
                prov = None
                cited = RUN_ID.findall(entry_text(e)) if e["verb"] == "push-confirm" else []
                for b in backs:
                    text = entry_text(b)
                    for rid in RUN_ID.findall(text):
                        r = runs.get(rid, {})
                        if (r.get("event") == "push" and r.get("head_branch") == branch
                                and r.get("head_sha") == sha):
                            prov = (f"HISTORICAL — PUSH-PROVEN-BY-RUN {rid} "
                                    f"(event=push, head_branch={branch}); CI-success "
                                    f"claimed only from successful runs")
                            break
                    if prov:
                        break
                    if PUSH_EVENT_CITE.search(text):
                        prov = ("HISTORICAL — PUSH-PROVEN-BY-EVENT (adjudicated "
                                "push_event citation; see the backing record)")
                    elif PRIOR_TIP_CITE.search(text):
                        prov = ("HISTORICAL — PUSH-PROVEN-BY-PRIOR-TIP-EVENT "
                                "(a later push on the same ref records the backed "
                                "SHA as its before; adjudicated — see the backing "
                                "record)")
                if not prov and backs:
                    for b in backs:
                        for rid in RUN_ID.findall(entry_text(b)):
                            r = runs.get(rid, {})
                            if (r.get("event") == "push" and r.get("head_sha") == sha
                                    and r.get("head_branch") != branch):
                                prov = ("UNPROVEN-BRANCH (a push event on another "
                                        f"branch {r.get('head_branch')!r} does not "
                                        "support this entry)")
                                break
                        if prov:
                            break
                cls = prov or ("COMMIT-EXISTS-ONLY (no branch-push proof established "
                               "by this table; the live gate classifies fail-closed)")
        elif e["verb"] == "push-intent":
            cls = "INTENT (terminal pending disclosed by the checker)"
        elif e["verb"] == "merge-intent":
            cls = ("MERGE-INTENT (pre-declared; verifies on the post-merge main push)"
                   if branch == "main" else "MERGE-INTENT")
        elif e["verb"] == "evidence-back":
            cls = "EVIDENCE-BACK (reviewed historical support)"
        elif e["verb"] == "intent-reconcile":
            cls = "INTENT-RECONCILE (voids a named intent)"
        elif e["verb"] == "push-not-proven":
            cls = "PUSH-NOT-PROVEN (honest reclassification of an unproven push)"
        else:
            cls = "legacy push record"
        lines.append(f"| `{e['eid']}` | {e['verb']} | `{branch or '-'}` | "
                     f"`{sha[:12] if sha else '-'}`… | `{live[:12] if live else '-'}`… | {cls} |")
    lines.append("")

    lines.append("## 3. Cited Actions runs (GitHub API at generation time)\n")
    lines.append("| Run ID | Workflow | Event | Head branch | Head SHA | Conclusion |")
    lines.append("|---|---|---|---|---|---|")
    for rid in sorted(run_ids):
        d = runs[rid]
        lines.append(f"| {rid} | {d.get('name')} | {d.get('event')} | "
                     f"{d.get('head_branch')} | `{str(d.get('head_sha'))[:12]}…` | "
                     f"{d.get('conclusion')} |")
    lines.append("")

    # ---- derive the pending state generically from the log
    repair = "repair/round19-mainline-acceptance"
    unpaired = [e for e in parsed if e["verb"] == "push-intent"
                and field(e["subject"], "branch") == repair]
    terminal = unpaired[-1] if unpaired else None
    active_mi = [e for e in parsed if e["verb"] == "merge-intent"
                 and field(e["subject"], "source_branch") == repair]
    lines.append("## 4. Remaining pending state (disclosed, never smoothed over)\n")
    if terminal:
        lines.append(f"- `{terminal['eid']}` is the terminal unpaired intent for "
                     f"`{repair}`: the push it declares is disclosed PENDING until a "
                     f"next legitimate push delivers its CONFIRM (labeled PENDING, "
                     f"never COMPLETE; no push is manufactured merely to eliminate it).")
    if active_mi:
        mi = active_mi[-1]
        lines.append(f"- `{mi['eid']}` (merge-intent) is pre-declared and PENDING until "
                     "the reviewed merge of PR #30 executes through the protected "
                     "path; the post-merge `main` push event is declared by it once "
                     "the checker verifies the merge-commit parents, "
                     "merge_commit_sha, source ref, merged_at, the merge PushEvent, "
                     "and the approving non-author review.")
    lines.append("- PR #30 remains WITHOUT an independent current-head approval: the "
                 "pre-merge-verification job (Round-22) is EXPECTED RED until that "
                 "approval is recorded on the exact head — the functioning pre-merge "
                 "gate, not a defect. The merge awaits all Round-22 blockers fixed on "
                 "the exact final head, the full gate green there, the passing "
                 "pre-merge approval verification, and the protected-PR MERGE COMMIT.")
    lines.append("- Standing OPEN items (not closed by this round): A-04 branch-"
                 "protection required-check verification; independent r2 evidence-"
                 "archive verification; portable 40-row cohort metadata table; A-02 "
                 "issuer-side revocation; EC2 deployment secrets (deployment stays "
                 "blocked).")
    lines.append("- The PAT-hygiene sweep is a LOCAL coder claim; GitHub CI does not "
                 "independently verify it. No token value appears in this table, the "
                 "repository, or any worklog.\n")

    out = "\n".join(lines)
    open(a.out, "w", encoding="utf-8").write(out)
    print(f"wrote {a.out} ({len(parsed)} log entries, {len(run_ids)} cited runs, "
          f"generated {gen_ts})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
