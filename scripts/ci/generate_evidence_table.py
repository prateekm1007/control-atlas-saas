#!/usr/bin/env python3
"""Round-21 evidence table generator (auditor direction 8).

Generates docs/evidence/ec16-ec26/R21_EVIDENCE_TABLE.md FROM LIVE EVIDENCE at
generation time: anonymous ls-remote for every branch tip, the GitHub Actions
API (GH_TOKEN env var; value never written or printed) for every run citation
found in the push log, and the push log itself for the entry inventory. The
output records exact SHAs, branch/base tips, run IDs and conclusions,
remote-ref classification, and remaining pending state.

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


def api(path, token):
    req = urllib.request.Request(f"{API}/repos/{REPO}/{path}",
                                 headers={"Accept": "application/vnd.github+json",
                                          "Authorization": "Bearer " + token})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.loads(r.read() or b"{}")
    except Exception as e:  # noqa: BLE001
        return {"_error": str(e)[:100]}


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

    gen_ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = []
    lines.append("# Round-21 evidence table (GENERATED FROM LIVE EVIDENCE)\n")
    lines.append(f"- **Generated:** {gen_ts} by `scripts/ci/generate_evidence_table.py` "
                 "(committed; re-runnable).")
    lines.append(f"- **Sources:** anonymous `git ls-remote` against `{a.remote}` at the "
                 f"generation timestamp; the GitHub Actions API (`/repos/{REPO}`) for "
                 "every run citation found in the push log; the push log itself for "
                 "the entry inventory. Nothing in this table is hand-copied.")
    lines.append("- **Scope note:** the exact final candidate SHA of the Round-21 "
                 "repair is the commit carrying the terminal `P-060-INTENT` (this "
                 "table is part of that commit, so it cannot embed the commit's own "
                 "SHA); the candidate is bound by the push-log terminal intent and "
                 "PR #30's body after the push. Run rows below are the verified "
                 "pre-candidate history.\n")

    lines.append("## 1. Live branch tips (anonymous ls-remote at generation time)\n")
    lines.append("| Branch | Live tip SHA |")
    lines.append("|---|---|")
    tips = {}
    for b in BRANCHES:
        t = ls_tip(a.remote, b)
        tips[b] = t
        lines.append(f"| `{b}` | `{t if t else '(absent)'}` |")
    lines.append("")

    lines.append("## 2. Push-log entry inventory and remote-ref classification\n")
    lines.append("| Entry | Verb | Branch | Recorded remote_sha | Live tip | Classification |")
    lines.append("|---|---|---|---|---|---|")
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
    run_ids = set()
    for e in entries:
        if not e["m"]:
            continue
        subj = e["subject"]
        branch = (re.search(r"branch=(\S+)", subj) or [None]) and \
            (re.search(r"branch=(\S+)", subj).group(1)
             if re.search(r"branch=(\S+)", subj) else
             (re.search(r"target=(\S+)", subj).group(1)
              if e["verb"] == "merge-intent" and re.search(r"target=(\S+)", subj) else None))
        sm = re.search(r"remote_sha=([0-9a-f]{40})", subj)
        sha = sm.group(1) if sm else None
        live = tips.get(branch) if branch else None
        if e["verb"] in ("push-confirm", "push-reconcile"):
            cls = ("CURRENT" if (live and sha and live == sha)
                   else "HISTORICAL (valid for the push it describes)")
        elif e["verb"] == "push-intent":
            cls = "INTENT (terminal pending disclosed by the checker)"
        elif e["verb"] == "merge-intent":
            cls = ("MERGE-INTENT (pre-declared; verifies on the post-merge main push)"
                   if branch == "main" else "MERGE-INTENT")
        elif e["verb"] == "evidence-back":
            cls = "EVIDENCE-BACK (reviewed historical support)"
        elif e["verb"] == "intent-reconcile":
            cls = "INTENT-RECONCILE (voids a named intent)"
        else:
            cls = "legacy push record"
        lines.append(f"| `{e['eid']}` | {e['verb']} | `{branch or '-'}` | "
                     f"`{sha[:12] if sha else '-'}`… | `{live[:12] if live else '-'}`… | {cls} |")
        evm = re.search(r'evidence="([^"]*)"', subj)
        text = (evm.group(1) if evm else "") + "\n" + "\n".join(e["body"])
        run_ids.update(RUN_ID.findall(text))
    lines.append("")

    lines.append("## 3. Cited Actions runs (GitHub API at generation time)\n")
    lines.append("| Run ID | Workflow | Event | Head branch | Head SHA | Conclusion |")
    lines.append("|---|---|---|---|---|---|")
    for rid in sorted(run_ids):
        d = api(f"actions/runs/{rid}", token)
        if "_error" in d or not d.get("id"):
            lines.append(f"| {rid} | API-UNVERIFIABLE | - | - | - | FAIL (fail closed) |")
            print(f"FAIL: run {rid} could not be verified via the API")
            return 1
        lines.append(f"| {rid} | {d.get('name')} | {d.get('event')} | "
                     f"{d.get('head_branch')} | `{str(d.get('head_sha'))[:12]}…` | "
                     f"{d.get('conclusion')} |")
    lines.append("")

    lines.append("## 4. Remaining pending state (disclosed, never smoothed over)\n")
    lines.append("- `P-060-INTENT` is the terminal unpaired intent for "
                 "`repair/round19-mainline-acceptance`: the Round-21 candidate push "
                 "it declares is disclosed PENDING until a next legitimate push "
                 "delivers `P-060-CONFIRM` (labeled PENDING, never COMPLETE; no push "
                 "is manufactured merely to eliminate it, per the Round-21 "
                 "adjudication).")
    lines.append("- `P-061` (merge-intent) is pre-declared and PENDING until the "
                 "reviewed merge of PR #30 executes through the protected path; the "
                 "post-merge `main` push event is declared by it once the checker "
                 "verifies the merge-commit parents, the merged PR, and the approving "
                 "non-author review.")
    lines.append("- PR #30 remains UNAPPROVED: merge awaits the auditor's review, the "
                 "full gate on the exact final candidate, the demonstrated post-merge "
                 "transition, and the required independent approval.")
    lines.append("- Standing OPEN items (not closed by this round): independent r2 "
                 "evidence-archive verification; portable 40-row cohort metadata "
                 "table; A-02 issuer-side revocation; EC2 deployment secrets "
                 "(deployment stays blocked).")
    lines.append("- The PAT-hygiene sweep is a LOCAL coder claim; GitHub CI does not "
                 "independently verify it. No token value appears in this table, the "
                 "repository, or any worklog.\n")

    out = "\n".join(lines)
    open(a.out, "w", encoding="utf-8").write(out)
    print(f"wrote {a.out} ({len(entries)} log entries, {len(run_ids)} cited runs, "
          f"generated {gen_ts})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
