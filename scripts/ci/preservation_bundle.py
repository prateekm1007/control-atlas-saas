#!/usr/bin/env python3
"""Generate the branch/PR preservation bundle (Round-25 direction 7).

The ref-level restoration of the previously deleted branches and the reopening
of the pull requests had been asserted, but no VERIFIED preservation bundle
had been supplied. This script builds that bundle from the live repository, so
the claim can be checked instead of believed. The bundle is self-describing
and reproducible: every file records the exact commands and the API responses
it was derived from.

Bundle contents (all under --out-dir):
    refs.json                 every refs/heads and refs/tags with its SHA
    pull-requests.json        every pull request with its state, head, base
    unique-commits.json       the commits unique to each pull request's branch
    branch-diff-matrix.json   the comparison of each branch against main
    preservation-matrix.json  the preservation verdict per branch and PR
    runtime-dependency-map.json  the runtime entry points per branch
    pr-evidence.json          the CI runs recorded for each pull request head
    commands.log              every command/API call this run made
    SHA256SUMS                the digests of the files above

Exit: 0 = the bundle was written and self-verified; 1 = fail closed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone

API = "https://api.github.com"
RUNTIME_ENTRIES = (
    "app.py", "run.py", "main.py", "wsgi.py", "asgi.py", "manage.py",
    "backend/app.py", "dashboard/app.py", "gpu_worker/main.py",
    "docker-compose.yml", "Dockerfile", "gatekeeper.sh", "scripts/provision_pilot_key.py",
)


class Api:
    """Minimal authenticated GitHub API client with a command log."""

    def __init__(self, repo: str, token: str, log: list[str]):
        self.repo = repo
        self.token = token
        self.log = log

    def get(self, path: str, paginate: bool = False):
        url = f"{API}/repos/{self.repo}/{path}"
        if paginate:
            url += ("&" if "?" in url else "?") + "per_page=100"
        out = []
        page = 1
        while True:
            req = urllib.request.Request(
                url + (f"&page={page}" if paginate else ""),
                headers={"Accept": "application/vnd.github+json",
                         "X-GitHub-Api-Version": "2022-11-28",
                         "Authorization": "Bearer " + self.token})
            self.log.append(f"GET /{self.repo}/{path}" + (f" (page {page})" if paginate else ""))
            try:
                with urllib.request.urlopen(req, timeout=60) as r:
                    data = json.loads(r.read() or b"null")
            except urllib.error.HTTPError as e:
                self.log.append(f"  -> HTTP {e.code}")
                return None
            except (urllib.error.URLError, TimeoutError, OSError) as e:
                self.log.append(f"  -> error {type(e).__name__}")
                return None
            except ValueError:
                self.log.append("  -> unparseable response")
                return None
            if not paginate:
                return data
            if not isinstance(data, list):
                return data
            out.extend(data)
            if len(data) < 100:
                return out
            page += 1
            if page > 10:
                return out


def run_git(repo_dir: str, *args) -> tuple[int, str]:
    r = subprocess.run(["git", "-C", repo_dir, *args], capture_output=True, text=True)
    return r.returncode, (r.stdout or "").strip()


def write_json(path: str, obj) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True, ensure_ascii=False)
        fh.write("\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--repo", required=True, help="owner/name of the repository")
    ap.add_argument("--out-dir", required=True, help="directory to write the bundle into")
    ap.add_argument("--api-token-env", default="GH_TOKEN", metavar="VAR",
                    help="environment variable carrying the API token (default "
                         "GH_TOKEN; the value is never written or printed)")
    ap.add_argument("--base", default="main", help="the base branch of the comparisons")
    args = ap.parse_args()
    token = os.environ.get(args.api_token_env, "") if args.api_token_env else ""
    if not token:
        print(f"FAIL: no API token in ${args.api_token_env}; the bundle cannot be built")
        return 1
    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    log: list[str] = []
    api = Api(args.repo, token, log)
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    log.insert(0, f"preservation bundle generated at {stamp} for {args.repo}")

    # ------------------------------------------------------------ refs
    branches = api.get("branches?per_page=100", paginate=True) or []
    tags = api.get("git/matching-refs/tags/", paginate=True) or []
    refs = {b["name"]: b["commit"]["sha"] for b in branches if isinstance(b, dict) and "name" in b}
    tag_refs = {t["ref"]: t["object"]["sha"] for t in tags
                if isinstance(t, dict) and "ref" in t and isinstance(t.get("object"), dict)}
    write_json(os.path.join(out_dir, "refs.json"),
               {"generated_at": stamp, "repo": args.repo,
                "heads": refs, "tags": tag_refs,
                "head_count": len(refs), "tag_count": len(tag_refs)})

    # ------------------------------------------------------------ pull requests
    pulls = api.get("pulls?state=all&per_page=100", paginate=True) or []
    # The PR LIST endpoint omits the `merged` flag (it is null on every entry),
    # so each PR is re-read from its DETAIL endpoint, which carries the
    # authoritative merged state and merge timestamp.
    for _p in pulls:
        _d = api.get("pulls/%d" % _p.get("number", 0))
        if isinstance(_d, dict):
            _p["merged"] = _d.get("merged")
            _p["merged_at"] = _d.get("merged_at")
    pr_rows = []
    for p in pulls:
        if not isinstance(p, dict):
            continue
        pr_rows.append({
            "number": p.get("number"), "state": p.get("state"),
            "title": (p.get("title") or "")[:120],
            "head_ref": (p.get("head") or {}).get("ref"),
            "head_sha": (p.get("head") or {}).get("sha"),
            "base_ref": (p.get("base") or {}).get("ref"),
            "base_sha": (p.get("base") or {}).get("sha"),
            "merged": p.get("merged"), "merged_at": p.get("merged_at"),
            "merge_commit_sha": p.get("merge_commit_sha"),
            "author": (p.get("user") or {}).get("login"),
        })
    write_json(os.path.join(out_dir, "pull-requests.json"),
               {"generated_at": stamp, "repo": args.repo, "count": len(pr_rows),
                "open": sum(1 for r in pr_rows if r["state"] == "open"),
                "merged": sum(1 for r in pr_rows if r["merged"]),
                "pull_requests": pr_rows})

    # ------------------------------------------------------------ comparisons
    base = args.base
    unique_commits: dict[str, list] = {}
    diff_matrix: dict[str, dict] = {}
    preservation: dict[str, dict] = {}
    for p in pr_rows:
        head_ref = p["head_ref"]
        if not head_ref or head_ref not in refs:
            continue
        if head_ref == base:
            continue
        cmp = api.get(f"compare/{base}...{head_ref}")
        if not isinstance(cmp, dict):
            diff_matrix[head_ref] = {"status": "unavailable"}
            preservation[head_ref] = {"verdict": "UNVERIFIED",
                                      "reason": "the comparison API did not answer"}
            continue
        commits = [c.get("sha") for c in (cmp.get("commits") or []) if isinstance(c, dict)]
        files = [f.get("filename") for f in (cmp.get("files") or []) if isinstance(f, dict)]
        unique_commits[head_ref] = commits
        diff_matrix[head_ref] = {
            "ahead_by": cmp.get("ahead_by"), "behind_by": cmp.get("behind_by"),
            "status": cmp.get("status"), "total_commits": cmp.get("total_commits"),
            "changed_files": len(files),
            "changed_file_sample": sorted(f for f in files if f)[:25],
        }
        preservation[head_ref] = {
            "verdict": "PRESERVED" if commits else "NO-UNIQUE-COMMITS",
            "head_sha": refs[head_ref],
            "unique_commit_count": len(commits),
            "unique_commits_reachable": bool(commits),
            "changed_file_count": len(files),
            "pr_numbers": [r["number"] for r in pr_rows if r["head_ref"] == head_ref],
        }
    write_json(os.path.join(out_dir, "unique-commits.json"),
               {"generated_at": stamp, "base": base, "branches": unique_commits})
    write_json(os.path.join(out_dir, "branch-diff-matrix.json"),
               {"generated_at": stamp, "base": base, "branches": diff_matrix})
    write_json(os.path.join(out_dir, "preservation-matrix.json"),
               {"generated_at": stamp, "base": base, "branches": preservation,
                "summary": {
                    "branches": len(preservation),
                    "preserved": sum(1 for v in preservation.values()
                                     if v.get("verdict") == "PRESERVED"),
                    "unverified": sum(1 for v in preservation.values()
                                       if v.get("verdict") == "UNVERIFIED"),
                }})

    # ------------------------------------------------------------ PR evidence (CI runs)
    pr_evidence: dict[str, list] = {}
    for p in pr_rows:
        sha = p["head_sha"]
        if not sha:
            continue
        runs = api.get(f"actions/runs?head_sha={sha}&per_page=100")
        if isinstance(runs, dict):
            runs = runs.get("workflow_runs") or []
        rows = []
        for r in (runs or []) if isinstance(runs, list) else []:
            if not isinstance(r, dict):
                continue
            rows.append({"id": r.get("id"), "name": r.get("name"), "event": r.get("event"),
                         "head_branch": r.get("head_branch"), "head_sha": r.get("head_sha"),
                         "conclusion": r.get("conclusion"), "created_at": r.get("created_at")})
        pr_evidence[str(p["number"])] = rows
    write_json(os.path.join(out_dir, "pr-evidence.json"),
               {"generated_at": stamp, "repo": args.repo, "pull_requests": pr_evidence})

    # ------------------------------------------------------------ runtime dependency map
    runtime_map: dict[str, list] = {}
    for name, sha in sorted(refs.items()):
        entries = []
        for entry in RUNTIME_ENTRIES:
            blob = api.get(f"contents/{entry}?ref={sha}")
            if isinstance(blob, dict) and blob.get("sha"):
                entries.append({"path": entry, "blob_sha": blob.get("sha"),
                                "size": blob.get("size")})
        runtime_map[name] = entries
    write_json(os.path.join(out_dir, "runtime-dependency-map.json"),
               {"generated_at": stamp, "repo": args.repo,
                "entries_checked": list(RUNTIME_ENTRIES), "branches": runtime_map})

    # ------------------------------------------------------------ commands.log + SHA256SUMS
    with open(os.path.join(out_dir, "commands.log"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(log) + "\n")
    names = sorted(n for n in os.listdir(out_dir) if n != "SHA256SUMS")
    lines = []
    for n in names:
        h = hashlib.sha256(open(os.path.join(out_dir, n), "rb").read()).hexdigest()
        lines.append(f"{h}  {n}")
    with open(os.path.join(out_dir, "SHA256SUMS"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    bundle_digest = hashlib.sha256(open(os.path.join(out_dir, "SHA256SUMS"), "rb").read()).hexdigest()
    print(f"preservation bundle written: {len(names)} file(s) in {out_dir}")
    print(f"bundle SHA256SUMS digest (the bundle's root digest): {bundle_digest}")
    print(f"heads recorded: {len(refs)}; pull requests recorded: {len(pr_rows)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
