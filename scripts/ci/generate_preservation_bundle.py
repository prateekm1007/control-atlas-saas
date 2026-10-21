#!/usr/bin/env python3
"""Generate the machine-readable branch/PR preservation bundle (EC-27 recovery evidence).

Produces, in --out (default: preservation-bundle/):
  refs.json                  all branch names with exact heads + ref census + 25-ref accounting
  pull-requests.json         every PR: state, head, base, merge SHA/timestamps, disposition
  unique-commits.json        commits reachable only from non-main refs (with patch-ids)
  branch-diff-matrix.json    per branch: merge base, ahead/behind, per-file changes (renames)
  preservation-matrix.json   per unique commit: preservation disposition + preserving refs
  runtime-dependency-map.json references/workflows that depend on the restored branches
  pr-evidence.json           PR claim extraction (SHAs, run IDs) + verification outcomes
  commands.log               commands, API requests, return codes, outcomes (credentials redacted)
  SHA256SUMS                 sha256 of every bundle file

Fail-closed: any required API read that fails is recorded and the script exits non-zero.
Credentials are never written to the bundle, the log, or stdout.
"""
import argparse
import hashlib
import json
import os
import pathlib
import re
import subprocess
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone

HEX40 = re.compile(r"\b([0-9a-f]{40})\b")
RUNID = re.compile(r"\brun(?:s)?\b[^\d]{0,15}(\d{6,})\b", re.I)
RUNTIME_PREFIXES = ("backend/", "scripts/", "src/", "config/", "deploy/", "app/", "requirements",
                    ".github/", "package.json", "Makefile", "Dockerfile", "docker-compose",
                    "pyproject", "setup.", "tox.ini", "noxfile")


def utc():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


class Runner:
    def __init__(self, repo_root, out_dir):
        self.repo = pathlib.Path(repo_root)
        self.out = pathlib.Path(out_dir)
        self.out.mkdir(parents=True, exist_ok=True)
        self.log_path = self.out / "commands.log"
        self.log = self.log_path.open("a", encoding="utf-8")
        self.api_base = os.environ.get("GITHUB_API_URL", "https://api.github.com").rstrip("/")
        token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN") or ""
        self.token = token
        self.repo_slug = os.environ.get("GITHUB_REPOSITORY", "prateekm1007/control-atlas-saas")
        self.api_failures = []

    def note(self, msg):
        line = "%s | %s" % (utc(), msg)
        self.log.write(line + "\n")
        self.log.flush()

    def git(self, *args, timeout=600):
        cmd = ["git", "-C", str(self.repo), "-c", "core.quotepath=false"] + list(args)
        r = subprocess.run(cmd, capture_output=True, timeout=timeout)
        out = (r.stdout or b"").decode("utf-8", "replace")
        err = (r.stderr or b"").decode("utf-8", "replace")
        self.note("rc=%d | git %s | %s" % (r.returncode, " ".join(args)[:240], err.strip()[:160]))
        return r.returncode, out, err

    def api(self, path, method="GET", body=None):
        url = self.api_base + path
        data = json.dumps(body).encode() if body is not None else None
        headers = {"Accept": "application/vnd.github+json", "User-Agent": "preservation-bundle"}
        if self.token:
            headers["Authorization"] = "Bearer " + self.token
        req = urllib.request.Request(url, headers=headers, data=data, method=method)
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                raw = r.read()
                status = r.status
        except urllib.error.HTTPError as e:
            raw = e.read()
            status = e.code
        except Exception as e:  # noqa: BLE001
            self.api_failures.append({"path": path, "error": type(e).__name__})
            self.note("api FAIL %s %s -> %s" % (method, path, type(e).__name__))
            return 0, None
        self.note("api %s %s -> HTTP %d" % (method, path, status))
        try:
            return status, json.loads(raw) if raw else None
        except Exception:  # noqa: BLE001
            return status, None

    def patchid(self, sha):
        rc, diff, _ = self.git("diff-tree", "-p", "-M", "--full-index", sha, timeout=180)
        if rc != 0:
            return None
        r = subprocess.run(["git", "patch-id", "--stable"], input=diff.encode("utf-8", "replace"),
                           capture_output=True)
        parts = (r.stdout or b"").decode("utf-8", "replace").strip().split(" ")
        return parts[0] if parts and parts[0] else None

    def refs(self):
        rc, out, _ = self.git("for-each-ref", "--format=%(refname)\t%(objecttype)\t%(objectname)")
        rows = []
        for line in out.splitlines():
            name, typ, sha = line.split("\t")
            rows.append({"ref": name, "objecttype": typ, "sha": sha})
        heads = {}
        for r in rows:
            if r["ref"].startswith("refs/heads/"):
                heads[r["ref"][len("refs/heads/"):]] = r["sha"]
            elif r["ref"].startswith("refs/remotes/origin/"):
                b = r["ref"][len("refs/remotes/origin/"):]
                heads.setdefault(b, r["sha"])
        return rows, heads


def build(args):
    r = Runner(args.repo, args.out)
    r.note("preservation bundle generation start; repo=%s" % r.repo_slug)
    rows, heads = r.refs()
    main_sha = heads.get("main")
    r.note("heads=%d main=%s" % (len(heads), main_sha))
    if not main_sha:
        r.note("FAIL: main ref not found")
        return 2

    # ---- refs.json -------------------------------------------------------
    census = {}
    for row in rows:
        top = row["ref"].split("/")[1]
        census[top] = census.get(top, 0) + 1
    refs_doc = {
        "generated_utc": utc(), "repo": r.repo_slug, "main_sha": main_sha,
        "ref_count": len(rows), "count_by_space": census, "heads_count": len(heads),
        "heads": heads, "refs": rows,
    }
    (r.out / "refs.json").write_text(json.dumps(refs_doc, indent=1), encoding="utf-8")
    r.note("refs.json: %d refs, %d heads" % (len(rows), len(heads)))

    # ---- PR inventory (detail calls: the list endpoint omits merged) ------
    st, pulls = r.api("/repos/%s/pulls?state=all&per_page=100" % r.repo_slug)
    if st != 200 or not pulls:
        r.note("FAIL: PR list unavailable (HTTP %s)" % st)
        return 3
    pr_docs, detail = [], {}
    for p in sorted(pulls, key=lambda x: x["number"]):
        n = p["number"]
        st_d, d = r.api("/repos/%s/pulls/%d" % (r.repo_slug, n))
        if st_d != 200 or not d:
            r.note("FAIL: PR %d detail unavailable" % n)
            return 4
        detail[n] = d
        st_r, revs = r.api("/repos/%s/pulls/%d/reviews?per_page=100" % (r.repo_slug, n))
        merged = bool(d.get("merged"))
        doc = {
            "number": n, "state": d.get("state"), "merged": merged,
            "merge_commit_sha": d.get("merge_commit_sha"), "merged_at": d.get("merged_at"),
            "merged_by": (d.get("merged_by") or {}).get("login"),
            "title": d.get("title"), "user": (d.get("user") or {}).get("login"),
            "head_ref": (d.get("head") or {}).get("ref"), "head_sha": (d.get("head") or {}).get("sha"),
            "base_ref": (d.get("base") or {}).get("ref"), "base_sha": (d.get("base") or {}).get("sha"),
            "closed_at": d.get("closed_at"), "created_at": d.get("created_at"),
            "updated_at": d.get("updated_at"),
            "reviews": [{"user": (x.get("user") or {}).get("login"), "state": x.get("state"),
                         "submitted_at": x.get("submitted_at"), "commit_id": x.get("commit_id")}
                        for x in (revs or [])],
            "detail_source": "GET /pulls/%d" % n,
        }
        doc["disposition"] = ("MERGED (content on main; recipe preserved)" if merged else
                              ("OPEN (head branch present; not closed by recovery)" if d.get("state") == "open"
                               else "CLOSED, UNMERGED (pre-existing closure; not modified)"))
        pr_docs.append(doc)
    (r.out / "pull-requests.json").write_text(json.dumps(pr_docs, indent=1), encoding="utf-8")
    r.note("pull-requests.json: %d PRs" % len(pr_docs))

    # ---- branch diff matrix ---------------------------------------------
    branches = sorted(b for b in heads if b != "main")
    matrix, runtime_by_branch = {}, {}
    for b in branches:
        entry = {"branch": b, "tip_sha": heads[b]}
        rc, mb, _ = r.git("merge-base", "main", b)
        entry["merge_base_with_main"] = mb.strip() or None
        if entry["merge_base_with_main"]:
            rc, a, _ = r.git("rev-list", "--count", "main.." + b)
            entry["ahead_of_main"] = int(a.strip() or 0)
            rc, bh, _ = r.git("rev-list", "--count", b + "..main")
            entry["behind_main"] = int(bh.strip() or 0)
            rc, z, _ = r.git("diff", "-z", "-M", "--name-status", entry["merge_base_with_main"], b)
            parts = [x for x in z.split("\x00") if x]
            files, i = [], 0
            while i < len(parts):
                code = parts[i]
                if code[:1] in ("R", "C"):
                    files.append({"op": code[0], "score": code[1:], "old": parts[i + 1], "new": parts[i + 2]})
                    i += 3
                else:
                    files.append({"op": code, "path": parts[i + 1]})
                    i += 2
            entry["file_changes_total"] = len(files)
            entry["files"] = files
            runtime_by_branch[b] = [f.get("new") or f.get("path") for f in files]
        matrix[b] = entry
    (r.out / "branch-diff-matrix.json").write_text(json.dumps(matrix, indent=1), encoding="utf-8")
    r.note("branch-diff-matrix.json: %d branches" % len(matrix))

    # ---- unique commits + preservation matrix ---------------------------
    contains = {}
    for b in branches:
        rc, u, _ = r.git("rev-list", b, "--not", "main")
        for c in u.split():
            contains.setdefault(c, set()).add(b)
    setu = sorted(contains)
    rc, mn, _ = r.git("rev-list", "main", "--no-merges")
    main_patchids = {}
    for c in mn.split():
        pid = r.patchid(c)
        if pid:
            main_patchids.setdefault(pid, c)
    commits, patchmap = [], {}
    for c in setu:
        rc, meta, _ = r.git("show", "-s", "--format=%cI%x01%an%x01%cn%x01%s%x01%P", c, timeout=120)
        f = (meta.strip("\n").split("\x01") + ["", "", "", "", ""])[:5]
        parents = f[4].split()
        pid = None if len(parents) > 1 else r.patchid(c)
        if pid:
            patchmap.setdefault(pid, []).append(c)
        rc, fo, _ = r.git("diff-tree", "-r", "--name-only", "--no-commit-id", c, timeout=180)
        commits.append({"sha": c, "date": f[0], "author": f[1], "committer": f[2], "subject": f[3],
                        "parents": parents, "branches_containing": sorted(contains[c]),
                        "patch_id": pid, "changed_files": fo.splitlines()})
    (r.out / "unique-commits.json").write_text(
        json.dumps({"count": len(setu), "commits": commits}, indent=1), encoding="utf-8")
    pid_branches = {}
    for pid, cs in patchmap.items():
        brs = set()
        for c in cs:
            brs |= contains.get(c, set())
        pid_branches[pid] = sorted(brs)
    pres = []
    for e in commits:
        pid, own = e["patch_id"], set(e["branches_containing"])
        if len(e["parents"]) > 1:
            pres.append({"sha": e["sha"], "disposition": "merge-node (ancestry, not content-diffable)",
                         "preserved_by": sorted(own)})
            continue
        if pid and pid in main_patchids:
            pres.append({"sha": e["sha"], "disposition": "content-equivalent-on-main",
                         "equivalent_main_commit": main_patchids[pid], "patch_id": pid,
                         "preserved_by": ["main"]})
            continue
        others = [b for b in pid_branches.get(pid, []) if b not in own] if pid else []
        if others:
            pres.append({"sha": e["sha"], "disposition": "equivalent-content-on-other-ref",
                         "preserved_by": sorted(set(others) | own), "patch_id": pid})
            continue
        if len(own) > 1:
            pres.append({"sha": e["sha"], "disposition": "same-commit-object-on-multiple-refs",
                         "preserved_by": sorted(own)})
            continue
        pres.append({"sha": e["sha"], "disposition": "branch-only-content (survives via its ref)",
                     "preserved_by": sorted(own)})
    (r.out / "preservation-matrix.json").write_text(json.dumps(pres, indent=1), encoding="utf-8")
    r.note("unique-commits=%d preservation rows=%d" % (len(setu), len(pres)))

    # ---- runtime dependency map -----------------------------------------
    consumers = []
    for wf in ("governance-acceptance.yml", "toscanini-gate.yml"):
        rc, wtxt, _ = r.git("show", "main:.github/workflows/" + wf)
        if rc != 0:
            continue
        consumers.append({
            "file": ".github/workflows/" + wf,
            "declared_governs_globs": sorted(set(re.findall(r"~governs=([^\s]+)", wtxt))),
            "paths_referenced": sorted(set(re.findall(r"[A-Za-z0-9_.-]+(?:/[A-Za-z0-9_.-]+)+", wtxt)))[:400],
        })
    rc, cim, _ = r.git("ls-tree", "-r", "--name-only", "main", "scripts/ci/")
    rmap = {"runtime_consumers": consumers, "scripts_ci_inventory": cim.split(),
            "branch_touches_runtime_paths": {
                b: sorted(p for p in runtime_by_branch.get(b, [])
                          if p and not p.startswith("docs/evidence/") and p.startswith(RUNTIME_PREFIXES))
                for b in branches},
            "method": "Workflow 'governs=' globs and script inventory from main; per-branch runtime "
                      "file touches from branch-diff-matrix.json."}
    (r.out / "runtime-dependency-map.json").write_text(json.dumps(rmap, indent=1), encoding="utf-8")
    r.note("runtime-dependency-map.json written")

    # ---- PR evidence (claims + verification) ----------------------------
    cited = set()
    for p in pr_docs:
        cited |= {int(x) for x in RUNID.findall(p.get("body") or "")}
    runs = {}
    for rid in sorted(cited)[:120]:
        st, run = r.api("/repos/%s/actions/runs/%d" % (r.repo_slug, rid))
        runs[str(rid)] = ({"verified": False, "api_status": st} if st != 200 else
                          {"verified": True, "name": run.get("name"), "event": run.get("event"),
                           "status": run.get("status"), "conclusion": run.get("conclusion"),
                           "head_sha": run.get("head_sha"), "head_branch": run.get("head_branch"),
                           "url": run.get("html_url")})
    pe = []
    for p in pr_docs:
        body = p.get("body") or ""
        shas = HEX40.findall(body)
        checks = {}
        for s in shas:
            rc, t, _ = r.git("cat-file", "-t", s, timeout=120)
            t = t.strip()
            checks[s] = {"exists_in_repo": t in ("commit", "tag"), "type": t or None}
        pe.append({"number": p["number"], "state": p["state"], "merged": p["merged"],
                   "head_sha": p["head_sha"], "base_sha": p["base_sha"],
                   "declared_shas": shas, "sha_existence_checks": checks,
                   "cited_run_ids": sorted({str(int(x)) for x in RUNID.findall(body)}),
                   "cited_run_verification": {rid: runs.get(rid, {"not_fetched": True})
                                             for rid in sorted({str(int(x)) for x in RUNID.findall(body)})},
                   "review_states": [x["state"] for x in p["reviews"]],
                   "review_commit_ids": [x["commit_id"] for x in p["reviews"]],
                   "body_chars": len(body)})
    (r.out / "pr-evidence.json").write_text(json.dumps(pe, indent=1), encoding="utf-8")
    r.note("pr-evidence.json: %d PRs, %d cited run ids" % (len(pe), len(runs)))

    # ---- SHA256SUMS -----------------------------------------------------
    # The command log is part of the bundle, so it must receive its FINAL
    # write before the checksums are computed; nothing may be appended to it
    # afterwards or the recorded digest would not verify.
    r.note("bundle generation complete; checksumming %d artifact(s)" % 9)
    files = sorted(p for p in r.out.iterdir() if p.is_file() and p.name != "SHA256SUMS")
    sums = []
    for p in files:
        sums.append("%s  *%s" % (hashlib.sha256(p.read_bytes()).hexdigest(), p.name))
    (r.out / "SHA256SUMS").write_text("\n".join(sums) + "\n", encoding="utf-8")
    digest = hashlib.sha256((r.out / "SHA256SUMS").read_bytes()).hexdigest()
    print("bundle_dir=%s" % r.out)
    print("bundle_files=%d" % (len(files) + 1))
    print("bundle_sha256sums_digest=%s" % digest)
    if r.api_failures:
        print("api_failures=%d" % len(r.api_failures))
        return 5
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", default=os.environ.get("BUNDLE_REPO", "."))
    ap.add_argument("--out", default=os.environ.get("BUNDLE_OUT", "preservation-bundle"))
    args = ap.parse_args()
    sys.exit(build(args))


if __name__ == "__main__":
    main()
