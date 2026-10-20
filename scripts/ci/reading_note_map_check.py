#!/usr/bin/env python3
"""EC-26 reading-note mapping check (inert until adopted).

Auditor revision 2026-10-20 (fail-closed redesign):
  - The required rule set is DERIVED from the machine-readable rule inventory
    table in docs/directives/INDEX.md (the BEGIN/END RULE INVENTORY block).
    It is NEVER taken from the PR body or from an arbitrary --rules list.
  - Every enforcement point named in the inventory must EXIST:
      script:PATH / procedural:PATH  -> PATH must exist (relative to repo root)
      workflow:FILE#JOB              -> JOB must exist in workflow FILE
      pytest:NODE                    -> the file part of NODE must exist
    For ACTIVE rules, a `workflow:` reference must additionally be LIVE: unpinned,
    under `.github/workflows/`, and present in the working tree (Round-19 fix --
    a `proposals/ci/...` reference is inert and must not satisfy an ACTIVE rule).
  - Round-20 trigger-liveness fix (auditor direction 5): liveness is more than
    path + job existence. A workflow under `.github/workflows/` whose `on:`
    block contains only `workflow_dispatch` would satisfy the Round-19 test
    without enforcing any push or pull_request. An ACTIVE rule's `workflow:`
    reference MUST therefore declare the events its enforcement claims, with
    the suffix `~events=push+pull_request`, and the live workflow MUST
    actually be triggered by every declared event with branch filters that
    include the enforcement branch (manual-only, branch-disabled, or otherwise
    untriggered workflows FAIL). PROPOSED rules may cite inert drafts (their
    activation path is the reviewed copy into the live path).
  - Every acceptance test named in the inventory must EXIST:
      AT-nn    -> key in scripts/ci/tests/ACCEPTANCE_TESTS.json
      pytest:NODE / workflow:FILE#JOB -> as above
      NONE     -> declared absence (allowed; reported as a warning, not hidden)
  - The reading note must (a) exist at the canonical path declared inside the
    INDEX block ("reading-note-path:"), (b) cite every pinned blob listed in
    INDEX ("pinned-blobs:"), and (c) map every required rule ID (a line
    containing the ID together with "->" or a markdown table row "|").
  - Missing input, missing notes, empty rule sets, duplicate rule IDs,
    malformed rows, unknown status tokens, or named-but-missing checks FAIL
    closed (exit 1).

Exit codes: 0 = pass; 1 = fail; 2 = usage/internal error.
"""
import argparse
import fnmatch
import json
import os
import subprocess
import re
import sys

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

BEGIN = "<!-- BEGIN RULE INVENTORY"
END = "<!-- END RULE INVENTORY"
STATUS_RE = re.compile(r"^\|\s*([^|]+?)\s*\|\s*([^|]+?)\s*\|\s*([^|]*?)\s*\|\s*([^|]*?)\s*\|$")
KEYVAL_RE = re.compile(r"^([a-z][a-z0-9-]*):\s*(.+)$", re.I)
VALID_STATUS = {"ACTIVE", "PROPOSED"}
ITEM_SPLIT = re.compile(r"[,;]\s*")
ITEM_RE = re.compile(r"^(NONE|AT-\d+|(script|procedural|workflow|pytest):(.+))$", re.S)


def _path_exists(repo: str, path: str) -> bool:
    """PATH may carry a @REF suffix for files living on another branch/ref."""
    if "@" in path:
        fpath, _, ref = path.rpartition("@")
        if fpath and ref and os.path.exists(os.path.join(repo, fpath)):
            return True
        return subprocess.run(["git", "cat-file", "-e", f"{ref}:{fpath}"],
                              capture_output=True, cwd=repo).returncode == 0
    return os.path.exists(os.path.join(repo, path))


def _yaml_doc(repo: str, filespec: str):
    """Load a workflow YAML from @REF when given (explicit pin wins), else the working tree."""
    fpath, _, ref = (filespec.rpartition("@") if "@" in filespec else (filespec, "", ""))
    if ref:
        t = subprocess.run(["git", "show", f"{ref}:{fpath}"], capture_output=True, text=True,
                           cwd=repo).stdout
        if not t:
            return None
        return yaml.safe_load(t) or {}
    fp = os.path.join(repo, fpath)
    if os.path.exists(fp):
        return yaml.safe_load(open(fp, encoding="utf-8")) or {}
    return None


TRIGGER_EVENTS = {"push", "pull_request", "pull_request_target"}


def _on_block(doc: dict):
    """Return the workflow's `on:` value, handling the YAML 1.1 `on`->True gotcha."""
    for key in (True, "on", "On", "ON"):
        if isinstance(doc, dict) and key in doc:
            return doc[key]
    return None


def _validate_triggers(filespec: str, doc: dict, events: list[str],
                       branch: str, findings: list[str]):
    """Round-20: the live workflow must actually run on every declared event.

    Fails closed for: manual-only workflows (workflow_dispatch and friends),
    missing declared events, and branch filters that exclude the enforcement
    branch (branches / branches-ignore).
    """
    on_val = _on_block(doc)
    if on_val is None:
        findings.append(f"workflow {filespec} declares no 'on:' triggers at all "
                        f"(cannot enforce {', '.join(events)})")
        return
    if isinstance(on_val, str):
        on_val = [on_val]
    if isinstance(on_val, list):
        on_val = {str(e): None for e in on_val}
    if not isinstance(on_val, dict):
        findings.append(f"workflow {filespec} has an unparseable 'on:' block")
        return
    for ev in events:
        if ev not in on_val:
            have = ", ".join(sorted(str(k) for k in on_val)) or "(none)"
            findings.append(f"workflow {filespec} is NOT triggered by {ev!r} "
                            f"(its triggers: {have}) -- the claimed enforcement "
                            f"would not run on that event")
            continue
        cfg = on_val.get(ev)
        if cfg is None or not isinstance(cfg, dict):
            continue  # event with no filter: runs for every branch
        pats = cfg.get("branches")
        if isinstance(pats, list) and pats:
            if not any(fnmatch.fnmatchcase(branch, str(p)) for p in pats):
                findings.append(f"workflow {filespec} trigger {ev} filters branches "
                                f"{pats}, which exclude the enforcement branch "
                                f"{branch!r} (disabled for the required branch)")
        pats = cfg.get("branches-ignore")
        if isinstance(pats, list) and pats:
            if any(fnmatch.fnmatchcase(branch, str(p)) for p in pats):
                findings.append(f"workflow {filespec} trigger {ev} ignores the "
                                f"enforcement branch {branch!r} "
                                f"(branches-ignore {pats})")


def resolve_refs(value: str, repo: str, require_live: bool = False,
                  enforcement_branch: str = "main"):
    """Validate each named reference in a cell. Returns (findings, none_declared).

    Round-19 EC-26 fix: when require_live is set (rule status ACTIVE), a
    `workflow:FILE#JOB` reference must be a LIVE workflow -- unpinned (no
    @<ref>), located under `.github/workflows/`, and present in the working
    tree. GitHub Actions runs only `.github/workflows/`; a `proposals/ci/...`
    draft is inert no matter how correct its contents are, and an @<ref> pin
    can name a branch whose version has the job while the live working-tree
    version does not. PROPOSED rules may still cite inert drafts (their
    activation path is the reviewed copy into the live path).

    Round-20 trigger-liveness fix: an ACTIVE rule's workflow reference must
    also declare its required events (`~events=push+pull_request` suffix) and
    the live workflow must be triggered by each declared event with branch
    filters covering the enforcement branch. Without the annotation the
    reference fails -- path + job existence alone no longer proves liveness.
    """
    findings = []
    items = [v.strip() for v in ITEM_SPLIT.split(value) if v.strip()]
    if not items:
        return [f"empty enforcement/test cell (fail closed)"], False
    none_declared = False
    for item in items:
        m = ITEM_RE.match(item)
        if not m:
            findings.append(f"unparseable reference {item!r}")
            continue
        if item.upper() == "NONE":
            none_declared = True
            continue
        if item.startswith("AT-"):
            continue  # registry membership validated separately
        kind, target = m.group(2), m.group(3).strip()
        if kind in ("script", "procedural"):
            if not _path_exists(repo, target):
                findings.append(f"named enforcement file does not exist: {target}")
        elif kind == "workflow":
            base, _, ev_spec = target.partition("~events=")
            events = ([e.strip() for e in ev_spec.split("+") if e.strip()]
                      if ev_spec else None)
            filespec_all, _, job = base.rpartition("#")
            filespec, _, ref = (filespec_all.rpartition("@")
                                if "@" in filespec_all else (filespec_all, "", ""))
            job = job.strip()
            if not filespec or not job:
                findings.append(f"malformed workflow reference {item!r} (need file#JOB)")
                continue
            if events is not None:
                bad = [e for e in events if e not in TRIGGER_EVENTS]
                if bad:
                    findings.append(f"workflow reference {item!r} declares unknown "
                                    f"event token(s) {bad} (allowed: "
                                    f"{sorted(TRIGGER_EVENTS)})")
                    events = None
            if require_live:
                if ref:
                    findings.append(f"ACTIVE rule pins its workflow to @{ref}; the live "
                                    f"working-tree version must be cited: {item}")
                    continue
                if not filespec.startswith(".github/workflows/"):
                    findings.append(f"ACTIVE rule cites a workflow outside the live CI path "
                                    f".github/workflows/ (inert draft?): {item}")
                    continue
                if not os.path.exists(os.path.join(repo, filespec)):
                    findings.append(f"ACTIVE rule's workflow is absent from the working tree "
                                    f"(not live on this candidate): {item}")
                    continue
                if events is None:
                    findings.append(f"ACTIVE rule's workflow reference must declare the "
                                    f"events its enforcement claims (append "
                                    f"'~events=push+pull_request'); path + job existence "
                                    f"alone does not prove liveness: {item}")
                    continue
            doc = _yaml_doc(repo, filespec_all)
            if doc is None:
                findings.append(f"named workflow file does not exist: {filespec}")
            elif yaml is None:
                findings.append("PyYAML unavailable; refusing to guess (fail closed)")
            else:
                jobs = doc.get("jobs") or {}
                if not any(str(k) == job or (isinstance(v, dict) and str(v.get("name")) == job)
                           for k, v in jobs.items()):
                    findings.append(f"workflow job not found: {job!r} in {filespec}")
                elif events:
                    _validate_triggers(filespec, doc, events,
                                       enforcement_branch, findings)
        elif kind == "pytest":
            node_file = target.split("::", 1)[0]
            if not _path_exists(repo, node_file):
                findings.append(f"pytest node file does not exist: {node_file}")
    return findings, none_declared


def parse_index(path: str):
    """Return (meta, rules). Fail closed on malformed structure."""
    text = open(path, encoding="utf-8").read()
    if BEGIN not in text or END not in text:
        raise ValueError(f"rule inventory block missing from {path}")
    block = text.split(BEGIN, 1)[1].split(END, 1)[0].strip("\n")
    meta, rules, seen_ids, in_table = {}, [], [], False
    for raw in block.splitlines():
        line = raw.strip()
        if not line:
            continue
        m = KEYVAL_RE.match(line)
        if m and not line.startswith("|"):
            meta[m.group(1).lower()] = m.group(2).strip()
            continue
        if line.startswith("|"):
            cells = [c.strip() for c in line.strip("|").split("|")]
            if all(re.fullmatch(r"-{2,}", c or "-") for c in cells):
                continue  # separator row
            if len(cells) != 4:
                raise ValueError(f"malformed inventory row (need 4 columns): {raw!r}")
            if not in_table:
                if [c.lower() for c in cells] != ["rule id", "status", "enforcement point(s)", "acceptance test(s)"]:
                    raise ValueError(f"unexpected inventory header: {raw!r}")
                in_table = True
                continue
            rid, status, enf, tests = cells
            if rid in seen_ids:
                raise ValueError(f"duplicate rule ID in inventory: {rid}")
            seen_ids.append(rid)
            if status not in VALID_STATUS:
                raise ValueError(f"unknown status token {status!r} for {rid} (allowed: {sorted(VALID_STATUS)})")
            rules.append({"id": rid, "status": status, "enforcement": enf, "tests": tests})
    if not in_table:
        raise ValueError("inventory block contains no rule table")
    if not rules:
        raise ValueError("inventory rule set is empty (fail closed)")
    ids = [r["id"] for r in rules]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate rule IDs (fail closed)")
    return meta, rules


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", default="docs/directives/INDEX.md")
    ap.add_argument("--note", default=None, help="override the canonical note path from INDEX")
    ap.add_argument("--repo-root", default=".")
    ap.add_argument("--registry", default="scripts/ci/tests/ACCEPTANCE_TESTS.json")
    ap.add_argument("--include-proposed", action="store_true",
                    help="also require PROPOSED rules to be mapped (self-test mode for this branch)")
    ap.add_argument("--require-strings", default="",
                    help="comma-separated strings that must appear in the note (series ID, base SHA, ...)")
    ap.add_argument("--enforcement-branch", default="main",
                    help="branch the ACTIVE-rule workflows must actually trigger for "
                         "(default: main)")
    args = ap.parse_args()
    repo = os.path.abspath(args.repo_root)
    findings = []

    try:
        meta, rules = parse_index(os.path.join(repo, args.index))
    except (OSError, ValueError) as e:
        print(f"FAIL: cannot derive rule inventory from INDEX ({e}); fail closed")
        return 1

    note_rel = args.note or meta.get("reading-note-path")
    if not note_rel:
        print("FAIL: INDEX declares no reading-note-path and no --note override; fail closed")
        return 1
    note_path = os.path.join(repo, note_rel)
    if not os.path.exists(note_path):
        print(f"FAIL: reading note missing at canonical path {note_rel}; fail closed")
        return 1
    note = open(note_path, encoding="utf-8").read()

    blobs = [b.strip() for b in meta.get("pinned-blobs", "").split(",") if b.strip()]
    if not blobs:
        findings.append("INDEX declares no pinned-blobs; fail closed")
    for b in blobs:
        if b not in note:
            findings.append(f"note does not cite pinned blob {b}")

    registry_path = os.path.join(repo, args.registry)
    try:
        registry = json.load(open(registry_path, encoding="utf-8"))["tests"]
    except (OSError, ValueError, KeyError) as e:
        print(f"FAIL: acceptance-test registry unusable ({e}); fail closed")
        return 1

    required = [r for r in rules if r["status"] == "ACTIVE" or args.include_proposed]
    if not required:
        print("FAIL: derived required rule set is empty; fail closed")
        return 1

    for r in required:
        require_live = r["status"] == "ACTIVE"
        for cell, kind in ((r["enforcement"], "enforcement"), (r["tests"], "test")):
            refs, none_declared = resolve_refs(cell, repo, require_live=require_live,
                                               enforcement_branch=args.enforcement_branch)
            for f in refs:
                findings.append(f"{r['id']} {kind}: {f}")
            if none_declared:
                print(f"WARN: {r['id']} declares NONE for {kind} (declared absence, not hidden)")
        # test IDs must exist in the registry
        for tok in re.split(r"[,;\s]+", r["tests"]):
            tok = tok.strip()
            if re.fullmatch(r"AT-\d+", tok) and tok not in registry:
                findings.append(f"{r['id']}: acceptance test {tok} not in registry")
        # the note must map the rule
        mapped = any((rid in line and ("->" in line or "|" in line))
                     for line in note.splitlines() for rid in [r["id"]])
        if not mapped:
            findings.append(f"note does not map required rule {r['id']}")

    for s in [x.strip() for x in args.require_strings.split(",") if x.strip()]:
        if s not in note:
            findings.append(f"note does not contain required string {s!r}")

    if findings:
        print("\n".join(sorted(set(findings))))
        print(f"FAIL: {len(set(findings))} finding(s); required set = "
              f"{len(required)} rule(s) from {args.index}")
        return 1
    print(f"PASS: {note_rel} cites {len(blobs)} pinned blob(s) and maps all {len(required)} "
          f"required rule(s) derived from {args.index}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
