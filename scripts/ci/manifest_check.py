#!/usr/bin/env python3
"""Governed-artifact manifest checker for the EC-16..EC-26 package (Round-18).

Direction 3C (auditor, 2026-10-20): the root `SHA256SUMS` manifest must cover the
COMPLETE intended package, be regenerated AFTER all content edits, and be verified
in a clean checkout. This tool defines the governed-artifact scope, verifies the
manifest against it, and can regenerate it.

Modes:
  manifest_check.py --verify            (default) verify the manifest
  manifest_check.py --update            regenerate the manifest over the current tree
  manifest_check.py --verify --manifest PATH --repo-root DIR

Checks performed by --verify (all fail closed, exit 1):
  1. every manifest line is well-formed: "<64 lowercase hex>  <path>" (two spaces);
  2. no duplicate paths and no duplicate hash+path lines;
  3. every listed path EXISTS and its actual SHA-256 matches the recorded digest
     (stale hash -> FAIL);
  4. every in-scope file present in the tree is listed (unlisted in-scope file -> FAIL);
  5. every REQUIRED scope group has at least one file;
  6. the manifest file itself is never listed (it cannot cover itself).
  7. every listed path is inside the governed scope (no out-of-scope entries -> FAIL).

Scope definition (single source of truth — edit here and nowhere else):
  REQUIRED groups (at least one file each must exist and be listed):
    docs/constitution/amendments/**   the PROPOSED amendment documents
    docs/directives/INDEX.md          the canonical directive index
    docs/evidence/ec16-ec26/**        the EC package evidence (receipts, logs, summaries)
    docs/evidence/push-log/**         push log, schema, incidents, reconciliation
    docs/evidence/reading-note/READING_NOTE_LATEST.md
    proposals/ci/**                   inert workflow drafts + template + move-map
    proposals/remediation/**          historical remediation fragments + fixtures
    scripts/ci/*.py                   checker scripts (enforcement points)
    scripts/ci/tests/*.py|*.json      acceptance suite + registry + runner
  OPTIONAL groups (must be listed IF a matching file exists in the tree):
    .gitleaksignore                   per-fingerprint secret-scan dispositions

Exit: 0 clean; 1 findings; 2 usage error.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys

MANIFEST_DEFAULT = "SHA256SUMS"

# (pattern, kind) — kind: "required" | "optional". Patterns use repo-relative
# glob-style suffixes implemented by simple prefix/dir matching below.
SCOPE = [
    ("docs/constitution/amendments/", "required"),
    ("docs/directives/INDEX.md", "required"),
    ("docs/evidence/ec16-ec26/", "required"),
    ("docs/evidence/push-log/", "required"),
    ("docs/evidence/reading-note/READING_NOTE_LATEST.md", "required"),
    ("proposals/ci/", "required"),
    ("proposals/remediation/", "required"),
    ("scripts/ci/*.py", "required"),
    ("scripts/ci/tests/", "required"),
    (".gitleaksignore", "optional"),
]


def iter_files(repo: str):
    for root, dirs, files in os.walk(repo):
        dirs[:] = [d for d in dirs if d != ".git"]
        for f in files:
            yield os.path.relpath(os.path.join(root, f), repo).replace(os.sep, "/")


def in_scope(path: str, pattern: str) -> bool:
    if pattern.endswith("/"):
        return path.startswith(pattern)
    if "*" in pattern:
        import fnmatch
        return fnmatch.fnmatch(path, pattern)
    return path == pattern


def scope_of(path: str):
    return [(p, k) for p, k in SCOPE if in_scope(path, p)]


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_manifest(repo: str, manifest: str):
    entries, errors = {}, []
    for n, raw in enumerate(open(os.path.join(repo, manifest), encoding="utf-8"), 1):
        line = raw.rstrip("\n")
        if not line.strip():
            continue
        parts = line.split("  ", 1)
        if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
            errors.append(f"manifest line {n}: malformed (expected '<64hex>  <path>'): {line!r}")
            continue
        digest, rel = parts[0].strip(), parts[1].strip()
        rel = rel.replace(os.sep, "/").lstrip("./")
        if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            errors.append(f"manifest line {n}: digest is not 64 lowercase hex: {digest!r}")
            continue
        if rel in entries:
            errors.append(f"manifest line {n}: duplicate path listed twice: {rel}")
        entries[rel] = digest
    return entries, errors


def verify(repo: str, manifest: str) -> int:
    findings: list[str] = []
    entries, errors = parse_manifest(repo, manifest)
    findings += errors

    if manifest.replace(os.sep, "/").lstrip("./") in entries:
        findings.append("the manifest file itself is listed; a manifest cannot cover itself")

    # listed entries must exist, hash-match and be in scope
    for rel, digest in sorted(entries.items()):
        p = os.path.join(repo, rel)
        if not os.path.isfile(p):
            findings.append(f"listed but MISSING from the tree: {rel}")
            continue
        actual = sha256_file(p)
        if actual != digest:
            findings.append(f"STALE hash for {rel}: manifest {digest[:12]}… != actual {actual[:12]}…")
        if not scope_of(rel):
            findings.append(f"listed path is OUTSIDE the governed scope: {rel}")

    # in-scope files must be listed
    listed = set(entries)
    seen_groups: dict[str, int] = {}
    for rel in iter_files(repo):
        groups = scope_of(rel)
        if not groups:
            continue
        seen_groups.update({p: seen_groups.get(p, 0) + 1 for p, _ in groups})
        if rel not in listed:
            findings.append(f"in-scope file NOT listed in {manifest}: {rel}")

    for pattern, kind in SCOPE:
        if kind == "required" and seen_groups.get(pattern, 0) == 0:
            findings.append(f"required scope group has no files in the tree: {pattern}")

    if findings:
        print("\n".join(f"FAIL: {f}" for f in findings))
        print(f"MANIFEST STATUS: INVALID ({len(findings)} finding(s))")
        return 1
    print(f"MANIFEST STATUS: COMPLETE — {len(entries)} entries cover the governed package; "
          f"all digests verified fresh against the working tree")
    return 0


def update(repo: str, manifest: str) -> int:
    entries, _ = parse_manifest(repo, manifest)
    lines = []
    for rel in sorted(iter_files(repo)):
        if rel == manifest.replace(os.sep, "/").lstrip("./"):
            continue
        if not scope_of(rel):
            continue
        lines.append(f"{sha256_file(os.path.join(repo, rel))}  {rel}")
    if not lines:
        print("FAIL: scope produced zero entries; refusing to write an empty manifest")
        return 1
    body = "\n".join(lines) + "\n"
    p = os.path.join(repo, manifest)
    old = open(p, encoding="utf-8").read() if os.path.exists(p) else ""
    open(p, "w", encoding="utf-8").write(body)
    print(f"MANIFEST UPDATED: {len(lines)} entries ({'changed' if old != body else 'unchanged'})")
    return verify(repo, manifest)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true", help="verify the manifest (default)")
    ap.add_argument("--update", action="store_true", help="regenerate the manifest")
    ap.add_argument("--manifest", default=MANIFEST_DEFAULT)
    ap.add_argument("--repo-root", default=".")
    a = ap.parse_args()
    repo = os.path.abspath(a.repo_root)
    if not os.path.isdir(repo):
        print(f"FAIL: repo root does not exist: {repo}")
        return 2
    if a.update:
        return update(repo, a.manifest)
    return verify(repo, a.manifest)


if __name__ == "__main__":
    sys.exit(main())
