#!/usr/bin/env python3
"""Package manifest coverage + integrity check (Round-19 direction 3C).

The root SHA256SUMS is the artifact manifest of the EC-16..EC-26 governance /
enforcement package. Direction 3C: define the manifest's exact scope, include
every intended package artifact, and FAIL when a required package file is
missing from the manifest or its digest is wrong.

Scope policy (exact, machine-enforced):

    docs/constitution/amendments/**        the amendment texts
    docs/directives/INDEX.md               the machine-readable rule inventory
    docs/evidence/ec16-ec26/**             package evidence
    docs/evidence/push-log/**              push log + schema + incident records
    docs/evidence/reading-note/*           the canonical reading note
    proposals/ci/**                        inert workflow drafts + template
    proposals/remediation/**               superseded remediation drafts + fixtures
    scripts/ci/**                          the checkers + acceptance tests

Deliberately OUT of scope (documented boundary): .github/workflows/ (live CI
infrastructure, validated by its own live jobs and the deploy-trigger lint),
docs/P0_STATUS.md and docs/decisions/ (repo status / register documents), and
SHA256SUMS itself (a manifest cannot contain its own digest; self-integrity
is provided by git's object model and the clean-checkout verification
recorded in the round evidence).

Rules enforced (all fail-closed):
  1. every tracked file matching the scope globs MUST appear in SHA256SUMS;
  2. every SHA256SUMS entry MUST be a scope file (no undocumented entries);
  3. every digest MUST match the file's actual SHA-256;
  4. the scope MUST be non-empty and SHA256SUMS MUST exist and be parseable.

Exit: 0 = manifest complete and correct; 1 = findings; 2 = usage/internal.
"""
import argparse
import hashlib
import re
import subprocess
import sys
from pathlib import Path

SCOPE_GLOBS = [
    "docs/constitution/amendments/**",
    "docs/directives/INDEX.md",
    "docs/evidence/ec16-ec26/**",
    "docs/evidence/push-log/**",
    "docs/evidence/reading-note/*",
    "proposals/ci/**",
    "proposals/remediation/**",
    "scripts/ci/**",
]
ENTRY_RE = re.compile(r"^([0-9a-f]{64})[ \t]+\*?(.+)$")


def tracked_files(repo: Path, glob: str) -> list[str]:
    r = subprocess.run(["git", "ls-files", "--", glob], capture_output=True,
                       text=True, cwd=repo, check=True)
    return [l.strip() for l in r.stdout.splitlines() if l.strip()]


def scope_files(repo: Path) -> set[str]:
    out: set[str] = set()
    for g in SCOPE_GLOBS:
        out.update(tracked_files(repo, g))
    return out


def load_manifest(path: Path) -> dict[str, str]:
    manifest: dict[str, str] = {}
    for n, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        m = ENTRY_RE.match(line)
        if not m:
            raise ValueError(f"unparseable manifest line {n}: {line!r}")
        digest, fpath = m.group(1), m.group(2).strip()
        if fpath in manifest:
            raise ValueError(f"duplicate manifest entry for {fpath}")
        manifest[fpath] = digest
    return manifest


def index_bytes(repo: Path, rel: str) -> bytes:
    """The COMMITTED (index) bytes of a tracked file.

    Falls back to the working tree when the path is not in the index (an
    untracked file can only be in the scope through a manual manifest entry,
    which the scope check rejects anyway).
    """
    r = subprocess.run(["git", "-C", str(repo), "cat-file", "blob", ":" + rel],
                       capture_output=True)
    if r.returncode == 0:
        return r.stdout
    fp = repo / rel
    return fp.read_bytes() if fp.is_file() else b""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", default=".")
    ap.add_argument("--manifest", default="SHA256SUMS")
    a = ap.parse_args()
    repo = Path(a.repo_root).resolve()
    findings: list[str] = []

    try:
        scope = scope_files(repo)
    except subprocess.CalledProcessError as e:
        print(f"FAIL: git ls-files unusable ({e}); fail closed")
        return 2
    if not scope:
        print("FAIL: derived package scope is empty; fail closed")
        return 1

    mpath = repo / a.manifest
    if not mpath.exists():
        print(f"FAIL: manifest {a.manifest} missing; fail closed")
        return 1
    try:
        manifest = load_manifest(mpath)
    except (OSError, ValueError) as e:
        print(f"FAIL: manifest unusable ({e}); fail closed")
        return 1

    missing = sorted(scope - set(manifest))
    for f in missing:
        findings.append(f"required package file missing from the manifest: {f}")
    orphans = sorted(set(manifest) - scope)
    for f in orphans:
        findings.append(f"manifest entry outside the defined package scope: {f}")
    for f in sorted(scope & set(manifest)):
        fp = repo / f
        if not fp.is_file():
            findings.append(f"manifest entry is not a file in the working tree: {f}")
            continue
        # The manifest pins the COMMITTED artifact bytes (the index entry),
        # not whatever a particular checkout materialises. Hashing the working
        # tree made the check platform-dependent: on a checkout whose line
        # endings differ from the stored blobs (Windows CRLF against LF blobs)
        # every digest mismatched, so a correct manifest was reported as
        # corrupt on one platform and accepted on another. The committed bytes
        # are the artifact the manifest is about.
        h = hashlib.sha256(index_bytes(repo, f)).hexdigest()
        if h != manifest[f]:
            findings.append(f"digest mismatch: {f} (manifest {manifest[f][:16]}... != committed {h[:16]}...)")

    if findings:
        print("\n".join(f"FAIL: {f}" for f in findings))
        print(f"FAIL: {len(findings)} finding(s); scope = {len(scope)} file(s), "
              f"manifest = {len(manifest)} entrie(s)")
        return 1
    print(f"PASS: manifest covers all {len(scope)} package-scope file(s) with correct digests")
    return 0


if __name__ == "__main__":
    sys.exit(main())
