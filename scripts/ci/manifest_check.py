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
  3. every digest MUST match the COMMITTED bytes (the HEAD tree) of the
     file's actual SHA-256 -- the manifest pins the committed artifact, so the
     check does not depend on the line endings a checkout materialises;
  4. the scope is derived from the SAME committed tree, and an index that
     carries staged changes is REJECTED (a staged tree is a different
     artifact and must never be reported as the committed one);
  5. the scope MUST be non-empty and SHA256SUMS MUST exist and be parseable.

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


def committed_tree(repo: Path, rev: str = "HEAD") -> list[str]:
    """Every path in the COMMITTED tree of `rev`."""
    r = subprocess.run(["git", "-C", str(repo), "ls-tree", "-r", "--name-only", rev],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise subprocess.CalledProcessError(r.returncode, r.args, r.stderr)
    return [l.strip() for l in r.stdout.splitlines() if l.strip()]


def _seg_match(pattern: str, text: str) -> bool:
    """One path segment against one segment pattern ('*' and '?' do not cross
    a '/'; '**' is handled at the path level)."""
    rx = ""
    for c in pattern:
        if c == "*":
            rx += "[^/]*"
        elif c == "?":
            rx += "[^/]"
        else:
            rx += re.escape(c)
    return re.fullmatch(rx, text) is not None


def path_matches(pattern: str, path: str) -> bool:
    """Does `path` match the glob `pattern` with git's path semantics?

    '*' matches inside one segment, '**' matches any number of segments.
    """
    pat = [s for s in pattern.split("/") if s != ""]
    seg = [s for s in path.split("/") if s != ""]

    def walk(i: int, j: int) -> bool:
        if i == len(pat):
            return j == len(seg)
        if pat[i] == "**":
            # '**' consumes zero or more whole segments
            if walk(i + 1, j):
                return True
            return j < len(seg) and walk(i, j + 1)
        if j >= len(seg):
            return False
        if not _seg_match(pat[i], seg[j]):
            return False
        return walk(i + 1, j + 1)

    return walk(0, 0)


def scope_files(repo: Path, rev: str = "HEAD") -> set[str]:
    """The package scope in the COMMITTED tree of `rev`.

    The scope is derived from the same committed tree the digests are taken
    from, so the two can never describe different states.
    """
    tree = committed_tree(repo, rev)
    return {p for p in tree if any(path_matches(g, p) for g in SCOPE_GLOBS)}


def index_matches_head(repo: Path) -> tuple[bool, str]:
    """Is the index identical to the commit under validation?

    The manifest validates a COMMITTED tree. If the index carries staged
    changes, the bytes a checker would read from the index are NOT the bytes
    of that commit: a manifest generated over the index would then be
    reported as "the committed artifact is verified" while the commit under
    validation says something else. Such a state is REJECTED instead.
    """
    r = subprocess.run(["git", "-C", str(repo), "diff", "--cached", "--name-only",
                        "HEAD"], capture_output=True, text=True)
    if r.returncode != 0:
        return False, "the index could not be compared with the commit under validation"
    changed = [l.strip() for l in r.stdout.splitlines() if l.strip()]
    if changed:
        return False, (f"{len(changed)} path(s) are staged and differ from the commit "
                       f"under validation (first: {changed[0]})")
    return True, ""


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


def committed_bytes(repo: Path, rel: str, rev: str = "HEAD") -> bytes:
    """The COMMITTED bytes of a tracked file at `rev`.

    NOT the index entry and NOT the materialised working tree: the manifest
    is about the committed artifact. Reading the index (`git cat-file blob
    :<path>`) would let a staged, uncommitted tree satisfy a check that
    claims to validate the commit.
    """
    r = subprocess.run(["git", "-C", str(repo), "cat-file", "blob", f"{rev}:{rel}"],
                       capture_output=True)
    if r.returncode == 0:
        return r.stdout
    return b""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", default=".")
    ap.add_argument("--manifest", default="SHA256SUMS")
    a = ap.parse_args()
    repo = Path(a.repo_root).resolve()
    findings: list[str] = []

    try:
        scope = scope_files(repo, "HEAD")
    except subprocess.CalledProcessError as e:
        print(f"FAIL: the committed tree could not be read ({e}); fail closed")
        return 2
    if not scope:
        print("FAIL: derived package scope is empty; fail closed")
        return 1
    # The index must not carry staged changes: the check validates the
    # COMMITTED tree, and a staged tree is a different artifact (Round-25).
    same, why = index_matches_head(repo)
    if not same:
        print(f"FAIL: the index differs from the commit under validation ({why}) -- "
              f"the manifest validates the COMMITTED tree, so a staged tree is rejected "
              f"(fail closed)")
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
        data = committed_bytes(repo, f, "HEAD")
        if not data:
            findings.append(f"manifest entry is not a file in the committed tree: {f}")
            continue
        # The manifest pins the COMMITTED artifact bytes, not whatever a
        # particular checkout materialises. Hashing the working tree made the
        # check platform-dependent: on a checkout whose line endings differ
        # from the stored blobs (Windows CRLF against LF blobs) every digest
        # mismatched, so a correct manifest was reported as corrupt on one
        # platform and accepted on another. Hashing the index instead would
        # let a staged tree satisfy a check that claims to validate the
        # commit -- the index is rejected above. The committed bytes are the
        # artifact the manifest is about.
        h = hashlib.sha256(data).hexdigest()
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
