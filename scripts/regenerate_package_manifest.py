#!/usr/bin/env python3
"""Deterministic regenerator for the repository package manifest (SHA256SUMS).

Reuses the exact scope policy enforced by ``scripts/ci/manifest_check.py``
(the EC-16..EC-26 governance/enforcement package scope) so the manifest can
never drift from what the checker requires. Output is deterministic: sorted by
relative POSIX path, one ``<sha256>  <path>`` line per scoped file. The file
's own digest is intentionally never included (self-integrity is provided by
git's object model and the clean-checkout verification).
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path


def _scope_files(repo: Path) -> set[str]:
    sys.path.insert(0, str(repo / "scripts" / "ci"))
    try:
        from manifest_check import scope_files as _sf
        return _sf(repo)
    finally:
        sys.path.pop(0)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root", default=".")
    ap.add_argument("--manifest", default="SHA256SUMS")
    a = ap.parse_args()
    repo = Path(a.repo_root).resolve()
    scope = sorted(_scope_files(repo))
    if not scope:
        print("FAIL: derived package scope is empty; abort (no manifest written)")
        return 1
    lines = []
    for rel in scope:
        fp = repo / rel
        if not fp.is_file():
            print(f"FAIL: scope file is not a regular file: {rel}")
            return 1
        digest = hashlib.sha256(fp.read_bytes()).hexdigest()
        lines.append(f"{digest}  {rel}")
    mpath = repo / a.manifest
    mpath.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"PASS: wrote {len(lines)} entries to {a.manifest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
