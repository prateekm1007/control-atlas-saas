#!/usr/bin/env python3
"""EC-21 declared base (inert until adopted). Verifies HEAD descends from the declared base SHA."""
import argparse, subprocess, sys

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--head", default="HEAD")
    a = ap.parse_args()
    r = subprocess.run(["git", "cat-file", "-e", f"{a.base}^{{commit}}"], capture_output=True)
    if r.returncode != 0:
        print(f"FAIL: EC-21 declared base not a commit on this remote: {a.base}")
        return 1
    r = subprocess.run(["git", "merge-base", "--is-ancestor", a.base, a.head], capture_output=True)
    if r.returncode != 0:
        print(f"FAIL: EC-21 HEAD does not descend from declared base {a.base}")
        return 1
    print(f"PASS: EC-21 HEAD descends from declared base {a.base}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
