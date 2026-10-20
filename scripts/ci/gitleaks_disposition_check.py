#!/usr/bin/env python3
"""EC-23 disposition cross-check (Round-18).

Guarantees the per-fingerprint discipline is airtight:
  1. every fingerprint in `.gitleaksignore` has a row in
     docs/evidence/ec16-ec26/GITLEAKS_FINDINGS_DISPOSITION.md (no unexplained ignore);
  2. every fingerprint in the disposition file exists in `.gitleaksignore`
     (no stale disposition rows);
  3. every disposition row carries a classification from the approved vocabulary
     (FALSE-POSITIVE-*, DOCUMENTED-HISTORICAL-A-02) — anything else fails closed.
Exit 0 = consistent; 1 = findings.
"""
from __future__ import annotations
import os
import re
import sys

IGNORE = ".gitleaksignore"
DISPO = "docs/evidence/ec16-ec26/GITLEAKS_FINDINGS_DISPOSITION.md"
FP = re.compile(r"^[0-9a-f]{40}:[^:]+:[^:]+:\d+$")
VALID_CLASSES = {"FALSE-POSITIVE-PROVENANCE-DIGEST", "FALSE-POSITIVE-IDENTIFIER-LIST",
                 "DOCUMENTED-HISTORICAL-A-02"}


def main() -> int:
    findings = []
    ignore_fps = []
    if os.path.exists(IGNORE):
        for n, raw in enumerate(open(IGNORE, encoding="utf-8"), 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if not FP.match(line):
                findings.append(f"{IGNORE}:{n}: not a valid fingerprint line: {line!r}")
                continue
            ignore_fps.append(line)
    else:
        findings.append(f"{IGNORE} missing (fail closed)")

    dispo_fps = []
    dispo_classes = {}
    current = None
    if os.path.exists(DISPO):
        for raw in open(DISPO, encoding="utf-8"):
            line = raw.strip()
            if line.startswith("## ") and line[3:] in VALID_CLASSES:
                current = line[3:]
            m = re.match(r"^[\|\s]*`([0-9a-f]{40}:[^`]+)`", line)
            if m and FP.match(m.group(1)):
                dispo_fps.append(m.group(1))
                dispo_classes[m.group(1)] = current
    else:
        findings.append(f"{DISPO} missing (fail closed)")

    for fp in ignore_fps:
        if fp not in dispo_fps:
            findings.append(f"fingerprint in {IGNORE} has NO disposition row: {fp}")
        elif dispo_classes.get(fp) not in VALID_CLASSES:
            findings.append(f"fingerprint {fp} disposition class missing/invalid: "
                            f"{dispo_classes.get(fp)!r}")
    for fp in dispo_fps:
        if fp not in ignore_fps:
            findings.append(f"disposition row not present in {IGNORE} (stale): {fp}")

    if findings:
        print("\n".join(f"FAIL: {f}" for f in findings))
        print(f"DISPOSITION STATUS: INVALID ({len(findings)} finding(s))")
        return 1
    print(f"DISPOSITION STATUS: CONSISTENT — {len(ignore_fps)} fingerprint(s), every "
          f".gitleaksignore entry justified, every disposition row active")
    return 0


if __name__ == "__main__":
    sys.exit(main())
