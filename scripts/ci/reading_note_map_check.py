#!/usr/bin/env python3
"""EC-26 reading-note mapping check (inert until adopted).

Verifies the note (a) cites each pinned blob SHA prefix and (b) maps every in-force rule ID to a
check (a line containing the ID and '->', or a table row containing the ID and a check reference).
"""
import argparse, re, sys

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--note", required=True)
    ap.add_argument("--blobs", default="", help="comma-separated pinned blob SHA prefixes")
    ap.add_argument("--rules", default="", help="comma-separated in-force rule IDs (e.g. EC-16,EC-17)")
    ap.add_argument("--rules-file", help="file with one rule ID per line (overrides --rules)")
    a = ap.parse_args()
    text = open(a.note, encoding="utf-8").read()
    findings = []
    for blob in [b.strip() for b in a.blobs.split(",") if b.strip()]:
        if blob not in text:
            findings.append(f"FAIL: EC-26 note does not cite pinned blob {blob}")
    rules = ([r.strip() for r in a.rules.split(",") if r.strip()] if not a.rules_file
             else [l.strip() for l in open(a.rules_file, encoding="utf-8") if l.strip()])
    for rule in rules:
        mapped = False
        for line in text.splitlines():
            if rule in line and ("->" in line or "|" in line):
                mapped = True
                break
        if not mapped:
            findings.append(f"FAIL: EC-26 rule {rule} not mapped to a check in {a.note}")
    print("\n".join(findings) if findings else f"PASS: {a.note} cites pinned blobs and maps all in-force rules")
    return 0 if not findings else 1

if __name__ == "__main__":
    sys.exit(main())
