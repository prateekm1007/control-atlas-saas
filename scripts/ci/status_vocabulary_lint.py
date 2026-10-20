#!/usr/bin/env python3
"""EC-22 status vocabulary lint (inert until adopted).

Vocabulary: NOT_RUN | RUN_FAIL | RUN_PASS(run <id>) | PUSHED | MERGED.
Flags: RUN_PASS without a run id on the same line; a percentage claim on a line without a run id;
"clean scan" or "pre-existing" without a run id / 7+hex citation on the same line.
"""
import argparse, re, sys

RUN_ID = re.compile(r"(?:run_id=|run\s)(\d+)", re.I)
PCT = re.compile(r"\d+(?:\.\d+)?\s?%")
HEXISH = re.compile(r"\b[0-9a-f]{7,40}\b", re.I)
BAD_PHRASES = ("clean scan", "pre-existing")

def check_line(path: str, n: int, line: str):
    out = []
    if "RUN_PASS" in line and not RUN_ID.search(line):
        out.append(f"FAIL: {path}:{n}: RUN_PASS without run id")
    if PCT.search(line) and not RUN_ID.search(line):
        out.append(f"FAIL: {path}:{n}: percentage claim without a run id")
    low = line.lower()
    for ph in BAD_PHRASES:
        if ph in low and not (RUN_ID.search(line) or HEXISH.search(line)):
            out.append(f"FAIL: {path}:{n}: '{ph}' without run id or SHA citation")
    return out

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--report-only", action="store_true")
    a = ap.parse_args()
    findings = []
    for p in a.files:
        try:
            with open(p, encoding="utf-8") as fh:
                for n, line in enumerate(fh, 1):
                    findings += check_line(p, n, line.rstrip("\n"))
        except OSError as e:
            findings.append(f"FAIL: {p}: unreadable ({e})")
    print("\n".join(findings) if findings else "PASS: status vocabulary clean")
    return 0 if (a.report_only or not findings) else 1

if __name__ == "__main__":
    sys.exit(main())
