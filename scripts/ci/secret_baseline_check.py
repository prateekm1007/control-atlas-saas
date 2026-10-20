#!/usr/bin/env python3
"""EC-23 all-refs scan baseline verification (Round-20 hardening).

Auditor finding (Round-20 direction 3): the previous all-refs classification
tolerated EVERY `generic-api-key` finding under `docs/` or `evidence/` as an
"evidence-digest false-positive class". That category exemption is not
sufficiently narrow: a newly introduced credential committed under those
directories would also pass.

Replacement policy (this module): the scan report is compared against an
EXPLICIT, REVIEWED ALLOWLIST of finding fingerprints
(`docs/evidence/ec16-ec26/secret_scan_baseline.json`), where every entry is
grounded in a documented source, classification, and reason. Rules:

  - every finding in the report MUST match an allowlisted fingerprint
    (unknown fingerprint => FAIL, wherever the file lives -- including
    docs/ and evidence/);
  - a finding whose fingerprint matches but whose rule/file contradict the
    allowlisted entry FAILS (a reused fingerprint is not blind trust);
  - the two historical H-1 findings (the former provider key, A-02) are a
    SEPARATE classification, disclosed on every run, and tolerated ONLY while
    the register keeps A-02 open (issuer-side revocation unverified);
  - the baseline itself is validated: unique fingerprints, required fields
    (source, classification, basis), known classification tokens -- an
    undocumented or malformed baseline entry FAILS (the allowlist is review
    material, not a rubber stamp);
  - duplicate report rows for the same fingerprint (the scanner reports the
    same commit once per reaching ref) are collapsed and counted honestly.

Exit: 0 = every finding classified, baseline well-formed; 1 = findings or
baseline defects; 2 = usage/internal error.
"""
import argparse
import json
import os
import re
import sys

REQUIRED_ENTRY_FIELDS = ("fingerprint", "rule", "file", "classification", "source", "basis")
CLASSIFICATIONS = {
    "DOCUMENTED-HISTORICAL-A-02",
    "FALSE-POSITIVE-PROVENANCE-DIGEST",
    "FALSE-POSITIVE-IDENTIFIER-LIST",
    "FALSE-POSITIVE-DOCUMENTATION-QUOTE",
    "FALSE-POSITIVE-PUSH-LOG-EVIDENCE-QUOTE",
}
FINGERPRINT_RE = re.compile(r"^[0-9a-f]{40}:\S+:\S+:\d+$")


def load_baseline(path: str):
    try:
        doc = json.load(open(path, encoding="utf-8"))
    except (OSError, ValueError) as e:
        return None, [f"baseline unusable ({e}); fail closed"]
    if not isinstance(doc, dict) or not isinstance(doc.get("entries"), list):
        return None, ["baseline must be a JSON object with an 'entries' array"]
    meta = doc.get("metadata")
    if not isinstance(meta, dict) or not meta.get("scanner_version") \
            or not meta.get("policy"):
        return None, ["baseline metadata must record scanner_version and policy"]
    findings, seen = [], {}
    for i, e in enumerate(doc["entries"]):
        missing = [k for k in REQUIRED_ENTRY_FIELDS if not e.get(k)]
        if missing:
            findings.append(f"baseline entry #{i} missing required field(s): {missing}")
            continue
        if e["classification"] not in CLASSIFICATIONS:
            findings.append(f"baseline entry {e['fingerprint'][:24]}... has unknown "
                            f"classification {e['classification']!r}")
        if not FINGERPRINT_RE.match(e["fingerprint"]):
            findings.append(f"baseline entry fingerprint malformed: {e['fingerprint']!r}")
        if e["fingerprint"] in seen:
            findings.append(f"baseline duplicate fingerprint: {e['fingerprint']}")
        seen[e["fingerprint"]] = e
    if "A-02" not in json.dumps(doc.get("a02_note", "")) and any(
            e.get("classification") == "DOCUMENTED-HISTORICAL-A-02"
            for e in doc["entries"]):
        findings.append("baseline carries H-1 entries but no a02_note recording the "
                        "open issuer-revocation blocker (A-02)")
    return doc, findings


def check(report_path: str, baseline_path: str) -> int:
    fails: list[str] = []
    try:
        report = json.load(open(report_path, encoding="utf-8"))
    except (OSError, ValueError) as e:
        print(f"FAIL: scan report unusable ({e}) -- scanner execution failure is "
              f"not a pass (fail closed)")
        return 1
    if not isinstance(report, list):
        print("FAIL: scan report is not a JSON array of findings")
        return 1
    baseline, baseline_fails = load_baseline(baseline_path)
    for f in baseline_fails:
        fails.append("baseline: " + f)
    if baseline is None:
        for f in fails:
            print("FAIL: " + f)
        return 1

    entries = {e["fingerprint"]: e for e in baseline["entries"]}
    rows_by_fp: dict[str, list] = {}
    for row in report:
        fp = row.get("Fingerprint")
        if not fp:
            fails.append("report row without Fingerprint: %s@%s"
                         % (row.get("RuleID"), row.get("File")))
            continue
        rows_by_fp.setdefault(fp, []).append(row)

    counts: dict[str, int] = {}
    unknown, contradicted = [], []
    for fp, rows in sorted(rows_by_fp.items()):
        entry = entries.get(fp)
        if entry is None:
            unknown.append((fp, rows[0]))
            continue
        for row in rows:
            if row.get("RuleID") != entry["rule"] or row.get("File") != entry["file"]:
                contradicted.append((fp, row, entry))
        counts[entry["classification"]] = counts.get(entry["classification"], 0) + 1

    for fp, row in unknown:
        fails.append("UNKNOWN FINGERPRINT (not in the reviewed allowlist): %s "
                     "[%s] %s:%s -- a new finding, wherever it lives, fails the scan"
                     % (fp, row.get("RuleID"), row.get("File"), row.get("StartLine")))
    for fp, row, entry in contradicted:
        fails.append("fingerprint %s matched but report says rule=%s file=%s while the "
                     "reviewed entry records rule=%s file=%s (contradiction fails)"
                     % (fp, row.get("RuleID"), row.get("File"),
                        entry["rule"], entry["file"]))

    print("all-refs scan: %d report row(s) / %d distinct fingerprint(s); baseline "
          "entries: %d" % (len(report), len(rows_by_fp), len(entries)))
    for cls in sorted(counts):
        print("  %s: %d fingerprint(s)" % (cls, counts[cls]))
    h1 = [e for e in baseline["entries"]
          if e["classification"] == "DOCUMENTED-HISTORICAL-A-02" and
          e["fingerprint"] in rows_by_fp]
    for e in h1:
        print("  KNOWN (A-02 open, issuer-side revocation unverified -- deployment "
              "stays blocked): %s" % e["fingerprint"])
    if baseline.get("a02_note"):
        print("  a02_note: %s" % baseline["a02_note"])
    if fails:
        for f in fails:
            print("FAIL: " + f)
        return 1
    print("PASS: every finding in the all-refs scan matches a reviewed, documented "
          "baseline fingerprint (no category exemption; unknown fingerprints fail)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--report", required=True, help="gitleaks JSON report (all-refs scan)")
    ap.add_argument("--baseline", required=True,
                    help="reviewed fingerprint allowlist (JSON)")
    a = ap.parse_args()
    if not os.path.exists(a.report):
        print(f"FAIL: report {a.report} does not exist (scanner produced no report; "
              f"fail closed)")
        return 1
    if not os.path.exists(a.baseline):
        print(f"FAIL: baseline {a.baseline} does not exist (no allowlist => no pass)")
        return 1
    return check(a.report, a.baseline)


if __name__ == "__main__":
    sys.exit(main())
