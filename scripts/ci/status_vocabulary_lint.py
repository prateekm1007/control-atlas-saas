#!/usr/bin/env python3
"""EC-22 status vocabulary lint (inert until adopted). Auditor revision 2026-10-20 (item 3.4).

Machine-checkable status grammar (one status per line):

    STATUS: NOT_RUN | RUN_FAIL | PUSHED | MERGED
    STATUS: RUN_PASS(run_id=<digits>, sha=<40hex>)          # needs a verifiable CI run
    STATUS: RUN_PASS(run_id=<digits>, sha=<40hex>) pct=<N>  # pct must equal computed progress

Rules (all fail closed):
  * Unknown status tokens are REJECTED (the vocabulary is exactly
    {NOT_RUN, RUN_FAIL, RUN_PASS, PUSHED, MERGED}; anything else fails).
  * RUN_PASS without a run_id AND a full 40-hex sha on the same line fails.
  * With --runs-manifest (JSON {"runs": {"<id>": {"head_sha": "...", "conclusion": "..."}}}):
      - a cited run_id absent from the manifest fails (fabricated / unverifiable);
      - a run whose manifest head_sha differs from the cited sha fails (wrong commit);
      - a run whose conclusion is not "success" fails;
      - WITHOUT a manifest, every RUN_PASS line fails (no unverified progress).
  * Any percentage claim (pct= or "N%") must equal the progress COMPUTED from the accepted
    RUN_PASS records in the same file (accepted RUN_PASS / total STATUS lines, floored);
    a percentage unsupported by run records fails.
  * "clean scan" without a scan=<citation> on the same line fails.
  * "pre-existing" without a base=<40hex> citation on the same line fails.

Exit: 0 clean (or --report-only), 1 findings, 2 usage.
"""
import argparse
import json
import re
import sys

# Round-19: PR-body Status-summary lines carry the same grammar (the PR template
# declares the vocabulary for them); both line forms are linted identically.
STATUS_LINE = re.compile(r"^\s*(?:STATUS|Status-summary):\s*(.+?)\s*$")
TOKEN = re.compile(r"^[A-Z_]+$")
RUNPASS_ARGS = re.compile(r"^RUN_PASS\s*\((.*)\)\s*$")
RUN_ID = re.compile(r"run_id=(\d+)")
SHA_ARG = re.compile(r"sha=([0-9a-f]{40})")
PCT = re.compile(r"(?:pct=|(\b)\s*)(\d+(?:\.\d+)?)\s?%")
HEX40 = re.compile(r"\b[0-9a-f]{40}\b", re.I)
BAD_PHRASES = {"clean scan": "scan=", "pre-existing": "base="}
ALLOWED = {"NOT_RUN", "RUN_FAIL", "RUN_PASS", "PUSHED", "MERGED"}


def lint_path(path: str, manifest: dict | None):
    findings = []
    try:
        lines = open(path, encoding="utf-8").read().splitlines()
    except OSError as e:
        return [f"{path}: unreadable ({e})"], {"total": 0, "accepted": 0}
    total_status, accepted_pass = 0, 0
    pct_claims = []
    for n, line in enumerate(lines, 1):
        low = line.lower()
        m = STATUS_LINE.match(line)
        if m:
            total_status += 1
            expr = m.group(1)
            bare = expr.split("(", 1)[0].strip()
            rest = expr
            if bare not in ALLOWED:
                findings.append(f"{path}:{n}: unknown status token {bare!r} "
                                f"(allowed: {sorted(ALLOWED)})")
            if bare == "RUN_PASS":
                rid = RUN_ID.search(expr)
                sha = SHA_ARG.search(expr)
                if not (rid and sha):
                    findings.append(f"{path}:{n}: RUN_PASS without run_id=<digits> and sha=<40hex> "
                                    f"on the same line")
                elif manifest is None:
                    findings.append(f"{path}:{n}: RUN_PASS accepted only against a --runs-manifest "
                                    f"(no unverified progress)")
                else:
                    rec = manifest.get("runs", {}).get(rid.group(1))
                    if rec is None:
                        findings.append(f"{path}:{n}: run_id={rid.group(1)} not found in the runs "
                                        f"manifest (fabricated or unverifiable)")
                    elif rec.get("head_sha") != sha.group(1):
                        findings.append(f"{path}:{n}: run_id={rid.group(1)} is for commit "
                                        f"{str(rec.get('head_sha'))[:12]}…, not the cited sha "
                                        f"{sha.group(1)[:12]}… (wrong-SHA citation)")
                    elif rec.get("conclusion") != "success":
                        findings.append(f"{path}:{n}: run_id={rid.group(1)} conclusion is "
                                        f"{rec.get('conclusion')!r}, not success")
                    else:
                        accepted_pass += 1
            pm = re.search(r"pct=(\d+(?:\.\d+)?)", expr)
            if pm:
                pct_claims.append((n, float(pm.group(1))))
        # phrase rules apply to every line, status or prose
        for phrase, cite in BAD_PHRASES.items():
            if phrase in low and cite not in line:
                findings.append(f"{path}:{n}: '{phrase}' without a {cite}<citation> on the same line")
        for pm in re.finditer(r"(\d+(?:\.\d+)?)\s?%", line):
            if "pct=" not in line:
                pct_claims.append((n, float(pm.group(1))))
    return findings, {"total": total_status, "accepted": accepted_pass, "pct_claims": pct_claims}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--runs-manifest", help="JSON mapping run IDs to {head_sha, conclusion}")
    ap.add_argument("--report-only", action="store_true")
    a = ap.parse_args()
    manifest = None
    if a.runs_manifest:
        try:
            manifest = json.load(open(a.runs_manifest, encoding="utf-8"))
        except (OSError, ValueError) as e:
            print(f"FAIL: runs manifest unusable ({e}); fail closed")
            return 1
    all_findings = []
    for p in a.files:
        f, stats = lint_path(p, manifest)
        all_findings += f
        computed = 0
        if stats["total"]:
            computed = int(stats["accepted"] * 100 // stats["total"])
        for n, claimed in stats.get("pct_claims", []):
            if claimed != computed:
                all_findings.append(
                    f"{p}:{n}: percentage claim {claimed:g}% is not supported by run records "
                    f"(computed progress: {computed}% from {stats['accepted']}/"
                    f"{stats['total']} accepted RUN_PASS)")
    if all_findings:
        print("\n".join(all_findings))
        print(f"FAIL: {len(all_findings)} finding(s)")
        return 0 if a.report_only else 1
    print("PASS: status vocabulary clean")
    return 0


if __name__ == "__main__":
    sys.exit(main())
