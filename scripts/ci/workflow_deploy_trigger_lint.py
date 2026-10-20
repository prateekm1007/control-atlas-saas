#!/usr/bin/env python3
"""EC-25 no-deploy-on-push lint (inert until adopted).

Scans .github/workflows/*.yml|*.yaml. A workflow FAILS when it contains a job whose id or
display name matches the deploy regex AND any of the forbidden triggers (push, workflow_run).
Exit 0 when clean; 1 when findings (suppress with --report-only for evidence capture on a
baseline that already violates, e.g. main).
"""
import argparse, os, re, sys
try:
    import yaml
except ImportError:
    yaml = None

DEPLOY_RE = re.compile(r"\b(deploy|publish|rollback|release)\b", re.I)
FORBIDDEN = {"push", "workflow_run"}

def findings_for(path: str):
    with open(path, encoding="utf-8") as fh:
        doc = yaml.safe_load(fh)
    if not isinstance(doc, dict):
        return []
    on = doc.get("on", doc.get(True, {})) or {}
    if isinstance(on, list):
        triggers = set(on)
    elif isinstance(on, dict):
        triggers = set(on.keys())
    else:
        triggers = set()
    bad_triggers = sorted(triggers & FORBIDDEN)
    if not bad_triggers:
        return []
    out = []
    for job_id, job in (doc.get("jobs") or {}).items():
        names = [str(job_id)]
        if isinstance(job, dict) and job.get("name"):
            names.append(str(job["name"]))
        if any(DEPLOY_RE.search(n) for n in names):
            out.append(f"FAIL: {os.path.basename(path)}: job '{job_id}' matches deploy regex "
                       f"and workflow triggers on {bad_triggers}")
    return out

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("path", nargs="?", default=".github/workflows")
    ap.add_argument("--report-only", action="store_true", help="always exit 0 (evidence capture)")
    a = ap.parse_args()
    if yaml is None:
        print("FAIL: PyYAML unavailable; refusing to guess (fail closed)")
        return 1
    all_findings, scanned = [], 0
    for fn in sorted(os.listdir(a.path)):
        if fn.endswith((".yml", ".yaml")):
            scanned += 1
            all_findings += findings_for(os.path.join(a.path, fn))
    print(f"# scanned {scanned} workflow file(s) under {a.path}")
    print("\n".join(all_findings) if all_findings else "PASS: no deploy job on push/workflow_run triggers")
    if all_findings and not a.report_only:
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())
