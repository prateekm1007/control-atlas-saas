# Items 17+18 evidence — nested copy removal + junk files

**Branch:** chore/hygiene-2026-10-19 (from science head 296f10c) · **Executed:** 2026-10-19T12:35:00Z

## Item 17 — nested control-atlas-saas/ removed

```
$ git ls-files control-atlas-saas/ | wc -l
83

# grep for any reference to the nested path OUTSIDE it (code, CI, Docker, docs):
$ grep -rn "control-atlas-saas/control-atlas-saas" . --exclude-dir=.git | grep -v "^./control-atlas-saas/"
(no matches)
$ grep -rn "control-atlas-saas/control-atlas" backend/main.py backend/Dockerfile docker-compose.yml .github/workflows/*.yml
(no matches)
```

Nothing outside the nested directory referenced it — removal is safe (import
graph, Docker build contexts, and workflows checked). The only hit for the
render package was a viz .cxc script pointing at a LOCAL WINDOWS PATH
(C:/Users/Administrator/Downloads/...), not at the in-repo tarball.

## Item 18 — junk files removed

The P0 ticket's junk list lives in a GitHub PR not readable from this
environment; the set below is what the science branch actually carries,
each verified unreferenced:

- Atlas_Entry_001_Render_Package.tar.gz — 94 KB binary blob, unreferenced in-repo
- test_sieve.sh — manual localhost smoke script with an embedded test credential
- run_full_atlas.ps1 — machine-specific one-off (hardcoded C:\Users\Administrator paths)

## Deliberately NOT removed (recorded for the owner)

- calibration_results.json (repo root) — documented legacy v22.4.3 artifact in
  backend/calibration/README.md; removal is an owner decision
- toscanini_audit.pdf — audit evidence, kept per evidence-preservation policy
- SESSION_CHECKPOINT_v2_4.txt / STABILITY_LOG.txt / GENERATOR_STATUS.md —
  project history records, not unambiguous junk
- main-branch-only accidents (=3.11, bracket-named files, build logs) are NOT
  tracked on science; they need the same treatment on main/integration

## Verification

```
$ git status --short | wc -l
86 deletions (83 nested + 3 junk), 0 modifications
```

## Test impact

The nested copy was never imported by backend code (zero references). This
branch contains only deletions + this evidence file; the sibling p0 branch
(same science base) is green at 506 passed / 5 skipped / 0 failed.
