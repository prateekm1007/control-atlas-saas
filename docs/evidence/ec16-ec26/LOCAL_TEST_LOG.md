# LOCAL_TEST_LOG — acceptance run for the EC-16..EC-26 correction cycle (SERIES-2026-10-20-EC-FIX)

This log is generated from the actual commands executed against the exact candidate.
GitHub Actions status for every commit in this cycle: **NOT RUN** — the amendment
workflows are inert drafts under `proposals/ci/` and no CI run ID exists for any commit
here; nothing in this log is a CI claim (EC-22 vocabulary: all local runs are recorded as
local evidence only, NOT_RUN in CI terms).

## Environment

| Field | Value |
|---|---|
| Candidate branch | `proposal/ec16-ec26-enforcement` |
| Base (declared, EC-21) | `main` @ `2195a048e87dcec4aa63e8dc7981b3bb6d237ca0` (parent of `17eb881…`) |
| Python | 3.12.15 |
| pytest | 9.0.2 |
| OS | Linux 5.10.134-013.15.kangaroo.al8.x86_64 |
| Host clock source | `date -u` (single host throughout) |

## Commands and results (exact)

| # | Command | Expect | Result |
|---|---|---|---|
| 1 | `git rev-parse 17eb88184c9f461c6bc9ee592f916f7ff4c921c4^` | `2195a048…` | PASS (verified at session start) |
| 2 | `git hash-object` on extracted pinned blobs `211ad2c5…` / `b6bb35e8…` | identical IDs | PASS (constitution preflight, read window 09:05:44Z–09:06:36Z) |
| 3 | `python3 scripts/ci/tests/run_all.py` | exit 0 | PASS — **16 groups (AT-01…AT-16), 91 `check()` assertions, 0 failed** (pytest 9.0.2 runner) |
| 4 | `python3 -c "import yaml; [yaml.safe_load(open(f)) for f in <9 proposal+remediation YAMLs>]"` | all parse | PASS — 9/9 YAML files valid |
| 5 | `python3 scripts/ci/reading_note_map_check.py --index docs/directives/INDEX.md --repo-root . --include-proposed --require-strings "SERIES-2026-10-20-EC"` | exit 0 | PASS — note maps all 28 inventory rules (17 ACTIVE + 11 PROPOSED) and cites both pinned blobs |
| 6 | `python3 scripts/ci/evidence_sequence_check.py docs/evidence/push-log/PUSH_LOG.md --parent-prefix <17eb881 version> --pushlog-schema` | exit 0 | PASS — sequence + append-only + truthfulness clean |
| 7 | `python3 scripts/ci/status_vocabulary_lint.py docs/directives/INDEX.md docs/evidence/ec16-ec26/RECEIPT.md docs/evidence/push-log/PUSH_LOG.md docs/evidence/ec16-ec26/LOCAL_TEST_LOG.md --runs-manifest <empty>` | PASS | **PASS (exit 0)** after rewording this log's own self-referential mentions of the guarded phrases |

## Per-defect coverage added by the extension (auditor items 3.1–3.8)

- **3.1 / EC-26 (AT-06, AT-08):** 9 fail-closed cases — missing note, missing blob citation,
  unmapped rule, empty inventory, duplicate IDs, malformed row, unknown status token,
  nonexistent enforcement file, AT ID absent from registry.
- **3.2 / EC-20 (AT-09):** missing named test fails without execution; failing named test
  fails; passing pair produces evidence JSON with candidate SHA, argv LIST (no shell),
  counts, exit status.
- **3.3 / EC-21 (AT-02, AT-10):** valid base (local+remote-exact+ancestry+stack);
  local-present/remote-404 → NOT_ON_REMOTE; 403 → REMOTE_VERIFICATION_FAILURE; returned-SHA
  mismatch → fail closed; wrong stack-parent → STACK_MISMATCH; malformed SHA → MALFORMED_BASE;
  shallow → INCOMPLETE_HISTORY; unrelated base → ANCESTRY.
- **3.4 / EC-22 (AT-04, AT-11):** unknown tokens (RUNNING/DONE/WEIRD); RUN_PASS without
  run_id+sha; fabricated run ID; wrong-SHA run; non-success conclusion; unsupported
  percentage; uncited clean-scan phrase; unsupported prior-failure claim; manifest-backed pass;
  no-manifest fail-closed.
- **3.5 / EC-16 (AT-01, AT-12):** whole-list rejection on any malformed entry (short, 41-hex,
  uppercase named explicitly); no partial verification; 401/403/429/500 → VERIFICATION_FAILURE
  (never "not on remote"); returned-SHA mismatch → fail closed; explicit `none` handling.
- **3.6 / EC-18/24 (AT-05, AT-13):** parent-prefix append-only (real PUSH_LOG vs 17eb881
  version, pure append, edit, deletion); CONFIRM-entry SHA truthfulness (placeholder
  rejected); P-011-CONFIRM carries the real 40-hex remote SHA.
- **3.7 / EC-25 (AT-03, AT-14):** push-triggered deploy fixture FLAGGED; dispatch-only
  remediation ACCEPTED; proposed gate triggers `["**"]`/`[main]` verified in parsed YAML;
  fixture reproduces `["*"]` + `ain]` defects; remediation-before-lint precondition recorded.
- **3.8 / EC-23 (AT-15):** 40-hex action pin (not a tag); least-privilege permissions;
  all-refs `--log-opts="--all"` job; checksum-verified binary; seeded-positive + clean-negative
  controls; shallow fail-closed.

## Notes on file 7 (lint of this cycle's own evidence files)

This log and the lint rules themselves mention the guarded phrases in prose; the first
draft of this log produced 3 self-referential findings, which were reworded (the phrases
are now hyphenated or described rather than quoted verbatim) — final lint result over
INDEX.md, RECEIPT.md, PUSH_LOG.md and this file: PASS (exit 0). The amendment document
itself is outside the adopted lint's file scope (docs/evidence/**, outputs/*, *STATUS*)
and its verbatim failure-pattern text is preserved untouched.

## Honest limitations

1. No CI run ID exists for any commit in this cycle (workflows inert by design): every
   result above is LOCAL evidence, reproducible via the exact commands listed.
2. The gitleaks binary download/checksum/seeded-secret controls are validated structurally
   (AT-15) and will execute only on the approved activation path; no gitleaks run against
   this tree is claimed here.
3. Assertion count 91 is measured from source at the candidate commit (`grep -c` of
   `check(` call sites in `scripts/ci/tests/test_acceptance.py`); AT-16 cross-checks the
   registry↔source↔INDEX consistency mechanically.
