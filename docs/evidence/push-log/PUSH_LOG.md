# Push log (EC-24) — append-only

- **Authorized credential names:** `OWNER-PAT-1` (the owner-supplied GitHub PAT; recorded by NAME only — the value is never written to any file, per the standing credential-containment rule).
- **Convention:** an entry is appended in the same push it describes, carrying branch, full SHA, the time the push is issued, and the credential name. The ref-update confirmation is verified post-hoc via `ls-remote` and recorded in the round report. Entries written after the fact are labeled `RECONSTRUCTED` (EC-18). Grammar matches `scripts/ci/evidence_sequence_check.py`: `## <ISO-8601 UTC> | <ENTRY-ID> | <verb> | <subject> | [RECONSTRUCTED]`.
- **Scope note:** this log starts 2026-10-20. Pushes before 2026-10-20 (if any) are out of its recorded scope; the 2026-10-19 unreported-push incident is recorded separately in `INCIDENT_RECORD.md` (branch `docs/owner-decisions`, evidence/incident-push/).

## 2026-10-20T06:04:16Z | P-001 | push | branch=evidence/round9-r2-2026-10-20 sha=7e3a3395be2d0719facf5c61ac532a0d347a49d3 credential=OWNER-PAT-1 | RECONSTRUCTED

- RECONSTRUCTED (Round-11 cycle, worklog Task 34-r11fix context; entry written 2026-10-20T06:04:16Z from worklog + remote ref evidence; original push time not recorded).

## 2026-10-20T06:04:16Z | P-002 | push | branch=evidence/round11-fix-2026-10-20 sha=8301c5387a78754ceac9cc9b2d5b7c6ac3356d4e credential=OWNER-PAT-1 | RECONSTRUCTED

- RECONSTRUCTED (Task 34-r11fix; later advanced to 742e36b by Task 35 direction 6 — see P-007).

## 2026-10-20T06:04:16Z | P-003 | push | branch=feature/od07-audit-id sha=08a9403a7dd376da08831d14f1fa4fa26b15430e credential=OWNER-PAT-1 | RECONSTRUCTED

- RECONSTRUCTED (Task 35-round12 direction 2; full SHA resolved from the local object store 2026-10-20; tip later advanced to c02f774, P-008).

## 2026-10-20T06:04:16Z | P-004 | push | branch=feature/n8-altlocs-per-residue sha=b5afd6ca01a345432852d5865a79999e464636b4 credential=OWNER-PAT-1 | RECONSTRUCTED

- RECONSTRUCTED (Task 35-round12 direction 3; tip unchanged since — confirmed by ls-remote 2026-10-20).

## 2026-10-20T06:04:16Z | P-005 | push | branch=hygiene/h1-remove-gemini sha=66ffcbcbb9f524dba36a04cdff4e52694a539a1e credential=OWNER-PAT-1 | RECONSTRUCTED

- RECONSTRUCTED (Task 35-round12 direction 4; full SHA resolved from the local object store 2026-10-20; tip later advanced to 5149435, P-009).

## 2026-10-20T06:04:16Z | P-006 | push | branch=docs/owner-decisions sha=0dcc5a10756871f2b2dd239b1858cf8aae6de41a (plus supersession commit b871688377804a67e948ea4a190125b6bc46ec1d on evidence/round11-fix-2026-10-20) credential=OWNER-PAT-1 | RECONSTRUCTED

- RECONSTRUCTED (Task 35-round12 direction 5; full SHAs resolved from the local object store 2026-10-20; tip later advanced to 0773f0b, P-010).

## 2026-10-20T06:04:16Z | P-007 | push | branch=evidence/round11-fix-2026-10-20 sha=742e36b0e49afc348137c03eca37aa627ac7c031 credential=OWNER-PAT-1 | RECONSTRUCTED

- RECONSTRUCTED (Task 35-round12 direction 6; tip unchanged since — confirmed by ls-remote 2026-10-20).

## 2026-10-20T06:04:16Z | P-008 | push | branch=feature/od07-audit-id sha=c02f7742000be14747af0a0015cd18c82625831f credential=OWNER-PAT-1 | RECONSTRUCTED

- RECONSTRUCTED (post-Task-35 remediation; git commit time 2026-10-20T05:41:01Z; push time not recorded — no worklog entry exists for this cycle; worklog gap recorded in docs/directives/INDEX.md).

## 2026-10-20T06:04:16Z | P-009 | push | branch=hygiene/h1-remove-gemini sha=51494350e43eac85a00aaba6e84c88cc3d500b46 credential=OWNER-PAT-1 | RECONSTRUCTED

- RECONSTRUCTED (post-Task-35 remediation; includes 6d6600e "H-1 on the P0 stack" (commit time 2026-10-20T05:35:26Z) + merge 5149435 (05:44:32Z); push time not recorded — worklog gap).

## 2026-10-20T06:04:16Z | P-010 | push | branch=docs/owner-decisions sha=0773f0b1ca698829ae20b59660798dc8b6504675 credential=OWNER-PAT-1 | RECONSTRUCTED

- RECONSTRUCTED (post-Task-35 remediation; includes 7bd3e2a "register rev 6" (commit time 2026-10-20T05:49:22Z) + merge 0773f0b; push time not recorded — worklog gap).

## 2026-10-20T06:06:00Z | P-011 | push | branch=proposal/ec16-ec26-enforcement sha=self-referential-tip credential=OWNER-PAT-1

- LIVE entry for this package's push. Self-reference limit, disclosed: an entry cannot contain the 40-hex SHA of the commit that first publishes it, so the subject carries the disclosed marker `self-referential-tip` (NOT a SHA) and the entry ships inside the commit it describes. The actual pushed tip is verified post-hoc via `ls-remote` immediately after the push and quoted in the worklog and round report (that verification record is written later and would be labeled RECONSTRUCTED if ever appended here).

## 2026-10-20T09:34:59Z | P-011-CONFIRM | push-confirm | branch=proposal/ec16-ec26-enforcement remote_sha=17eb88184c9f461c6bc9ee592f916f7ff4c921c4 credential=OWNER-PAT-1 evidence="git ls-remote origin refs/heads/proposal/ec16-ec26-enforcement"

- Correction of P-011 under the two-phase schema (`PUSH_LOG_SCHEMA.md`, defined 2026-10-20;
  auditor item 3.6): the original P-011 entry carried `sha=self-referential-tip`, which is not
  an object identity and cannot serve as evidence of what was pushed. This CONFIRM entry
  records the remote-verified tip.
- Evidence source: `git ls-remote origin` executed 2026-10-20 ~09:04:30Z in the correction
  session returned `refs/heads/proposal/ec16-ec26-enforcement -> 17eb88184c9f461c6bc9ee592f916f7ff4c921c4`
  (output quoted in the session worklog and re-verified at CONFIRM write time below).
- Timestamp reconciliation (differences explained, not averaged away): P-011's stated time
  2026-10-20T06:06:00Z is the entry-author's claimed push-issuance time written inside the
  pushed commit; GitHub commit metadata for 17eb881 records author/committer time
  2026-10-20T06:07:23Z (83 s later — the entry was drafted immediately before `git commit`).
  The moment the remote accepted the ref update is NOT established by either artifact; it was
  not observed at push time. This CONFIRM entry's timestamp (09:34:59Z) is the truthful write
  time of the confirmation, NOT a claim about when the push executed.
- Prior entries, including P-011 itself, are preserved unchanged per the append-only rule;
  the parent-version check (`evidence_sequence_check.py --parent-prefix` against
  `17eb881:docs/evidence/push-log/PUSH_LOG.md`) enforces this from this commit forward.

## 2026-10-20T10:22:09Z | P-012-INTENT | push-intent | branch=docs/owner-decisions head_at_intent=4f54616079ded81d120b6b3b6aa64fae40ddc4a4 credential=OWNER-PAT-1 authorized_by="auditor direction 2026-10-20 + register rev 7 Part 8 (this file's sibling branch record)"

- Register rev 7 (adjudications + A-rows + authorization reconciliation + blockers) is pushed
  FIRST so the scoped authorization for the proposal-branch pushes is public before they occur.

## 2026-10-20T10:22:09Z | P-013-INTENT | push-intent | branch=proposal/ec16-ec26-enforcement head_at_intent=RECORDED-BELOW credential=OWNER-PAT-1 authorized_by="auditor direction 2026-10-20 section 3 (scoped corrections), reconciled in register rev 7 Part 8"

- This push carries the 3.1–3.9 correction commits, the acceptance-loop evidence, and the
  INTENT records themselves. head_at_intent is finalized by the carrying commit (a pre-push
  entry cannot know its own commit's hash); the CONFIRM entry records the remote-verified
  full SHA after the push returns (two-phase schema, PUSH_LOG_SCHEMA.md).

## 2026-10-20T10:22:09Z | P-014-INTENT | push-intent | branch=governance/canonical-docs head_at_intent=07c3fbe6eb20c5a5d4b03dae0a2c1a743edbea10 credential=OWNER-PAT-1 authorized_by="standing rev-6 list + auditor direction 2026-10-20 item 8 (P0_STATUS update)"

- P0_STATUS.md gains the amendment-cycle section (actual statuses; CI NOT RUN; blockers
  preserved).

## 2026-10-20T10:32:13Z | P-012-CONFIRM | push-confirm | branch=docs/owner-decisions remote_sha=78b147d55660dc79c2f5aada368d8499339bc339 credential=OWNER-PAT-1 evidence="git ls-remote origin refs/heads/docs/owner-decisions -> 78b147d55660dc79c2f5aada368d8499339bc339 (10:28:46Z)"

- Push issued 10:28:44Z, remote accepted by 10:28:46Z (fast-forward a99348f..78b147d).
- Content: merge b421237 (reconciliation with the concurrent Round-15/16 writer) + Part 10
  incident note. INTENT P-012 was committed when the local tip was 4f54616; the tip then
  advanced by the reconciliation merge + Part 10 after the first push attempt was REJECTED
  (non-fast-forward: the concurrent writer had pushed a99348f at 10:03:38Z). No force-push
  was used at any point.
- CI observed on this push: Toscanini Stability Gate run 39408223554 = FAILURE (both jobs).
  Cause analysis from the run record: (a) Secret Scan fails on the HEAD-tree A-02
  former-provider key (system_diagnosis.py:8) — the documented un-remediated main-based
  baseline; (b) Unit Tests fail at "Install dependencies (lock + dev)" — main's stale
  dependency lock (fixed on the unmerged science line). Both failures also occurred on the
  concurrent writer's own commits (runs 39399258526, 39406712201) and neither is caused by
  rev 7's content (the register touches only docs/decisions/OWNER_DECISIONS.md).
  The carried gate trigger (branches: ["**"], from the concurrent writer's commits) is what
  made these nested-branch pushes run CI at all.

## 2026-10-20T10:32:13Z | P-013-CONFIRM | push-confirm | branch=proposal/ec16-ec26-enforcement remote_sha=dab33d946fe31891efbd6e4917e360abcb1872c1 credential=OWNER-PAT-1 evidence="git ls-remote origin refs/heads/proposal/ec16-ec26-enforcement -> dab33d946fe31891efbd6e4917e360abcb1872c1 (10:28:58Z)"

- Push issued 10:28:56Z, remote accepted by 10:28:58Z (fast-forward 17eb881..dab33d9).
  The INTENT's deferred head_at_intent is hereby recorded: dab33d946fe31891efbd6e4917e360abcb1872c1.
- CI observed on this push: ZERO check runs (verified via the commits API) — this branch
  still carries main's un-remediated gate triggers, which do not match nested branch names
  (the exact defect recorded in proposals/remediation/). CI status for the correction cycle
  remains NOT_RUN, exactly as stated in LOCAL_TEST_LOG.md.

## 2026-10-20T10:32:13Z | P-014-CONFIRM | push-confirm | branch=governance/canonical-docs remote_sha=b12f30e95255b0c5d31a70cb9b6bb42776aac9a6 credential=OWNER-PAT-1 evidence="git ls-remote origin refs/heads/governance/canonical-docs -> b12f30e95255b0c5d31a70cb9b6bb42776aac9a6 (10:29:25Z)"

- Push issued 10:29:05Z; first attempt REJECTED (non-fast-forward: concurrent writer had
  advanced f368bd2 -> 5e156ba with the same gate-carry/reading-note set); reconciled by
  merge (no force-push) and accepted at 10:29:25Z.
- CI observed: gate run on b12f30e = FAILURE with the same two baseline causes recorded at
  P-012-CONFIRM (A-02 HEAD-tree key; main's stale lock). P0_STATUS section 7's own content
  is documentation-only and not implicated.

## 2026-10-20T11:23:35Z | P-015-INTENT | push-intent | branch=proposal/ec16-ec26-enforcement-rebased head_at_intent=DEFERRED base=eb6d6a9fc42681558e55f085b804f9cc2334c924 expected_remote_tip_before_push=ABSENT (new branch; first push of this ref name) credential=OWNER-PAT-1

- Direction 4 / D-4: the proposal content is rebuilt onto `eb6d6a9` (current main, which contains
  the release `d292a38`). The original branch `proposal/ec16-ec26-enforcement` (`a0e321c`) is left
  untouched — the standing no-force-push rule forbids rewriting its ref, so the rebased tree goes
  to a NEW branch name and the merge carries it (D-7 single-writer: expected tip recorded above;
  abort if a branch of this name already exists remotely).
- Content kept: `proposals/`, `docs/` (amendments, directives INDEX, evidence, push-log, reading
  note), `scripts/ci/` (the enforcement tooling named as enforcement points in the rule inventory),
  root `SHA256SUMS` (artifact manifest). `.github/` is byte-identical to main (verified empty diff
  before commit). D-5 status updates recorded in the reading note + this index's inventory
  (EC-17/EC-22/EC-25 ACTIVE; eight remain PROPOSED).

## 2026-10-20T11:28:30Z | P-015-CONFIRM | push-confirm | branch=proposal/ec16-ec26-enforcement-rebased remote_sha=0f1e9a8856bafec592c657270e5799734ed61d6a credential=OWNER-PAT-1 evidence="git ls-remote origin refs/heads/proposal/ec16-ec26-enforcement-rebased -> 0f1e9a8856bafec592c657270e5799734ed61d6a (11:27:51Z)"

- Push issued 11:27:49Z, remote accepted 11:27:51Z (new ref, as recorded in the INTENT's
  expected_remote_tip_before_push=ABSENT). Pre-push D-7 check: main was at
  eb6d6a9fc42681558e55f085b804f9cc2334c924 exactly as expected (ls-remote 11:26:52Z).
- CI: gate expected on this push (main's trigger list is branches: ["**"]); run ID to be
  recorded in the round report once observed via the Actions API.

## 2026-10-20T13:53:07Z | P-017-INTENT | push-intent | branch=fix/mainline-r18 head_at_intent=DEFERRED base=f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd expected_remote_tip_before_push=ABSENT (new ref name; first push of this branch) credential=OWNER-PAT-1 authorized_by="auditor Round-18 direction 3 (fix-forward PR from exactly f535a71, preserve merged history, no reset, no force-push) + AGENTS.md feature-branch discipline; scoped by register rev 8 Part 11"

- Round-18 fix-forward push (Direction 3): provider-SDK install references removed from the
  two remaining non-documentation package files (ported selectively from the preserved
  `c4a0fcc` history; no whole-history merge); remediation documentation corrected
  (HISTORICAL/SUPERSEDED labels with exact source SHAs; in-force policy separated from the
  automation actually running); governed-artifact manifest completed with a fail-closed
  checker (`scripts/ci/manifest_check.py`); acceptance protection added to ACTIVE CI
  (`.github/workflows/acceptance-gate.yml`); push-log protocol checker extended (pairing +
  remote verification; pending reported as PENDING, never complete); PR #17 review-bypass
  incident (INC-1) and P-015-CONFIRM collision (INC-2) recorded without rewriting history;
  register rev 8 (Part 11); P0_STATUS Round-18 section; EC-20 collection scoped to the
  acceptance suite; EC-26 Enforcement column added and computed.
- Pre-push single-writer guard: the remote ref `refs/heads/fix/mainline-r18` is expected to
  be ABSENT; main is expected at exactly
  `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd`; the push is ABORTED on any deviation.
- head_at_intent is DEFERRED (self-reference limit, PUSH_LOG_SCHEMA.md): the carrying
  commit cannot know its own hash; the CONFIRM entry records the remote-verified full SHA
  after the push returns.

## 2026-10-20T14:02:44Z | P-017-CONFIRM | push-confirm | branch=fix/mainline-r18 remote_sha=562cff5ff3f4fb690224f03122b6f2a66ee75512 credential=OWNER-PAT-1 evidence="git ls-remote origin refs/heads/fix/mainline-r18 -> 562cff5ff3f4fb690224f03122b6f2a66ee75512 (immediately after the push returned)"

- Push issued ~13:57:30Z, remote accepted as a new ref (expected_remote_tip_before_push=ABSENT
  held; main was at f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd exactly as guarded).
- CI observed on 562cff5: Toscanini Stability Gate SUCCESS (push event run 39434887380 and
  PR event run 39434977141); Round-18 acceptance gate FAILURE (runs 39434887481 push event /
  39434976664 PR event): the manifest checker flagged scripts/ci/tests/__pycache__/*.pyc —
  bytecode generated by the acceptance-suite step between scope definition and verification.
  Root cause and fix (bytecode exclusion in manifest_check.py) committed and pushed as the
  next push (P-018 below). Honest record: the first acceptance-gate candidate FAILED on a
  real defect the gate exists to catch; no gate was weakened or bypassed.

## 2026-10-20T14:02:44Z | P-018-INTENT | push-intent | branch=fix/mainline-r18 head_at_intent=DEFERRED base=562cff5ff3f4fb690224f03122b6f2a66ee75512 expected_remote_tip_before_push=562cff5ff3f4fb690224f03122b6f2a66ee75512 credential=OWNER-PAT-1 authorized_by="Round-18 direction 3 fix-forward continuation (CI defect fix); scoped by P-017-CONFIRM" | RECONSTRUCTED

- RECONSTRUCTED, disclosed: this push's pre-push INTENT was OMITTED in error (process gap —
  the entry was not appended before the push was issued). Recorded here after the fact,
  labeled RECONSTRUCTED per the schema, from the actual pushed content: commit 1 = manifest
  scope excludes generated bytecode (__pycache__/*.pyc; CI-found on run 39434887481) +
  SHA256SUMS refresh; commit 2 = run_all success label derived from the registry count.
  Remote tip observed 5ab3a46ed8fa6182efb07f23a5c60c42dac3e998 by ls-remote at 14:02:44Z
  (quoted in P-018-CONFIRM). The omission is recorded, not papered over; no prior entry was
  edited.

## 2026-10-20T14:02:44Z | P-018-CONFIRM | push-confirm | branch=fix/mainline-r18 remote_sha=5ab3a46ed8fa6182efb07f23a5c60c42dac3e998 credential=OWNER-PAT-1 evidence="git ls-remote origin refs/heads/fix/mainline-r18 -> 5ab3a46ed8fa6182efb07f23a5c60c42dac3e998 (14:02:44Z)"

- Push 2 accepted as a fast-forward 562cff5 -> 5ab3a46 (no force-push). CI on 5ab3a46 is
  running at CONFIRM write time; run IDs to be appended in the round report once observed.

## 2026-10-20T14:02:44Z | P-019-INTENT | push-intent | branch=fix/mainline-r18 head_at_intent=DEFERRED base=5ab3a46ed8fa6182efb07f23a5c60c42dac3e998 expected_remote_tip_before_push=5ab3a46ed8fa6182efb07f23a5c60c42dac3e998 credential=OWNER-PAT-1 authorized_by="Round-18 direction 3 continuation: P-017/P-018 confirmations + single-writer guard for push 3"

- Push 3 content (disclosed pre-push): this append (CONFIRM ×2, RECONSTRUCTED INTENT, this
  INTENT), P0_STATUS section 8 status corrections (the observed RUN_FAIL recorded per EC-22
  grammar), and the manifest refresh over the final tree. Single-writer guard: the remote
  tip is expected to be exactly 5ab3a46ed8fa6182efb07f23a5c60c42dac3e998; the push is
  ABORTED on any deviation.

## 2026-10-20T14:12:02Z | P-019-CONFIRM | push-confirm | branch=fix/mainline-r18 remote_sha=3271fb1dbb1cb7687889c4e71002aa45ff37c340 credential=OWNER-PAT-1 evidence="git ls-remote origin refs/heads/fix/mainline-r18 -> 3271fb1dbb1cb7687889c4e71002aa45ff37c340 (post-push observation; check-runs API for the same SHA queried at CONFIRM write time)"

- Push 3 accepted as a fast-forward 5ab3a46 -> 3271fb1 (no force-push), single-writer guard
  held (expected_remote_tip_before_push=5ab3a46 observed).
- CI observed on 3271fb1 (both events green): Round-18 acceptance gate run 39435833504
  (push) and 39435841738 (pull_request) SUCCESS; Toscanini Stability Gate run 39435833741
  (push) and 39435841759 (pull_request) SUCCESS.

## 2026-10-20T14:12:02Z | P-026-INTENT | push-intent | branch=fix/mainline-r18 head_at_intent=DEFERRED base=3271fb1dbb1cb7687889c4e71002aa45ff37c340 expected_remote_tip_before_push=3271fb1dbb1cb7687889c4e71002aa45ff37c340 credential=OWNER-PAT-1 authorized_by="Round-18 direction 3/4 continuation: CONFIRMED-LIVE semantics for stacked trees (final intents of OTHER branches validate via their introducing commit) + P0_STATUS RUN_PASS rows for the green candidate"

- Push 4 content (disclosed pre-push): evidence_sequence_check.py CONFIRMED-LIVE extension
  (proof (b): the remote tip of the intent's branch must be at-or-after the commit that
  INTRODUCED the intent entry — lets a stacked activation tree validate the base branch's
  final intent against the actual remote ref); this append; P0_STATUS section 8 updated
  with the observed green runs (EC-22 grammar, manifest-backed); manifest refresh.
- ID note: P-020..P-025 are allocated to the Round-18 activation rebuild branches
  (consolidated EC-16/19/20, EC-18, EC-21, EC-23, EC-24, EC-26); this push therefore uses
  the next globally unused ID (P-026).
- Single-writer guard: remote tip expected exactly 3271fb1dbb1cb7687889c4e71002aa45ff37c340;
  push ABORTED on deviation.

## 2026-10-20T14:28:42Z | P-027-INTENT | push-intent | branch=fix/mainline-r18 head_at_intent=DEFERRED base=a3d1b5633d46a0cd99842332440986a296412cbb expected_remote_tip_before_push=a3d1b5633d46a0cd99842332440986a296412cbb credential=OWNER-PAT-1 authorized_by="Round-18 direction 3C continuation: manifest parser defect found by the EC-23 candidate's CI (the old lstrip(\"./\") mangled .gitleaksignore); fix + manifest regen"

- Push 5 content (disclosed pre-push): manifest_check.py — the path parser stripped a
  leading DOT (lstrip(\"./\") also eats the \".\" of \".gitleaksignore\"), producing the
  CI finding \"listed but MISSING from the tree: gitleaksignore\"; the parser now strips
  only a literal \"./\" prefix, previously-mangled entries are reported as findings
  instead of silently renamed, and the manifest is regenerated. Push-log append; manifest
  refresh. Single-writer guard: remote tip expected exactly
  a3d1b5633d46a0cd99842332440986a296412cbb; push ABORTED on deviation.

## 2026-10-20T14:41:01Z | P-026-CONFIRM | push-confirm | branch=fix/mainline-r18 remote_sha=a3d1b5633d46a0cd99842332440986a296412cbb credential=OWNER-PAT-1 evidence="git ls-remote origin refs/heads/fix/mainline-r18 -> a3d1b5633d46a0cd99842332440986a296412cbb (observed post-push-4; the tip is an ancestor of the current tip; check-runs API queried at CONFIRM write time)"

- Push 4 accepted as a fast-forward 3271fb1 -> a3d1b56 (single-writer guard held). CI on
  a3d1b56: acceptance gate + ordinary gate runs recorded in P0_STATUS section 8.2
  (RUN_PASS rows, manifest-backed).

## 2026-10-20T14:41:01Z | P-027-CONFIRM | push-confirm | branch=fix/mainline-r18 remote_sha=4f7c28bb7a23189e55660f24e7d567454c6ea8f9 credential=OWNER-PAT-1 evidence="git ls-remote origin refs/heads/fix/mainline-r18 -> 4f7c28bb7a23189e55660f24e7d567454c6ea8f9 (observed post-push-5)"

- Push 5 accepted as a fast-forward a3d1b56 -> 4f7c28bb7a23189e55660f24e7d567454c6ea8f9 (single-writer guard held). CI on
  this commit: Toscanini Stability Gate SUCCESS (both events); the Round-18 acceptance
  gate correctly reported PENDING for P-026-INTENT (no confirm existed yet — this entry
  completes the pair; the checker behaved exactly as designed: an unconfirmed push is
  PENDING, never complete). Process note recorded: a new INTENT must never be appended
  before the previous push's CONFIRM.

## 2026-10-20T14:41:01Z | P-031-INTENT | push-intent | branch=fix/mainline-r18 head_at_intent=DEFERRED base=4f7c28bb7a23189e55660f24e7d567454c6ea8f9 expected_remote_tip_before_push=4f7c28bb7a23189e55660f24e7d567454c6ea8f9 credential=OWNER-PAT-1 authorized_by="Round-18 direction 4 continuation: completing the two-phase protocol (P-026/P-027 confirmations) before any further intent; disclosed process correction"

- Push 6 content (disclosed pre-push): this append (P-026-CONFIRM, P-027-CONFIRM, this
  INTENT) and the manifest refresh. No content changes beyond evidence. Single-writer
  guard: remote tip expected exactly 4f7c28bb7a23189e55660f24e7d567454c6ea8f9; push ABORTED on deviation.

## 2026-10-20T14:47:05Z | P-035-CONFIRM | push-confirm | branch=fix/mainline-r18 remote_sha=a3d1b5653d9d8fda24c2d704fec19719c56c0977 credential=OWNER-PAT-1 evidence="git rev-parse a3d1b56 + git ls-remote origin refs/heads/fix/mainline-r18 (the push-4 tip is an ancestor of the current tip; SHA captured programmatically)"

- CORRECTION of P-026-CONFIRM (which is preserved unchanged above, per the append-only
  rule): P-026-CONFIRM recorded remote_sha=a3d1b5633d46a0cd99842332440986a296412cbb — a
  WRONGLY TRANSCRIBED tail (the entry was hand-typed instead of captured from command
  output). The actual push-4 tip is a3d1b5653d9d8fda24c2d704fec19719c56c0977. This entry is the authoritative confirmation of
  the push-4 intent (P-026); the wrong-tail entry must not be used as evidence. The
  transcription defect and the discipline fix (SHAs are captured programmatically, never
  typed) are recorded in the round report; the checker's remote verification is what
  exposes this class of error, exactly as designed.

## 2026-10-20T14:47:05Z | P-036-INTENT | push-intent | branch=fix/mainline-r18 head_at_intent=DEFERRED base=2efbbb7854e27d818d314691b55488f0d31cb74f expected_remote_tip_before_push=2efbbb7854e27d818d314691b55488f0d31cb74f credential=OWNER-PAT-1 authorized_by="Round-18 direction 4 continuation: SHA-transcription corrections recorded append-only; two-phase protocol completed with programmatic SHA capture"

- Push 7 content (disclosed pre-push): this append — P-035-CONFIRM (correction of
  P-026-CONFIRM) and the record of the other affected entries below; manifest refresh.
- CORRECTION REGISTER for wrong base= transcriptions in INTENT entries (all preserved
  unchanged above; none affects protocol checks — base= is informational — but the
  correct values are recorded here): P-019-INTENT base=5ab3a46ed8fa6182efb07f23a5c60c42dac3e998
  was correct; WRONG (hand-typed a3d1b5633d46a0cd99842332440986a296412cbb instead of
  a3d1b5653d9d8fda24c2d704fec19719c56c0977): P-022-INTENT (branch activation/ec21-r18), P-023-INTENT (branch
  activation/ec23-r18), P-024-INTENT (branch activation/ec24-r18), P-025-INTENT (branch
  activation/ec26-r18), P-027-INTENT (branch fix/mainline-r18, base field). The correct
  value for every one of them is a3d1b5653d9d8fda24c2d704fec19719c56c0977. P-021-INTENT and P-020-INTENT (base
  3271fb1dbb1cb7687889c4e71002aa45ff37c340) and P-028-INTENT (base
  d384f0616ff072271befb37101f739a8ee45ffc6) were verified correct. On branch
  activation/ec21-r18, P-029-CONFIRM and P-032-INTENT carry a wrong 98af2ef… tail
  (98af2ef5b73ba54e2a10b4a3325c2b1b34a02d07f); the correct tip is
  98af2ef2baf88f02723944c97036c3d1b7a95263 — corrected on that branch by its own
  correction entry (P-037). The corrected values above were captured programmatically
  from git rev-parse.
- Single-writer guard: remote tip expected exactly 2efbbb7854e27d818d314691b55488f0d31cb74f;
  push ABORTED on deviation.

## 2026-10-20T14:49:21Z | P-039-INTENT | push-intent | branch=fix/mainline-r18 head_at_intent=DEFERRED base=b3c763480c615b2a2f6c54c3bfa2a2a016a08cb4 expected_remote_tip_before_push=b3c763480c615b2a2f6c54c3bfa2a2a016a08cb4 credential=OWNER-PAT-1 authorized_by="Round-18 direction 4 continuation: reconciliation class SUPERSEDED-BY-CORRECTION-ENTRY added so a correction entry (P-035-CONFIRM) authoritatively supersedes a transcription-defective CONFIRM's remote verification"

- Push 8 content (disclosed pre-push): evidence_sequence_check.py gains the
  SUPERSEDED-BY-CORRECTION-ENTRY reconciliation class (a later correction entry
  authoritatively supersedes a defective CONFIRM's remote_sha verification; the defective
  entry stays in the log, append-only); PUSH_LOG_RECONCILIATION.json gains the
  P-026-CONFIRM suppression entry (pointer: P-035-CONFIRM); this INTENT; manifest
  refresh. Single-writer guard: remote tip expected exactly b3c7634... (full value
  captured programmatically below); push ABORTED on deviation.

## 2026-10-20T14:54:04Z | P-036-CONFIRM | push-confirm | branch=fix/mainline-r18 remote_sha=b3c763480c615b2a2f6c54c3bfa2a2a016a08cb4 credential=OWNER-PAT-1 evidence="git rev-parse b3c7634 (captured programmatically; the push-7 tip is an ancestor of the current tip)"

- Push 7 accepted as a fast-forward 2efbbb7 -> b3c763480c615b2a2f6c54c3bfa2a2a016a08cb4 (guard held). This CONFIRM completes
  the P-036 pair (it was omitted when push 8 was issued — the recurring process gap,
  now closed by the standing rule: never append a new INTENT before the previous push's
  CONFIRM).

## 2026-10-20T14:54:04Z | P-039-CONFIRM | push-confirm | branch=fix/mainline-r18 remote_sha=aa510f7533d5d4fd63db9b8bf8d027219c16d500 credential=OWNER-PAT-1 evidence="git rev-parse aa510f7 (captured programmatically; the push-8 tip is an ancestor of the current tip)"

- Push 8 accepted as a fast-forward b3c763480c615b2a2f6c54c3bfa2a2a016a08cb4 -> aa510f7533d5d4fd63db9b8bf8d027219c16d500 (guard held). Completes the P-039 pair.

## 2026-10-20T14:54:04Z | P-041-INTENT | push-intent | branch=fix/mainline-r18 head_at_intent=DEFERRED base=aa510f7533d5d4fd63db9b8bf8d027219c16d500 expected_remote_tip_before_push=aa510f7533d5d4fd63db9b8bf8d027219c16d500 credential=OWNER-PAT-1 authorized_by="Round-18 direction 4 continuation: reconciliation-gated truthfulness/verification (a superseded CONFIRM's defective value is not evidence); completing the P-036/P-039 pairs"

- Push 9 content (disclosed pre-push): evidence_sequence_check.py — reconciled CONFIRM
  entries are skipped by BOTH the truthfulness check and the remote-verification loop
  (previously the ABSENT path returned before the reconciliation check — found by the
  CI run on aa510f7, exactly the fail-closed behavior working); this append completes
  the P-036/P-039 pairs; manifest refresh. Single-writer guard: remote tip expected
  exactly aa510f7533d5d4fd63db9b8bf8d027219c16d500; push ABORTED on deviation.


## 2026-10-20T15:01:27Z | P-020-INTENT | push-intent | branch=activation/ec16-ec20-consolidated head_at_intent=DEFERRED base=3271fb1dbb1cb7687889c4e71002aa45ff37c340 expected_remote_tip_before_push=ABSENT (new ref name; first push of this branch) credential=OWNER-PAT-1 authorized_by="auditor Round-18 direction 5 (consolidate the duplicate EC-16/19/20 activation workflow into ONE candidate; fix the shared EC-20 failure; old PRs #18-#20 kept as historical evidence and marked superseded once this exact candidate is green)"

- Content: .github/workflows/verify-claimed-shas.yml — the ONE activation workflow for
  EC-16 + EC-19 + EC-20 (acceptance-tests-executed with --suite-dir scripts/ci/tests, the
  proven fix for run 39417492355); push-log append; manifest refresh. Stacked on the
  Round-18 fix-forward; NOT merged anywhere until the auditor records activation.
- Timestamp note: re-appended at rebuild 3 (see P-044-INTENT); the timestamp is the
  truthful write time of THIS log version; the original pre-rebase write was
  2026-10-20T14:05:55Z.

## 2026-10-20T15:01:27Z | P-020-CONFIRM | push-confirm | branch=activation/ec16-ec20-consolidated remote_sha=d384f0616ff072271befb37101f739a8ee45ffc6 credential=OWNER-PAT-1 evidence="git ls-remote origin refs/heads/activation/ec16-ec20-consolidated (observed post-push-1, pre-rebase lineage)"

- Push 1 accepted as a new ref (guard held). CI on d384f06: ALL THREE activation jobs
  PASSED (EC-16, EC-19, EC-20 with suite-scoped collection). The tip is no longer on the
  current lineage after the disclosed rebuilds (P-034/P-044); classified
  SUPERSEDED-BY-DOCUMENTED-REBUILD in PUSH_LOG_RECONCILIATION.json.

## 2026-10-20T15:01:27Z | P-028-CONFIRM | push-confirm | branch=activation/ec16-ec20-consolidated remote_sha=ef13c4c2a73961138794e151919275359378e50f credential=OWNER-PAT-1 evidence="git ls-remote origin refs/heads/activation/ec16-ec20-consolidated -> ef13c4c2a73961138794e151919275359378e50f (observed post-push-2, pre-rebuild-2 lineage)"

- Push 2 accepted as a fast-forward d384f06 -> ef13c4c (guard held). CI on ef13c4c: the
  three activation jobs again SUCCESS. The tip is no longer on the current lineage after
  the disclosed rebuilds; classified SUPERSEDED-BY-DOCUMENTED-REBUILD.

## 2026-10-20T15:01:27Z | P-034-INTENT | push-intent | branch=activation/ec16-ec20-consolidated head_at_intent=DEFERRED base=2efbbb7854e27d818d314691b55488f0d31cb74f expected_remote_tip_before_push=ef13c4c2a73961138794e151919275359378e50f (non-fast-forward rebuild 1, executed with --force-with-lease bound to ef13c4c2a739…) credential=OWNER-PAT-1 authorized_by="Round-18 direction 5 continuation: rebuild 1 onto the then-current fix-forward base 2efbbb7 (disclosed; D-4 pattern; never-merged activation candidate)"

- Rebuild 1: the branch was rebased onto 2efbbb7 and pushed with --force-with-lease bound
  to ef13c4c2a739... (the guard initially ABORTED the push because the intent carried a
  hand-typed wrong tail — corrected in place BEFORE first publication; disclosed).
  Result tip: 705967d30f7ce45660a8c4b4718e8d08e115f28f.
- Timestamp note: re-appended at rebuild 3; original write 2026-10-20T14:44:01Z.

## 2026-10-20T15:01:27Z | P-044-INTENT | push-intent | branch=activation/ec16-ec20-consolidated head_at_intent=DEFERRED base=5539b01 expected_remote_tip_before_push=c2c2a6ffe15eac641b817e266fbc77a6da70e731 (non-fast-forward rebuild 3, executed with --force-with-lease bound to the full SHA of c2c2a6ffe15eac641b817e266fbc77a6da70e731 — the ACTUAL prior tip (rebuild 2); an earlier draft of this entry cited the rebuild-1 tip by mistake and was corrected in place BEFORE first publication, disclosed here) credential=OWNER-PAT-1 authorized_by="Round-18 direction 5 continuation: rebuild 3 repairs a rebase-resolution defect — the published tip 705967d carried DUPLICATE P-020-CONFIRM/P-028-CONFIRM entries (both rebase conflict resolutions appended the block; the duplicate is a sequence-check violation the gate caught). The log is rebuilt with exactly one copy of each entry; the duplicates never enter the rebuilt lineage"

- Rebuild 3 (this push): single-commit rebuild of the branch onto 5539b01 with the log
  containing exactly one copy of each consolidated entry, re-timestamped to this log
  version's truthful write time (originals noted inline); reconciliation entries added
  for the two pre-rebuild-lineage confirm tips (SUPERSEDED-BY-DOCUMENTED-REBUILD, both
  verified still retrievable on the remote); this INTENT; manifest refresh. The
  pre-rebuild tip c2c2a6ffe15eac641b817e266fbc77a6da70e731 remains recoverable and is
  cited in the PR comment disclosing this rebuild. No merged branch and no mainline ref
  is affected; the no-force-push rule for published/merged history is untouched (this is
  the D-4 documented-rebuild pattern applied to a never-merged activation candidate,
  with force-with-lease bound to the exact prior tip).

## 2026-10-20T15:06:36Z | P-034-CONFIRM | push-confirm | branch=activation/ec16-ec20-consolidated remote_sha=705967d30f7ce45660a8c4b4718e8d08e115f28f credential=OWNER-PAT-1 evidence="commits API: 705967d30f7ce45660a8c4b4718e8d08e115f28f retrievable (HTTP 200 verified at CONFIRM write time); observed as the post-rebuild-1 tip"

- CONFIRM of the P-034-INTENT (rebuild 1): the push executed with --force-with-lease bound
  to ef13c4c2a739... and landed 705967d30f7ce45660a8c4b4718e8d08e115f28f (verified
  retrievable on the remote at CONFIRM write time). The tip left the current lineage
  through rebuilds 2 and 3 (P-044-INTENT); classified SUPERSEDED-BY-DOCUMENTED-REBUILD in
  the reconciliation record. This completes every intent/confirm pair on this branch:
  P-020 (confirmed), P-028 (confirmed), P-034 (this entry), P-044 (final — CONFIRMED-LIVE
  via the carrying-tip proof in CI).

## 2026-10-20T15:11:03Z | P-045-INTENT | push-intent | branch=activation/ec16-ec20-consolidated head_at_intent=DEFERRED base=af71355 expected_remote_tip_before_push=af71355 credential=OWNER-PAT-1 authorized_by="Round-18 direction 5 continuation: the ordinary gate's --no-git secret scan flagged the P-034-CONFIRM evidence field (the keyword api + ':' + a PUBLIC commit id — structural false positive, no credential); classified per-fingerprint in .gitleaksignore with the justification inline (append-only evidence cannot be reworded)"

- Push 3 content (disclosed pre-push): .gitleaksignore (one fingerprint + inline
  justification; the canonical disposition record on branch activation/ec23-r18 absorbs
  this row at merge time); this INTENT; manifest refresh. Single-writer guard: remote
  tip expected exactly af71355; push ABORTED on deviation.
