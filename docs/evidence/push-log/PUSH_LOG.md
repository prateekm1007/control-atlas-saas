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

## 2026-10-20T11:45:50Z | P-INTENT-ec20 | push-intent | branch=activation/ec20 base=f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd expected_remote_tip_before_push=ABSENT (new branch) credential=OWNER-PAT-1

- Direction 6 activation branch for EC20: verify-claimed-shas.yml activated verbatim from proposals/ci/
  (EC-24: derived workflow as documented in the file header) on a feature branch; never
  merged until D-5 allows. Expected remote tip before push: ABSENT (new ref name).
