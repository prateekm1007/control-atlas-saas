# Directive index — single source of series status (EC-19)

- **Started:** 2026-10-20 by the EC-16..EC-26 proposal package (branch `proposal/ec16-ec26-enforcement`, base `main` @ `2195a048e87dcec4aa63e8dc7981b3bb6d237ca0`).
- **Status of this file:** rows marked RECONSTRUCTED were backfilled at index creation from the worklog and remote git metadata, not written contemporaneously (EC-18 discipline). Going forward, rows are appended at series open/close time.
- This index is the machine-readable source for: the ACTIVE series (R-DL1), in-force rule IDs and their enforcement points (EC-26), and superseded directives.

## Series register

| Series ID | Source | Opened (UTC) | Status | Items closed | Notes |
|---|---|---|---|---|---|
| SERIES-2026-10-20-R12 | Round-12 auditor directions 1–6 | 2026-10-20 | COMPLETE (RECONSTRUCTED) | 6/6 | Worklog Task 35-round12: reading note; OD-07 (`feature/od07-audit-id`); N8 (`feature/n8-altlocs-per-residue`); H-1 (`hygiene/h1-remove-gemini`); register (`docs/owner-decisions`); evidence-README retraction (`evidence/round11-fix-2026-10-20`). Tips later advanced by SERIES-2026-10-20-R13R14. |
| SERIES-2026-10-20-R13R14 | Round-13/14 audit remediation (commits only; directive text not present in this environment) | 2026-10-20 (commit times 05:35–05:49Z) | COMPLETE on remote (RECONSTRUCTED) | OD-07 canonical triple (`c02f774`); H-1 re-based onto P0 stack (`5149435`); register rev 6 (`0773f0b`) | **Worklog gap:** no worklog entry exists for this cycle; recorded 2026-10-20 from remote evidence (ls-remote + fetched commits). CI: `c02f774` success (runs 118475032870/118475032461), `5149435` success (runs 118475776494/118475776358), `0773f0b` no runs (main trigger defect). |
| SERIES-2026-10-20-EC | Auditor message, 2026-10-20: EC-16…EC-26 amendment proposals + two meta-amendments + adoption order | 2026-10-20 | **ACTIVE** — packaging; PROPOSALS PENDING OWNER APPROVAL | — | Proposal bundle = this branch (`proposal/ec16-ec26-enforcement`). Nothing in force until owner approval. Adoption order (auditor): Batch 1 {EC-17, EC-25}; Batch 2 {EC-22, EC-23}; Batch 3 {rest}. |

## In-force rule set (EC-26)

Until the owner approves the EC-16..EC-26 amendments, the in-force set is unchanged: the EC sections enforced by the pinned constitutions (`ENGINEERING_CONSTITUTION.md` blob `b6bb35e8…`, `constitution.md` blob `211ad2c5…`) plus OBS-153 / CRIT-154 as explicitly out-of-scope — as mapped in the Round-12 reading note (worklog Task 35, first commit of `feature/od07-audit-id`).

Upon approval, each adopted EC-16…EC-26 rule is appended here with its enforcement point (file/command/CI job) per the rule-authoring standard.

## Superseded directives

None recorded in this index. Historical supersessions prior to 2026-10-20 (register revs 1/3/4 rows, D-cycle, etc.) are recorded in `docs/decisions/OWNER_DECISIONS.md` (register rev 6, branch `docs/owner-decisions`) and are not duplicated here.
