# INTERNAL CROSSWALK — retired label spaces → OD numbers

> **INTERNAL.** This is the only file in the current delivery that may
> carry retired D-xx labels (U5). No user-facing document uses them.
> The mappings below state only what the auditor's checklists or
> directions actually stated; where no mapping was stated, the internal
> label is retained with "not stated" and no OD number is guessed.

## 1. r1 register internal labels (D-1 … D-12) → OD numbers

| Internal label | Subject (r1 register) | OD mapping |
|---|---|---|
| D-1 | Verdict precedence | not stated by auditor |
| D-2 | LAW-120 / V4 containment and threshold treatment | OD-14 (auditor checklist) |
| D-3 | LAW-130/155/200 scope; LAW-190/210/220/230 status | OD-03/OD-04 for the pLDDT pair (checklist); LAW-190…230 → OD-05 (consolidated table); earlier drafts left this not stated |
| D-4 | Cohort and split prespecification | OD-10 (checklist); elaborated by OD-18/OD-19 |
| D-5 | P0DTC2 / P62988 identity records | OD-11 and OD-12 (checklist) |
| D-6 | AlphaFold version pinning | OD-13 (checklist) |
| D-7 | Certificate and audit record policy | OD-06 (certificate) and OD-07 (audit-ID) per consolidated table |
| D-8 | Audit IDs / tenant-scoped history | OD-07 and OD-08 per consolidated table |
| D-9 | Runner provenance | not stated as an OD; delivered this round as a U9 unblocked patch |
| D-10 | CSS semantics and residue ordering | not stated by auditor |
| D-11 | v23 cohort metadata quarantine | OD-17 (auditor T2 / consolidated table) |
| D-12 | Credentials, cloud secrets, branch protection | owner actions A-01…A-04; OD mapping not stated |

## 2. Auditor Round-9-continuation label collisions (historical)

The auditor's Round-9-continuation message used "D-4" (cohort
eligibility) and "D-11" (stratum population) as ad-hoc labels. Those two
uses are **OD-18** and **OD-19** respectively. They are NOT the r1
register's D-4 (cohort and split prespecification → OD-10) or D-11
(metadata quarantine → OD-17). Both label spaces are preserved here to
prevent conflation; neither is used going forward.

## 3. OD-16 collision history

A prior coder-side register draft briefly used "OD-16" for metadata
integrity, colliding with the auditor's OD-16 (merge and deploy
sequence). Per the auditor's renumbering: **OD-16 = merge and deploy
sequence (auditor's meaning, unchanged)**; metadata integrity is
**OD-17**. The collision existed only in superseded drafts (sealed
inside the r1/r2 evidence archives as historical records); no current
document uses OD-16 with any other meaning.

## 4. Superseded register drafts (historical, sealed)

- r1 evidence archive `04-decisions/OWNER_DECISIONS.md` — rev 1, D-1…D-12
  internal labeling; predates the OD numbering.
- r2 evidence archive `04-decisions/OWNER_DECISIONS.md` — rev 3, OD-01…21
  with "not stated" placeholders for OD-05…OD-09; superseded by rev 4
  (the consolidated table filled those entries verbatim).
- Both archives are byte-identical sealed records (SHA-256 values in
  `SHA256SUMS` and `DELIVERY_README.md`) and are not amended.
