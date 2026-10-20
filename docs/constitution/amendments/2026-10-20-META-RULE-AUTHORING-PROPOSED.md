# META amendment — Rule-authoring standard (PROPOSED — PENDING OWNER APPROVAL)

- **Status:** PROPOSED. Not in force. Adopting this changes how every future rule is written and classified; it does not by itself reclassify any existing rule except as explicitly proposed for EC-24 in the main amendment document.
- **Source:** auditor message, 2026-10-20, section "What I'd change in how rules are written" (verbatim below).

## Auditor text (verbatim)

> Prose rules failed because nothing checked them. Each rule should name the file, command, or CI job that enforces it. A rule without an enforcement point should be marked advisory, and it shouldn't count toward compliance.

## Proposed constitutional text

**R-AS1 Enforcement point required.** Every rule in `constitution.md` and `ENGINEERING_CONSTITUTION.md` must name at least one enforcement point: a file path, a command, or a CI job ID. A rule may name several.

**R-AS2 Advisory classification.** A rule (or clause) with no implementable enforcement point is marked `ADVISORY` at the point of definition. `ADVISORY` rules do not count toward compliance percentages and cannot gate acceptance of any work item.

**R-AS3 Classification register.** The in-force set, each rule's enforcement points, and each ADVISORY marking live in `docs/directives/INDEX.md` (single index, EC-19). Changes to classification go through the same approval path as the rule itself.

**R-AS4 Retrofit, not rewrite.** On adoption, existing rules are annotated in place with their enforcement points (or marked `ADVISORY`); their wording is otherwise unchanged. No retrospective re-writing of rule text.

## Coder notes

- Applied to this package itself: EC-24's "only authorized credential names can push" clause is the first candidate for `ADVISORY` (no server-side enforcement point on a personal repo — see the EC-24 section of the main amendment document). Recorded as a proposal; requires owner approval.
- This document was written by the coder from the auditor's verbatim paragraph above; the constitutional wording (R-AS1..R-AS4) is the coder's draft for owner review, not auditor-issued text. The verbatim source governs in case of conflict.
