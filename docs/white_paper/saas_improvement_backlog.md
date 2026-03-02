# Toscanini SaaS Improvement Backlog
# Source: White paper data analysis, 2026-03-02
# Canon hash at time of analysis: 6a9cd4b4349b81de

---

## TICKET-WP-01: Dual Score Output
Priority: HIGH
Source: Section 2.2 — advisory_experimental ceiling
File: backend/main.py, ToscaniniResponse model

Add to verdict output:
  det_score_raw:        current score (existing)
  det_score_normalized: score over 10-law source-invariant subset

Laws excluded from normalized subset: LAW-100, LAW-160
No canon change. No threshold change. Response schema additive only.
Backward compatible — existing consumers see same det_score_raw.

---

## TICKET-WP-02: Resolution-Conditional Rotamer Tier
Priority: MEDIUM
Source: MYO_1MBN VETO (+1.60%), HBA_1HHO VETO (+0.35%) at 2.0-2.1 A

Add to law result output for LAW-150:
  resolution_context:   stated resolution from PDB header
  tolerance_band:       expected rotamer range at this resolution
  veto_type:            "geometric" | "resolution_expected"

Canon threshold unchanged at 20%.
Verdict unchanged.
Interpretation layer only — adds context to existing VETO.

---

## TICKET-WP-03: INDETERMINATE Subcategory
Priority: MEDIUM
Source: p53_AF (59.8% cov, disordered), INS_AF (12.7% cov, prepro sequence)

Add to verdict output when binary == INDETERMINATE:
  indeterminate_reason: "low_confidence" | "disordered_region" | "sequence_mismatch"

Classification logic:
  coverage < 30%:                        sequence_mismatch (check residue count vs canonical)
  30% <= coverage < 70%, predicted:      low_confidence
  30% <= coverage < 70%, known disorder: disordered_region (requires disorder DB lookup)

Default to low_confidence if disorder DB unavailable.

---

## TICKET-WP-04: LAW-200 Domain-Normalized Packing
Priority: LOW
Source: EGFR_AF SD 105.277 driven by 1210-residue structure

Current: bounding box / total atoms (global)
Proposed: median packing per chain segment (rolling 100-residue window)
Report: median_packing alongside existing observed value

No threshold change. Adds stability for large proteins.

---

## TICKET-WP-05: Statistical Summary Endpoint
Priority: LOW
Source: White paper descriptive statistics section

New endpoint: GET /audit/population_stats
  Input:  list of audit_ids
  Output: per-law mean, SD, min, max across submitted audits
          grouped by source_type (experimental / predicted)

Enables researchers to compute their own distributions
without exporting raw audit records.

---

## TICKET-WP-06: Functional Law Tier (v25.0 planning)
Priority: FUTURE — requires governance review
Source: Section 6.6 — KRAS G12D passes geometric audit

Scope: Tier 2 functional laws require:
  - Sequence-to-structure mapping
  - Reference structure database
  - Active site annotation source (e.g. UniProt, CSA)

Not buildable without external data dependencies.
Governance review required before any canon addition.
Design doc to be written separately.

---

## TICKET-WP-07: Violin Plot Endpoint
Priority: LOW
Source: Visualization gap identified in white paper review

Add to dashboard: per-law violin plot for population audits
Input: list of audit_ids, grouped by source_type
Output: matplotlib violin figure served as PNG

Complements existing PDF certificate.
No audit logic change.

