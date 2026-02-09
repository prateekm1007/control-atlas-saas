# Asset Dossier: TwinRod-v2

## Identity
- **Sequence:** KAWAKKAAAKAEAAKAEAAKYWPTG-PAEAAKP-KAWAKKAAAKAEAAKAEAAKYWPTG
- **Length:** 57 aa (symmetric bipod)
- **Architecture:** Two CHAMP-005 units linked by PAEAAKP

## Verified Metrics
**Source:** structures/twinrod_v2/audit.json (SHA256: 4530118...)

| Metric | Value |
|--------|-------|
| Contact Density (ρ) | 336 |
| Minimum Clearance | 2.32 Å |
| Steric Clashes | 1 (at 2.5Å threshold) |
| Target Heavy Atoms | 1924 |
| Binder Heavy Atoms | 417 |

**Chai-1 Confidence:**
- Aggregate Score: 0.263
- PTM: 0.669
- iPTM: 0.162

## Geometric Audit Status
**❌ CLASH_VETO** (Tier 1 Failure)

One atom pair at 2.32Å violates the 2.5Å hard clash threshold per LAW-148.

## Disposition
Moved to Negative Knowledge Graph. The bipodal architecture shows promise 
(ρ=336 is high contact density) but requires clash resolution before 
progressing to Tier 2 validation.

## Artifact Provenance
- **File:** structures/twinrod_v2/structure.cif
- **SHA256:** 4530118486fd8a746a6cf8cc32a22f8952fbc88c89c5dd3b15bdb22dc4938271
- **Generated:** 2026-01-14 (Kaggle)
- **Extracted:** 2026-01-16 (Local)
- **Audited:** 2026-01-16 (standard_physics_audit v1.0)
