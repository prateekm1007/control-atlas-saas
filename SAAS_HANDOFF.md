# 🛡️ Sovereign Sieve (Falsify.bio) Handoff: v9.2-Industrial

## 1. System Identity & Mission
*   **Project:** Sovereign Sieve (Toscanini Platform)
*   **Status:** Beta v9.2.0 (Operational / Productized)
*   **Mission:** To provide a **Falsification-as-a-Service (FaaS)** layer for the generative AI protein design industry. We don't predict success; we **veto impossibility** using industrial physics (Amber14/OBC2).
*   **Axiom:** "The filter is more valuable than the generator."

## 2. Directory Architecture
The platform is built as a multi-container Docker ecosystem:
*   `backend/app/core/mdi_codex.json`: The 64-law "Constitution."
*   `backend/app/services/physics_engine.py`: The logic kernel enforcing the MDI.
*   `backend/dashboard.py`: The "Forensic Station" UI (Streamlit/3Dmol.js).
*   `tools/`: Native physics kernels from the research cycles (e.g., `native_audit.py`).

## 3. Enforced Law Summary (corrected to the pinned registry, WP-8)
The pinned engine governs **15 unique forensic laws — 14 enforced and LAW-155 retired** (see
`governance/law_registry.yaml`; runtime source
`backend/tos/governance/station_sop.py::LAW_CANON`, canon hash
`6a9cd4b4349b81de`):
*   **GEOMETRIC:** LAW-100 Bond Integrity, LAW-105 Reliability Coverage, LAW-110 Backbone Gap, LAW-120 Bond Angle, LAW-125 Ramachandran, LAW-130 Clashscore, LAW-135 Omega Planarity.
*   **CHEMICAL:** LAW-145 Chirality, LAW-150 Rotamer Audit, LAW-170 Residue Identity, LAW-182 Hydrophobic Burial, LAW-195 Disulfide Geometry.
*   **HEURISTIC:** LAW-200 Packing Quality (heuristic failures do not by themselves trigger VETO).
*   **RETIRED:** LAW-155 Voxel Occupancy (ID retained for traceability; measurement is not enforced), plus LAW-190/210/220/230.
*   LAW-160 Chain Integrity (deterministic).

> **WP-8 label corrections.** This document previously cited laws that do
> not exist in the pinned set and are hereby corrected — LAW-138
> "Deformation", LAW-151 "Rho", LAW-152 "Ghost Decoupling", LAW-155.1
> "Retention", LAW-161 "BSA", LAW-171 "Grease Limit", LAW-180 "R/W Tax",
> and LAW-185 "Self-Healing / 1ns Kinetic Burst" are not enforced laws.
> LAW-155 is **Voxel Occupancy** (retired from enforced semantics), not "Steric Clash".
> LAW-170 is **Residue Identity**, not "pI Balance".
> LAW-190/210/220/230 are **retired** (never measured).
> IDs are never renumbered or reused.

## 4. Current Operational Baseline
1.  **Forensic Trace:** Every audit generates a ledger showing which laws were Passed, Vetoed, or Blocked by causality.
2.  **Visual Proof:** 3Dmol.js renders clashing atom pairs with yellow dashed measurement lines.
3.  **Analytics:** Interaction Matrix (Heatmap) tracks co-occurring AI failure modes.

## 🚀 Priority Directives for Next Session
1.  **Monetization:** Implement Stripe subscription for "Industrial Tier" access.
2.  **Tier 3 Logic:** Build the "1ns Endurance Burst" using OpenMM (no pinned law ID attaches to it — see WP-8 note above).
3.  **Market Deployment:** Transition from Localhost to Production (Railway/AWS).
