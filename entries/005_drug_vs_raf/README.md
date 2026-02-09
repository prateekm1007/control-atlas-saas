# Atlas Entry #005: The Blockade (AMG 510 vs RAF)

**"Two Objects, One Space"**

## 1. Biological Intent
This entry visualizes the mechanism of action (MoA) for KRAS G12C inhibitors.
It proves **mutual exclusivity**: The conformation locked by AMG 510 is sterically incompatible with RAF-RBD binding.

* **Primary State:** KRAS G12C + AMG 510 (PDB: `6OIM`)
* **Competitor:** RAF-RBD (from PDB: `4G0N`)
* **Key Insight:** Drug binding enforces a topology that physically rejects effector recruitment.

## 2. Technical Implementation
* **Superposition:** Effector-bound Ras aligned to drug-bound Ras
* **Filtering:** Ras from 4G0N hidden; RAF retained as phantom
* **Visual Language:**
  * Gray — Drug-bound KRAS
  * Red — AMG 510
  * Tan — RAF-RBD
  * Orange — Steric clash zone

## 3. Automation
* **Script:** `entry_005_overlay_competition_auto.cxc`
* **Logic:** Load → Match → Filter → Clash Detection → Render

---
*Provenance: Control-Atlas v0.1*
