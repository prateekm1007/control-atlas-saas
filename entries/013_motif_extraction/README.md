# Atlas Entry #013: Motif Extraction (M-001)

**"The Component Library"**

## 1. Intent
To transition from analyzing proteins to **harvesting control parts**.
We isolate the conserved "Switch" geometry from the KRAS scaffold to create a standalone, searchable geometric object.

## 2. Technical
*   **Source:** `6OIM` (Chain A)
*   **Definition:** Switch I (30-38) + Switch II (60-76) [Vetter 2001].
*   **Operation:** Inverse deletion (Scaffold removal).
*   **Asset Generated:** `library/motifs/M001_Switch_GTPase.pdb`

## 3. Insight
By separating the control logic from the protein identity, we create a **search query**. We can now scan the PDB for *other* proteins that contain this exact 3D arrangement of atoms, identifying hidden targets by mechanism rather than homology.

---
*Provenance: Control-Atlas v0.2 (Phase 2)*
