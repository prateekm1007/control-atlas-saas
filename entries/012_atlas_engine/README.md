# Atlas Entry #012: The Engine (Batch Automation)

**"Scaling the Truth"**

## 1. Intent
To demonstrate **Infrastructure Scalability**.
An Atlas is only valuable if it is living code. This entry provides the `run_full_atlas.ps1` engine that regenerates every scientific claim in the repository from source.

## 2. Technical
*   **Artifact:** `run_full_atlas.ps1` (Root Level)
*   **Logic:** Recursive discovery of `*_auto.cxc` scripts.
*   **Execution:** Serial execution via ChimeraX CLI.

## 3. Insight
We do not store images. We store the **instructions to create truth**.
This allows us to update the entire Atlas when ChimeraX upgrades, or when we refine our styling, with zero manual effort.

---
*Provenance: Control-Atlas v0.1*
