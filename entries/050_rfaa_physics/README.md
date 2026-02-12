# Entry 050 — Atomic Physics Upgrade (RFAA)

Integration of RoseTTAFold All-Atom for explicit ligand-protein complex validation.

## Objectives
1.  **Atomic Resolution:** Move beyond residue-level pockets to atomic interactions.
2.  **Bound-State Validation:** Falsify pockets that collapse upon ligand binding.
3.  **Covalent Modeling:** Explicitly model warhead geometry (e.g., KRAS G12C).

## Requirements (GPU Node)
- 400GB+ Database (UniRef30, BFD)
- NVIDIA GPU (A100/H100 recommended)
- Conda env `RFAA`

## Integration Point
This module hooks into `Entry 049` (Screening) as the final "Physics Check" before Chemistry.

## Metrics
- **PAE_interaction:** < 10 Å (Tight binding)
- **pLDDT_interface:** > 85 (Confident structure)
