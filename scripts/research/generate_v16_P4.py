from pathlib import Path
from generator_guard import guard

# --- Candidate definition (LAW CLEAN) ---
candidate = {
    "buffer": "GG",
    "anchor": "YWPG",
    "scaffold": "helical_4turn"
}

assert guard(candidate), "⛔ Candidate blocked by generator guard"

outdir = Path("outputs/v16_P4")
outdir.mkdir(parents=True, exist_ok=True)
pdb_path = outdir / "final_model.pdb"

print("🧬 Running inference for v16_P4")

# === REQUIRED: REAL INFERENCE CALL ===
# You MUST replace the function below with your actual generator.
# It must WRITE a valid PDB with ATOM records to pdb_path.

from chai_inference import run_inference  # adjust import to your repo

run_inference(
    buffer=candidate["buffer"],
    anchor=candidate["anchor"],
    scaffold=candidate["scaffold"],
    hinge="PPPP",
    output_pdb=pdb_path
)

print(f"✅ v16_P4 written to {pdb_path}")
