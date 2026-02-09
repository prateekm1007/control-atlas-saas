from pathlib import Path
import sys

# --- Canonical project root resolution ---
PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from md_pipeline import run_production_md

input_pdb = PROJECT_ROOT / "outputs/v16_P4/final_model.pdb"
output_dir = PROJECT_ROOT / "entries/097_tier3_md/v16_P4"

print("🚀 Tier-3 MD: v16_P4 (Helical Needle)")
print(f"📁 Project root: {PROJECT_ROOT}")

run_production_md(
    pdb=input_pdb,
    output_dir=output_dir,
    ns=10,
    implicit_solvent="OBC2",
    platform="OpenCL"
)
