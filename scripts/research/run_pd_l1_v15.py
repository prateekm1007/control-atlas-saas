from pathlib import Path
from generator_guard import guard
import sys

def generate(candidate, out_pdb: Path):
    if not guard(candidate):
        raise RuntimeError("⛔ Illegal candidate blocked by guard")

    print("🧬 Generating with Chai-1:")
    print(candidate)

    # === REAL INFERENCE MUST HAPPEN HERE ===
    try:
        from chai_model import infer_structure  # <-- replace with real backend
    except ImportError:
        raise RuntimeError(
            "⛔ Chai inference backend not wired. "
            "Replace infer_structure import with your actual Chai call."
        )

    structure = infer_structure(
        buffer=candidate["buffer"],
        anchor=candidate["anchor"],
        scaffold=candidate["scaffold"]
    )

    out_pdb = Path(out_pdb)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    structure.write_pdb(out_pdb)

    print(f"✅ PDB written to {out_pdb}")

if __name__ == "__main__":
    if "--out" not in sys.argv:
        raise RuntimeError("⛔ --out <pdb> is required")

    candidate = {
        "buffer": sys.argv[sys.argv.index("--buffer") + 1],
        "anchor": sys.argv[sys.argv.index("--anchor") + 1],
        "scaffold": sys.argv[sys.argv.index("--scaffold") + 1],
    }

    out_pdb = sys.argv[sys.argv.index("--out") + 1]
    generate(candidate, out_pdb)
