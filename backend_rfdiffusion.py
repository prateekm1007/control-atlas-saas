from pathlib import Path
import subprocess

def generate_structure(candidate, output_pdb: Path):
    output_pdb = Path(output_pdb)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)

    # External RFdiffusion call (example minimal invocation)
    subprocess.run(
        [
            "python3",
            "rfdiffusion/run_inference.py",
            f"inference.output_prefix={output_pdb.with_suffix('')}",
            "inference.num_designs=1"
        ],
        check=True
    )

    # Normalize RFdiffusion output to expected PDB name
    candidates = list(output_pdb.parent.glob(output_pdb.stem + "*_0.pdb"))
    if not candidates:
        raise RuntimeError("RFdiffusion did not produce a PDB")

    candidates[0].rename(output_pdb)
