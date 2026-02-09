from pathlib import Path
from backend_rfdiffusion import generate_structure

def run_inference(
    buffer: str,
    anchor: str,
    scaffold: str,
    hinge: str,
    output_pdb: Path
):
    candidate = {
        "buffer": buffer,
        "anchor": anchor,
        "scaffold": scaffold,
        "hinge": hinge
    }
    generate_structure(candidate, output_pdb)
