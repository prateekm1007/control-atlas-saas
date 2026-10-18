from pathlib import Path
from typing import Dict

class BackendNotReady(Exception):
    pass

def generate_structure(
    candidate: Dict,
    output_pdb: Path
):
    """
    Must:
    - generate atomic coordinates
    - write a valid PDB with ATOM records to output_pdb
    - raise BackendNotReady if backend is not wired
    """
    raise BackendNotReady("Backend not wired")
