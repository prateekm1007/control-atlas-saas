"""
Local cache backend for structure predictions.
"""

import json
import hashlib
from pathlib import Path

DEFAULT_CACHE_DIR = Path(__file__).parent.parent / "cache"

def _key(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()[:16]

def get_cached(sequence: str, cache_dir: Path = DEFAULT_CACHE_DIR):
    cache_dir = Path(cache_dir)
    if not cache_dir.exists():
        return None

    k = _key(sequence)
    pdb = cache_dir / f"{k}.pdb"
    meta = cache_dir / f"{k}.json"

    if pdb.exists() and meta.exists():
        with open(meta) as f:
            m = json.load(f)
        return {
            "status": "SUCCESS",
            "pdb_path": str(pdb),
            "confidence_global": m["confidence_global"],
            "confidence_per_residue": m["confidence_per_residue"],
            "source": "cache",
            "cached": True,
            "error": None
        }
    return None

def save_to_cache(sequence, pdb_string, confidence_global, confidence_per_residue, source, cache_dir=DEFAULT_CACHE_DIR):
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    k = _key(sequence)
    pdb = cache_dir / f"{k}.pdb"
    meta = cache_dir / f"{k}.json"

    pdb.write_text(pdb_string)
    meta.write_text(json.dumps({
        "sequence": sequence,
        "confidence_global": confidence_global,
        "confidence_per_residue": confidence_per_residue,
        "source": source
    }, indent=2))

    return str(pdb)
