#!/usr/bin/env python3
"""
Entry 028 — Structure Prediction Ingress
"""

import argparse
import json
from pathlib import Path

try:
    from backends.esmfold import predict_structure as esmfold_predict
    from backends.cache import get_cached, save_to_cache
except ImportError:
    from .backends.esmfold import predict_structure as esmfold_predict
    from .backends.cache import get_cached, save_to_cache


class StructureProvider:
    def __init__(self, backend: str = "esmfold", cache_dir: Path = None):
        self.backend = backend
        self.cache_dir = cache_dir or (Path(__file__).parent / "cache")

    def predict(self, sequence: str) -> dict:
        sequence = sequence.upper().strip()

        cached = get_cached(sequence, self.cache_dir)
        if cached:
            print(f"[*] Cache hit")
            return cached

        print(f"[*] Calling {self.backend} API...")

        if self.backend == "esmfold":
            result = esmfold_predict(sequence)
        else:
            return {
                "status": "ERROR",
                "pdb_path": None,
                "confidence_global": 0.0,
                "confidence_per_residue": [],
                "source": self.backend,
                "cached": False,
                "error": f"Unknown backend: {self.backend}"
            }

        if result["status"] == "ERROR":
            return {
                "status": "ERROR",
                "pdb_path": None,
                "confidence_global": 0.0,
                "confidence_per_residue": [],
                "source": self.backend,
                "cached": False,
                "error": result["error"]
            }

        pdb_path = save_to_cache(
            sequence=sequence,
            pdb_string=result["pdb_string"],
            confidence_global=result["confidence_global"],
            confidence_per_residue=result["confidence_per_residue"],
            source=self.backend,
            cache_dir=self.cache_dir
        )

        return {
            "status": "SUCCESS",
            "pdb_path": pdb_path,
            "confidence_global": result["confidence_global"],
            "confidence_per_residue": result["confidence_per_residue"],
            "source": self.backend,
            "cached": False,
            "error": None
        }


def main():
    parser = argparse.ArgumentParser(description="Structure Prediction (Entry 028)")
    parser.add_argument("--sequence", required=True)
    parser.add_argument("--backend", default="esmfold")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    provider = StructureProvider(backend=args.backend)
    result = provider.predict(args.sequence)

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        if result["status"] == "SUCCESS":
            print(f"[+] Structure: {result['pdb_path']}")
            print(f"[+] Confidence: {result['confidence_global']:.1%}")
        else:
            print(f"[!] Failed: {result['error']}")


if __name__ == "__main__":
    main()
