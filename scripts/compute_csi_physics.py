#!/usr/bin/env python3
import json
import argparse
import numpy as np
from pathlib import Path

def compute_physics_csi(ledger, chassis_id):
    relevant = [
        e for e in ledger
        if chassis_id in str(e.get("motif", "")) or chassis_id in str(e.get("id", ""))
    ]
    targets = {e.get("target") for e in relevant if e.get("target")}
    if len(targets) < 2:
        raise ValueError("Physics-CSI requires ≥2 targets")

    energies = [float(e["energy"]) for e in relevant if isinstance(e.get("energy"), (int, float))]

    # Energetic deviation (normalized heuristic: std_dev / 5000)
    e_dev = np.std(energies) / 5000 if len(energies) >= 2 else 0.75

    return {
        "chassis_id": chassis_id,
        "targets": sorted(list(targets)),
        "energy_std": round(np.std(energies), 2) if energies else None,
        "energy_norm": round(e_dev, 3) if e_dev is not None else None,
        "note": "Physics-derived CSI component; LAW-131 gate remains authoritative"
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ledger", required=True)
    parser.add_argument("--chassis_id", required=True)
    args = parser.parse_args()

    with open(args.ledger) as f:
        content = f.read()
        import re
        matches = re.findall(r'\{.*?\}', content, re.DOTALL)
        ledger = [json.loads(m) for m in matches]

    result = compute_physics_csi(ledger, args.chassis_id)
    print(json.dumps(result, indent=2))
