#!/usr/bin/env python3

import argparse
import json
import os
import statistics
import sys

def energy_proxy(name):
    # Deterministic placeholder (replace with MM-GBSA later)
    return (len(name) % 17) / 10.0

def strain_proxy(name):
    # RMSD / clash proxy
    return (len(name) % 11) / 10.0

def interface_proxy(name):
    # Epitope coverage proxy
    return (len(name) % 7) / 10.0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", required=True)
    parser.add_argument("--chassis_id", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if not os.path.isdir(args.inputs):
        print(f"❌ Input directory not found: {args.inputs}")
        sys.exit(1)

    files = sorted(os.listdir(args.inputs))
    if len(files) < 2:
        print("❌ CSI requires ≥2 motifs evaluated on the same chassis")
        sys.exit(1)

    energies, strains, interfaces = [], [], []

    for f in files:
        energies.append(energy_proxy(f))
        strains.append(strain_proxy(f))
        interfaces.append(interface_proxy(f))

    def norm(vals):
        return min(1.0, statistics.mean(vals) / 2.0)

    E_dev = norm(energies)
    S_strain = norm(strains)
    I_penalty = norm(interfaces)

    CSI = (0.4 * E_dev) + (0.4 * S_strain) + (0.2 * I_penalty)

    result = {
        "law": "LAW-131",
        "chassis_id": args.chassis_id,
        "targets_evaluated": len(files),
        "E_dev": round(E_dev, 3),
        "S_strain": round(S_strain, 3),
        "I_penalty": round(I_penalty, 3),
        "CSI": round(CSI, 3),
        "interpretation": (
            "GREEN/BLUE" if CSI < 1.2 else
            "YELLOW" if CSI < 1.5 else
            "RED"
        )
    }

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(result, f, indent=2)

    print(f"✅ CSI computed: {result['CSI']} ({result['interpretation']})")

if __name__ == "__main__":
    main()
