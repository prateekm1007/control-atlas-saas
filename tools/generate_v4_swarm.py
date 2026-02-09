# v4.0 Ionic-Claw Saturation Swarm Generator

scaffold = "KAAAKKAAAKPPPP"  # Validated Ionic Needle spine
warheads = ["YWPG", "YWPA", "YWPTG", "YWPSG"]
side_hooks = ["", "W", "S", "D", "K"]

v4_swarm = []

for w in warheads:
    for h in side_hooks:
        # Insert hook into the spine
        motif = f"{scaffold[:5]}{h}{scaffold[6:]}{w}"
        v4_swarm.append(motif)

with open("artifacts/v4_swarm_manifest.txt", "w") as f:
    for m in v4_swarm:
        f.write(m + "\n")

print(f"✅ v4.0 Swarm Manifest Created: {len(v4_swarm)} motifs")
