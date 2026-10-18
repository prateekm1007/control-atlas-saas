variants = []

# Helical variants (length tuning)
for turns in [2, 3, 4]:
    variants.append({
        "id": f"PDL1_V15_H{turns}",
        "buffer": "GG",
        "anchor": "YWPG",
        "scaffold": f"helical_{turns}turn"
    })

# Beta hairpin variants (turn chemistry)
for turn in ["GP", "NG", "DG"]:
    variants.append({
        "id": f"PDL1_V15_B{turn}",
        "buffer": "GG",
        "anchor": "YWPG",
        "scaffold": f"beta_hairpin_{turn}"
    })

for v in variants:
    print(v)
