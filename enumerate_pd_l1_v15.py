anchors = ["YWPG"]
scaffolds = ["helical", "beta_hairpin"]
buffers = ["GG"]

library = []

for b in buffers:
    for a in anchors:
        for s in scaffolds:
            library.append({
                "target": "PD_L1",
                "buffer": b,
                "anchor": a,
                "scaffold": s
            })

for i, m in enumerate(library, 1):
    print(f"{i:03d}", m)
