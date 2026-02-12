import csv

def load_edges(path):
    d = {}
    with open(path) as f:
        for r in csv.reader(f):
            i, j, v = int(r[0]), int(r[1]), float(r[2])
            d[(i, j)] = v
    return d

g12c = load_edges("data/processed/atlas_entry_001_edges.csv")
wt = load_edges("data/processed/atlas_entry_001_wt_edges.csv")

delta = {}
for k, v in g12c.items():
    delta[k] = v - wt.get(k, 0.0)

with open("data/processed/atlas_entry_001_delta.csv", "w") as f:
    for (i, j), v in sorted(delta.items(), key=lambda x: -x[1]):
        if v > 0:
            f.write(f"{i},{j},{v:.3f}\n")
