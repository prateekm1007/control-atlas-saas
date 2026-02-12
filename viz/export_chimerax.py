def export_edges(edges, outfile):
    with open(outfile, "w") as f:
        for (i, j), ccs in edges.items():
            f.write(f"{i},{j},{ccs:.3f}\n")
