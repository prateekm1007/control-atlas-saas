def compute_ccs(mi_vals, path_lengths, rmsf, cfg):
    raw = {}
    median_rmsf = sorted(rmsf.values())[len(rmsf)//2]

    for (i, j), mi in mi_vals.items():
        if (i, j) not in path_lengths:
            continue
        if rmsf[j] / median_rmsf < cfg["thresholds"]["rigid_rmsf_ratio"]:
            continue
        raw[(i, j)] = (
            cfg["weights"]["w_motion"] * mi +
            cfg["weights"]["w_path"] * (1.0 / path_lengths[(i, j)])
        )

    if not raw:
        return {}

    mn, mx = min(raw.values()), max(raw.values())
    return {k: (v - mn) / (mx - mn) for k, v in raw.items()}
