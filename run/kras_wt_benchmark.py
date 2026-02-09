import json, yaml
import MDAnalysis as mda
import networkx as nx

from signals.distances import compute_distance_series
from signals.mutual_info import mutual_information
from signals.rmsf import compute_rmsf
from control.contact_graph import build_contact_graph
from control.ccs import compute_ccs
from control.negative_controls import scrambled_null
from atlas.entry import AtlasEntry
from viz.export_chimerax import export_edges

with open("config/regions_kras.json") as f:
    regions = json.load(f)

with open("config/ccs_v0_1.yaml") as f:
    cfg = yaml.safe_load(f)

u = mda.Universe("data/raw/kras_wt.pdb", "data/raw/kras_wt.xtc")

pairs = []
for i in range(regions["Mutation"][0], regions["Switch_II"][1] + 1):
    for j in range(regions["Switch_II"][0], regions["Switch_II"][1] + 1):
        pairs.append((i, j))

dist_series = compute_distance_series(u, pairs)

mi_vals = {}
for (i, j), sig in dist_series.items():
    ref = dist_series[(i, j)]
    mi = mutual_information(sig, ref, bins=cfg["mutual_information"]["bins"])
    noise = scrambled_null(mutual_information, sig, ref, cfg["mutual_information"]["shuffles"])
    if mi > noise:
        mi_vals[(i, j)] = mi

G = build_contact_graph(u, cfg["contact_graph"]["cutoff_angstrom"])
path_lengths = dict(nx.all_pairs_shortest_path_length(G))

flat_paths = {}
for (i, j) in mi_vals:
    if i in path_lengths and j in path_lengths[i]:
        flat_paths[(i, j)] = path_lengths[i][j]

rmsf = compute_rmsf(u)

ccs_map = compute_ccs(mi_vals, flat_paths, rmsf, cfg)

entry = AtlasEntry(
    target="KRAS",
    variant="WT",
    ccs_map=ccs_map,
    provenance={
        "Algorithm": "CCS v0.1",
        "ForceField": "AMBER",
        "Sampling": "Benchmark"
    }
)

export_edges(entry.top_edges(cfg["thresholds"]["min_ccs_render"]),
             "data/processed/atlas_entry_001_wt_edges.csv")
