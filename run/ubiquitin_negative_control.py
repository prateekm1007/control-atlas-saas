import yaml
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

with open("config/ccs_v0_1.yaml") as f:
    cfg = yaml.safe_load(f)

u = mda.Universe("data/raw/ubq.pdb", "data/raw/ubq.xtc")

# all-to-all pairs (small protein)
ca = u.select_atoms("name CA")
resids = list(ca.resids)
pairs = [(i, j) for i in resids for j in resids if i < j]

dist_series = compute_distance_series(u, pairs)

mi_vals = {}
for (i, j), sig in dist_series.items():
    mi = mutual_information(sig, sig, bins=cfg["mutual_information"]["bins"])
    noise = scrambled_null(mutual_information, sig, sig, cfg["mutual_information"]["shuffles"])
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
    target="UBIQUITIN",
    variant="WT",
    ccs_map=ccs_map,
    provenance={
        "Algorithm": "CCS v0.1",
        "ForceField": "AMBER",
        "Sampling": "NegativeControl"
    }
)

export_edges(entry.top_edges(cfg["thresholds"]["min_ccs_render"]),
             "data/processed/negative_control_ubq_edges.csv")
