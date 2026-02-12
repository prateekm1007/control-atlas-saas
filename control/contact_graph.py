import networkx as nx
from scipy.spatial.distance import euclidean

def build_contact_graph(universe, cutoff=8.0):
    ca = universe.select_atoms("name CA")
    pos = ca.positions
    resids = ca.resids

    G = nx.Graph()
    for i in range(len(resids)):
        for j in range(i + 1, len(resids)):
            if euclidean(pos[i], pos[j]) <= cutoff:
                G.add_edge(resids[i], resids[j])
    return G
