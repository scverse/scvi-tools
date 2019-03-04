import numpy as np
import igraph as ig
import louvain
from sklearn.neighbors import kneighbors_graph



def louvain_clusters(latent, k=10, rands=0, mutual=False):
    nn_matrix = kneighbors_graph(latent, k)
    rows, cols = nn_matrix.nonzero()
    if mutual == True:
        edges = [[row, col] if row < col else (col, row) for row, col in zip(rows, cols)]
        edges = np.asarray(edges)
        unique_edges, edges_count = np.unique(edges, return_counts=True, axis=0)
        edges = unique_edges[edges_count == 2]
    else:
        edges = [(row, col) for row, col in zip(rows, cols)]
    g = ig.Graph()
    g.add_vertices(latent.shape[0])
    g.add_edges(edges)
    louvain.set_rng_seed(rands)
    res = louvain.find_partition(g, louvain.ModularityVertexPartition)
    clusters = np.asarray(res.membership)
    return clusters

