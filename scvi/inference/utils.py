import numpy as np
import leidenalg
import igraph as ig

from sklearn.neighbors import kneighbors_graph


def louvain_clusters(latent, k=10, rands=0):
    nn_matrix = kneighbors_graph(latent, k)
    rows, cols = np.where(nn_matrix.todense() == 1)
    edges = [(row, col) for row, col in zip(rows, cols)]
    g = ig.Graph()
    g.add_vertices(latent.shape[0])
    g.add_edges(edges)
    res = leidenalg.find_partition(g, leidenalg.ModularityVertexPartition, seed=rands)
    clusters = np.asarray(res.membership)
    return clusters
