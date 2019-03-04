import numpy as np

import leidenalg
import igraph as ig

"""
try:
    import igraph as ig
except DeprecationWarning:
    print("Deprecation Warning casted as error by louvain")
    pass
"""
#import louvain

from sklearn.neighbors import kneighbors_graph


def louvain_clusters(latent, k=10, rands=0):
    nn_matrix = kneighbors_graph(latent, k)
    rows, cols = np.where(nn_matrix.todense() == 1)
    edges = [(row, col) for row, col in zip(rows, cols)]
    g = ig.Graph()
    g.add_vertices(latent.shape[0])
    g.add_edges(edges)
    #louvain.set_rng_seed(rands)
    res = leidenalg.find_partition(g, leidenalg.ModularityVertexPartition, seed=rands)
    #print(res)
    #print(res.__class__)
    #res = louvain.find_partition(g, louvain.ModularityVertexPartition)
    clusters = np.asarray(res.membership)
    return clusters
