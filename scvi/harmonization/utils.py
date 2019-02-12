from numpy import loadtxt
from scipy.io import mmread
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scvi.dataset.dataset import GeneExpressionDataset
import igraph as ig
import louvain
from sklearn.neighbors import kneighbors_graph


def get_matrix_from_dir(dirname, storage='../'):
    geneid = loadtxt(storage + dirname + '/genes.tsv', dtype='str', delimiter="\t")
    cellid = loadtxt(storage + dirname + '/barcodes.tsv', dtype='str', delimiter="\t")
    count = mmread(storage + dirname + '/matrix.mtx')
    return count, geneid, cellid


class ExpressionMatrix(GeneExpressionDataset):
    def __init__(self, filename, save_path):
        self.save_path = save_path
        self.filename = filename
        count, barcode, genenames = self.preprocess()
        super(ExpressionMatrix, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(count),
            gene_names=genenames)

    def preprocess(self):
        data = pd.read_csv(self.save_path + self.filename, sep="\t", low_memory=False)
        barcode = data.columns.values
        genenames = np.asarray(data.index)
        count = np.asarray(data).astype('int64').T
        return count, barcode, genenames


def louvain_clusters(latent, k=10, rands=0):
    nn_matrix = kneighbors_graph(latent, k)
    rows, cols = np.where(nn_matrix.todense() == 1)
    edges = [(row, col) for row, col in zip(rows, cols)]
    g = ig.Graph()
    g.add_vertices(latent.shape[0])
    g.add_edges(edges)
    louvain.set_rng_seed(rands)
    res = louvain.find_partition(g, louvain.ModularityVertexPartition)
    clusters = np.asarray(res.membership)
    return clusters


def plot_marker_genes(latent_u, count, genenames, markers):
    nrow = (len(markers) // 3 + 1)
    figh = nrow * 4
    plt.figure(figsize=(10, figh))
    for i, x in enumerate(markers):
        exprs = count[:, genenames == x].ravel()
        idx = (exprs > 0)
        plt.subplot(nrow, 3, (i + 1))
        plt.scatter(latent_u[:, 0], latent_u[:, 1], c='lightgrey', edgecolors='none', s=5)
        plt.scatter(latent_u[idx, 0], latent_u[idx, 1], c=exprs[idx], cmap=plt.get_cmap('viridis_r'),
                    edgecolors='none', s=3)
        plt.title(x)
        plt.tight_layout()
    plt.show()


def plot_marker_genes_compare(latent_u, count, genenames, markers, subset):
    nrow = len(markers)
    figh = nrow * 4
    plt.figure(figsize=(8, figh))
    for i, x in enumerate(markers):
        exprs = count[:, genenames == x].ravel()
        idx = (exprs > 0)
        plt.subplot(nrow, 2, (i*2 + 1))
        plt.scatter(latent_u[subset, 0], latent_u[subset, 1], c='lightgrey', edgecolors='none', s=5)
        plt.scatter(latent_u[idx, 0][subset[idx]], latent_u[idx, 1][subset[idx]], c=exprs[idx][subset[idx]],
                    cmap=plt.get_cmap('viridis_r'), edgecolors='none', s=3)
        plt.title(x+' control')
        plt.tight_layout()
        plt.subplot(nrow, 2, (i * 2 + 2))
        subset = np.asarray([not x for x in subset])
        plt.scatter(latent_u[subset, 0], latent_u[subset, 1], c='lightgrey', edgecolors='none', s=5)
        plt.scatter(latent_u[idx, 0][subset[idx]], latent_u[idx, 1][subset[idx]], c=exprs[idx][subset[idx]],
                    cmap=plt.get_cmap('viridis_r'),  edgecolors='none', s=3)
        plt.title(x+' stimulated')
    plt.show()


def find_markers(deres, absthres, relthres, ngenes, clusterid):
    allgenes = []
    for i, x in enumerate(deres):
        markers = x.loc[(x['mean1'] > absthres) & (x['mean2'] / x['mean2'] < relthres)]
        ngenes = np.min([len(markers), ngenes])
        markers = markers[:ngenes]
        markers['clusterid'] = np.repeat(clusterid[i], ngenes)
        allgenes.append(markers)
    markers = pd.concat(allgenes)
    return markers
