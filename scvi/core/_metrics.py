import logging
from typing import Tuple, Union

import numpy as np
import scipy
import torch
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
)
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors

logger = logging.getLogger(__name__)


def nearest_neighbor_overlap(x1, x2, k=100):
    """
    Compute the overlap between the k-nearest neighbor graph of x1 and x2.

    Using Spearman correlation of the adjacency matrices.
    Compute the overlap fold enrichment between the protein and mRNA-based cell 100-nearest neighbor
        graph and the Spearman correlation of the adjacency matrices.
    """
    if len(x1) != len(x2):
        raise ValueError("len(x1) != len(x2)")
    n_samples = len(x1)
    k = min(k, n_samples - 1)
    nne = NearestNeighbors(n_neighbors=k + 1)  # "n_jobs=8
    nne.fit(x1)
    kmatrix_1 = nne.kneighbors_graph(x1) - scipy.sparse.identity(n_samples)
    nne.fit(x2)
    kmatrix_2 = nne.kneighbors_graph(x2) - scipy.sparse.identity(n_samples)

    # 1 - spearman correlation from knn graphs
    spearman_correlation = scipy.stats.spearmanr(
        kmatrix_1.A.flatten(), kmatrix_2.A.flatten()
    )[0]
    # 2 - fold enrichment
    set_1 = set(np.where(kmatrix_1.A.flatten() == 1)[0])
    set_2 = set(np.where(kmatrix_2.A.flatten() == 1)[0])
    fold_enrichment = (
        len(set_1.intersection(set_2))
        * n_samples ** 2
        / (float(len(set_1)) * len(set_2))
    )
    return spearman_correlation, fold_enrichment


def unsupervised_clustering_accuracy(
    y: Union[np.ndarray, torch.Tensor], y_pred: Union[np.ndarray, torch.Tensor]
) -> tuple:
    """Unsupervised Clustering Accuracy."""
    if len(y_pred) != len(y):
        raise ValueError("len(y_pred) != len(y)")
    u = np.unique(np.concatenate((y, y_pred)))
    n_clusters = len(u)
    mapping = dict(zip(u, range(n_clusters)))
    reward_matrix = np.zeros((n_clusters, n_clusters), dtype=np.int64)
    for y_pred_, y_ in zip(y_pred, y):
        if y_ in mapping:
            reward_matrix[mapping[y_pred_], mapping[y_]] += 1
    cost_matrix = reward_matrix.max() - reward_matrix
    row_assign, col_assign = linear_sum_assignment(cost_matrix)

    # Construct optimal assignments matrix
    row_assign = row_assign.reshape((-1, 1))  # (n,) to (n, 1) reshape
    col_assign = col_assign.reshape((-1, 1))  # (n,) to (n, 1) reshape
    assignments = np.concatenate((row_assign, col_assign), axis=1)

    optimal_reward = reward_matrix[row_assign, col_assign].sum() * 1.0
    return optimal_reward / y_pred.size, assignments


def knn_purity(latent, label, n_neighbors=30):
    nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(latent)
    indices = nbrs.kneighbors(latent, return_distance=False)[:, 1:]
    neighbors_labels = np.vectorize(lambda i: label[i])(indices)

    # pre cell purity scores
    scores = ((neighbors_labels - label.reshape(-1, 1)) == 0).mean(axis=1)
    res = [
        np.mean(scores[label == i]) for i in np.unique(label)
    ]  # per cell-type purity

    return np.mean(res)


@torch.no_grad()
def clustering_scores(
    self, adata, latent, labels, prediction_algorithm: str = "knn"
) -> Tuple:
    if adata.uns["scvi_summary_stats"]["n_labels"] > 1:
        if prediction_algorithm == "knn":
            labels_pred = KMeans(
                self.dataset.adata.uns["scvi_summary_stats"]["n_labels"],
                n_init=200,
            ).fit_predict(latent)
        elif prediction_algorithm == "gmm":
            gmm = GaussianMixture(
                self.dataset.adata.uns["scvi_summary_stats"]["n_labels"]
            )
            gmm.fit(latent)
            labels_pred = gmm.predict(latent)

        asw_score = silhouette_score(latent, labels)
        nmi_score = normalized_mutual_info_score(labels, labels_pred)
        ari_score = adjusted_rand_score(labels, labels_pred)
        uca_score = unsupervised_clustering_accuracy(labels, labels_pred)[0]
        logger.debug(
            "Clustering Scores:\nSilhouette: %.4f\nNMI: %.4f\nARI: %.4f\nUCA: %.4f"
            % (asw_score, nmi_score, ari_score, uca_score)
        )
        return asw_score, nmi_score, ari_score, uca_score
