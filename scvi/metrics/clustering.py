import numpy as np
import scipy
import torch
from scipy.stats import itemfreq, entropy
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score as ARI
from sklearn.metrics import normalized_mutual_info_score as NMI
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture as GMM
from sklearn.neighbors import NearestNeighbors
from sklearn.utils.linear_assignment_ import linear_assignment


def unsupervised_clustering_accuracy(y, y_pred):
    """
    Unsupervised Clustering Accuracy
    """
    assert len(y_pred) == len(y)
    u = np.unique(np.concatenate((y, y_pred)))
    n_clusters = len(u)
    mapping = dict(zip(u, range(n_clusters)))
    reward_matrix = np.zeros((n_clusters, n_clusters), dtype=np.int64)
    for y_pred_, y_ in zip(y_pred, y):
        if y_ in mapping:
            reward_matrix[mapping[y_pred_], mapping[y_]] += 1
    cost_matrix = reward_matrix.max() - reward_matrix
    ind = linear_assignment(cost_matrix)
    return sum([reward_matrix[i, j] for i, j in ind]) * 1.0 / y_pred.size, ind


def clustering_scores(latent, labels, prediction_algorithm='knn', n_labels=None):
    if n_labels is not None:
        n_labels = len(np.unique(labels))

    if prediction_algorithm == 'knn':
        labels_pred = KMeans(n_labels, n_init=200).fit_predict(latent)  # n_jobs>1 ?
    elif prediction_algorithm == 'gmm':
        gmm = GMM(n_labels, covariance_type='diag', n_init=200)
        gmm.fit(latent)
        labels_pred = gmm.predict(latent)

    return {
        'asw': silhouette_score(latent, labels),
        'nmi': NMI(labels, labels_pred),
        'ari': ARI(labels, labels_pred),
        'uca': unsupervised_clustering_accuracy(labels, labels_pred)[0]
    }


def get_latent(vae, data_loader):
    latent = []
    batch_indices = []
    labels = []
    for tensors in data_loader:
        sample_batch, local_l_mean, local_l_var, batch_index, label = tensors
        latent += [vae.sample_from_posterior_z(sample_batch, y=label)]
        batch_indices += [batch_index]
        labels += [label]
    return np.array(torch.cat(latent)), np.array(torch.cat(batch_indices)), np.array(torch.cat(labels)).ravel()


def select_indices_evenly(n_samples, categorical_indices):
    r'''
    In t_sne or entropy_batch_mixing, unbalance between the size of the different batches can impair the results.
    Instead of uniformly random selection of indices over the whole dataset, we split this selection in each
    categorical subsets, so as to provide equally weighted indices for each cluster.
    :param n_samples: The number of samples to select per category.
    :param categorical_indices: a numpy array of shape (n_samples, 1) or n_samples, giving the categories.
    :return: a numpy array of indices to select
    '''
    categorical_indices = categorical_indices.ravel()
    indices = []
    for i in np.unique(categorical_indices):
        indices_i = np.where(categorical_indices == i)[0]
        indices += [indices_i[np.random.permutation(len(indices_i))][:n_samples]]
    return np.concatenate(indices, axis=0)


def nn_overlap(X1, X2, k=100):
    nne = NearestNeighbors(n_neighbors=k + 1, n_jobs=8)
    assert len(X1) == len(X2)
    n_samples = len(X1)
    nne.fit(X1)
    kmatrix_1 = nne.kneighbors_graph(X1) - scipy.sparse.identity(n_samples)
    nne.fit(X2)
    kmatrix_2 = nne.kneighbors_graph(X2) - scipy.sparse.identity(n_samples)

    # 1 - spearman correlation from knn graphs
    spearman_correlation = scipy.stats.spearmanr(kmatrix_1.A.flatten(), kmatrix_2.A.flatten())[0]
    # 2 - fold enrichment
    set_1 = set(np.where(kmatrix_1.A.flatten() == 1)[0])
    set_2 = set(np.where(kmatrix_2.A.flatten() == 1)[0])
    fold_enrichment = len(set_1.intersection(set_2)) * n_samples ** 2 / (float(len(set_1)) * len(set_2))
    return spearman_correlation, fold_enrichment


def entropy_from_indices(indices):
    return entropy(np.array(itemfreq(indices)[:, 1].astype(np.int32)))


def entropy_batch_mixing(latent_space, batches, n_neighbors=50):
    batches = batches.ravel()
    nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(latent_space)
    indices = nbrs.kneighbors(latent_space, return_distance=False)[:, 1:]
    batch_indices = np.vectorize(lambda i: batches[i])(indices)
    return np.mean(np.apply_along_axis(entropy_from_indices, axis=1, arr=batch_indices))
