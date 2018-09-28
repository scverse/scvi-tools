import numpy as np
import scipy
import torch
from scipy.stats import itemfreq, entropy
from sklearn.cluster import KMeans
from sklearn.neighbors import KNeighborsClassifier
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
    if prediction_algorithm == 'KMeans':
        labels_pred = KMeans(n_labels, n_init=200).fit_predict(latent)  # n_jobs>1 ?
    elif prediction_algorithm == 'knn':
        neigh = KNeighborsClassifier(n_neighbors=10)
        neigh = neigh.fit(latent, labels)
        labels_pred = neigh.predict(latent)
    elif prediction_algorithm == 'gmm':
        gmm = GMM(n_labels, covariance_type='diag', n_init=200)
        gmm.fit(latent)
        labels_pred = gmm.predict(latent)
    else:
        print('algorithm not included: choose from KMeans, knn, or gmm')
    return {
        'asw': silhouette_score(latent, labels),
        'nmi': NMI(labels, labels_pred),
        'ari': ARI(labels, labels_pred),
        'uca': unsupervised_clustering_accuracy(labels, labels_pred)[0]
    }


def get_latent_mean(vae, data_loader):
    latents, batch_indices, labels = get_latents(vae, data_loader)
    return latents[0], batch_indices, labels


get_latent = get_latent_mean


def get_latents(vae, data_loader):
    latents = [[] for _ in range(vae.n_latent_layers)]
    batch_indices = []
    labels = []
    for tensors in data_loader:
        sample_batch, local_l_mean, local_l_var, batch_index, label = tensors
        latents_ = vae.get_latents(sample_batch, y=label)
        latents = [l + [l_] for l, l_ in zip(latents, latents_)]
        labels += [label]
        batch_indices += [batch_index]

    latents = [np.array(torch.cat(l)) for l in latents]
    labels = np.array(torch.cat(labels)).ravel()
    batch_indices = np.array(torch.cat(batch_indices))
    return latents, batch_indices, labels


def get_latents_with_predictions(vae, data_loader):
    latents = [[]] * vae.n_latent_layers
    batch_indices = []
    predictions = []
    labels = []
    for tensorlist in data_loader:
        sample_batch, local_l_mean, local_l_var, batch_index, label = tensorlist
        sample_batch = sample_batch.type(torch.float32)
        latents_ = vae.get_latents(sample_batch, label)
        predictions += [vae.classify(sample_batch)]
        latents = [l + [l_] for l, l_ in zip(latents, latents_)]
        labels += [label]
        batch_indices += [batch_index]

    latents = [np.array(torch.cat(l)) for l in latents]
    labels = np.array(torch.cat(labels)).ravel()
    batch_indices = np.array(torch.cat(batch_indices))
    predictions = np.array(torch.cat(predictions))
    return latents, batch_indices, labels, predictions


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
    jaccard_index = len(set_1.intersection(set_2))/len(set_1.union(set_2))
    return spearman_correlation, fold_enrichment,jaccard_index


def entropy_from_indices(indices):
    return entropy(np.array(itemfreq(indices)[:, 1].astype(np.int32)))


def entropy_batch_mixing(latent_space, batches, n_neighbors=50, n_pools=50, n_samples_per_pool=100, max_number=500):
    n_samples = len(latent_space)
    keep_idx = np.random.choice(np.arange(n_samples), size=min(len(latent_space), max_number), replace=False)
    latent_space, batches = latent_space[keep_idx], batches[keep_idx]

    batches = batches.ravel()
    nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(latent_space)
    indices = nbrs.kneighbors(latent_space, return_distance=False)[:, 1:]
    batch_indices = np.vectorize(lambda i: batches[i])(indices)
    entropies = np.apply_along_axis(entropy_from_indices, axis=1, arr=batch_indices)

    # average n_pools entropy results where each result is an average of n_samples_per_pool random samples.
    if n_pools == 1:
        score = np.mean(entropies)
    else:
        score = np.mean([
            np.mean(entropies[np.random.choice(len(entropies), size=n_samples_per_pool)])
            for _ in range(n_pools)
        ])

    return score


def entropy_from_indices(indices):
    return entropy(np.array(itemfreq(indices)[:, 1].astype(np.int32)))
