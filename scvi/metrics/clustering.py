import numpy as np
import scipy
import torch
from scipy.stats import itemfreq, entropy
from sklearn.neighbors import NearestNeighbors


def get_latent_mean(vae, data_loader):
    return get_latent(vae, data_loader)


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


def entropy_batch_mixing(latent_space, batches, n_neighbors=50, n_pools=50, n_samples_per_pool=100):
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
