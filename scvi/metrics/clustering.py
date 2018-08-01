import numpy as np
import scipy
import torch
from sklearn.neighbors import NearestNeighbors
from sklearn.manifold import TSNE


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


# CLUSTERING METRICS
def entropy_batch_mixing(latent_space, batches, max_number=500):
    # latent space: numpy matrix of size (number_of_cells, latent_space_dimension)
    # with the encoding of the different inputs in the latent space
    # batches: numpy vector with the batch indices of the cells
    n_samples = len(latent_space)
    keep_idx = np.random.choice(np.arange(n_samples), size=min(len(latent_space), max_number), replace=False)
    latent_space, batches = latent_space[keep_idx], batches[keep_idx]

    def entropy(hist_data):
        n_batches = len(np.unique(hist_data))
        if n_batches > 2:
            raise ValueError("Should be only two clusters for this metric")
        frequency = np.mean(hist_data == 1)
        if frequency == 0 or frequency == 1:
            return 0
        return -frequency * np.log(frequency) - (1 - frequency) * np.log(1 - frequency)

    n_samples = latent_space.shape[0]
    distance = np.zeros((n_samples, n_samples))
    neighbors_graph = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(i, n_samples):
            distance[i, j] = distance[j, i] = sum((latent_space[i] - latent_space[j]) ** 2)

    for i, d in enumerate(distance):
        neighbors_graph[i, d.argsort()[:51]] = 1
    kmatrix = neighbors_graph - np.identity(latent_space.shape[0])

    score = 0
    for t in range(50):
        indices = np.random.choice(np.arange(latent_space.shape[0]), size=100)
        score += np.mean([entropy(
            batches[kmatrix[indices].nonzero()[1][kmatrix[indices].nonzero()[0] == i]]
        ) for i in range(100)])
    return score / 50


def get_data(vae, data_loader, mode):
    """ Specific to smFISH/scRNA-seq datasets.
        Yields positional information, and other useful informations
        such as expected frequencies, real values...
        Easier to use a new method for this as vae_fish has some particularities
        that would require to modify a lot the other functions if we didn't create this one
    """
    vae.eval()
    latent = []
    batch_indices = []
    labels = []
    expected_frequencies = []
    values = []
    x_coords = []
    y_coords = []
    for tensors in data_loader:
        if mode == "scRNA":
            sample_batch, local_l_mean, local_l_var, batch_index, label = tensors
            batch_index = torch.zeros_like(batch_index)
        if mode == "smFISH":
            sample_batch, local_l_mean, local_l_var, batch_index, label, x_coord, y_coord = tensors
            batch_index = torch.ones_like(batch_index)
            x_coords += [x_coord]
            y_coords += [y_coord]
        latent += [vae.sample_from_posterior_z(sample_batch, y=label, mode=mode)]
        batch_indices += [batch_index]
        labels += [label]
        expected_frequencies += [vae.get_sample_scale(sample_batch, mode=mode, batch_index=batch_index)]
        values += [sample_batch]
    if mode == "scRNA":
        return np.array(torch.cat(latent)), np.array(torch.cat(batch_indices)), np.array(torch.cat(labels)).ravel(), \
               np.array(torch.cat(expected_frequencies)), np.array(torch.cat(values))
    if mode == "smFISH":
        return np.array(torch.cat(latent)), np.array(torch.cat(batch_indices)), np.array(torch.cat(labels)).ravel(), \
               np.array(torch.cat(expected_frequencies)), np.array(torch.cat(values)),\
               np.array(torch.cat(x_coords)), np.array(torch.cat(y_coords))


def get_common_t_sne(latent_seq, latent_fish, n_samples=1000):
    idx_t_sne_a = np.random.permutation(len(latent_seq))[:n_samples]
    idx_t_sne_b = np.random.permutation(len(latent_fish))[:n_samples]
    full_latent = np.concatenate((latent_seq[idx_t_sne_a, :], latent_fish[idx_t_sne_b, :]))
    if full_latent.shape[1] != 2:
        latent = TSNE().fit_transform(full_latent)
    if latent.shape[0] != len(idx_t_sne_a) + len(idx_t_sne_b):
        print("Be careful! There might be a mistake in the downsampling of the data points")
    return latent[:len(idx_t_sne_a), :], latent[len(idx_t_sne_a):, :], idx_t_sne_a, idx_t_sne_b
