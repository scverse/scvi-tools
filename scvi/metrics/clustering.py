import numpy as np
import torch

from scvi.utils import to_cuda, no_grad, eval_modules


@no_grad()
def get_latents(vae, data_loader):
    latents = [[]] * vae.n_latent_layers
    batch_indices = []
    labels = []
    for tensorlist in data_loader:
        if vae.use_cuda:
            tensorlist = to_cuda(tensorlist)
        sample_batch, local_l_mean, local_l_var, batch_index, label = tensorlist
        sample_batch = sample_batch.type(torch.float32)
        latents_ = vae.get_latents(sample_batch, label)
        latents = [l + [l_] for l, l_ in zip(latents, latents_)]
        labels += [label]
        batch_indices += [batch_index]

    latents = [np.array(torch.cat(l)) for l in latents]
    labels = np.array(torch.cat(labels)).ravel()
    batch_indices = np.array(torch.cat(batch_indices))
    return latents, batch_indices, labels


@eval_modules()
def get_latent_mean(vae, data_loader):
    return get_latent(vae, data_loader)


def get_latent(vae, data_loader):
    latents, batch_indices, labels = get_latents(vae, data_loader)
    return latents[0], batch_indices, labels


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
