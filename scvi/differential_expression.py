import numpy as np
import torch


def get_statistics(vae, data_loader, M_sampling=2, M_permutation=100, permutation=False):
    """
    Output average over statistics in a symmetric way (a against b)
    forget the sets if permutation is True
    :param vae: The vae model
    :param data_loader:
    :param M_sampling: 200 - default value in Romain's code
    :param M_permutation: 10000 - default value in Romain's code
    :param permutation:
    :return: A 1-d vector of statistics of size n_genes
    """
    # Compute sample rate for the whole dataset ?
    px_scales = []
    all_labels = []
    for sample_batch, _, _, batch_index, labels in data_loader:
        sample_batch = sample_batch.repeat(1, M_sampling).view(-1, sample_batch.size(
            1))  # sample_batch.repeat(1, sample_batch)
        batch_index = batch_index.repeat(1, M_sampling).view(-1, 1)
        labels = labels.repeat(1, M_sampling).view(-1, 1)
        if torch.cuda.is_available():
            sample_batch = sample_batch.cuda(async=True)
            batch_index = batch_index.cuda(async=True)
        px_scales += [vae.get_sample_rate(sample_batch, batch_index)]
        all_labels += [labels]

    couple_celltypes = (2, 4)  # the couple types on which to study DE
    px_scale = torch.cat(px_scales)
    all_labels = torch.cat(all_labels)
    sample_rate_a = px_scale[all_labels == couple_celltypes[0]].view(-1, px_scale.size(1)).cpu().data.numpy()
    sample_rate_b = px_scale[all_labels == couple_celltypes[1]].view(-1, px_scale.size(1)).cpu().data.numpy()

    # agregate dataset
    samples = np.vstack((sample_rate_a, sample_rate_b))

    # prepare the pairs for sampling
    list_1 = list(np.arange(sample_rate_a.shape[0]))
    list_2 = list(sample_rate_a.shape[0] + np.arange(sample_rate_b.shape[0]))
    if not permutation:
        # case1: no permutation, sample from A and then from B
        u, v = np.random.choice(list_1, size=M_permutation), np.random.choice(list_2, size=M_permutation)
    else:
        # case2: permutation, sample from A+B twice
        u, v = (np.random.choice(list_1 + list_2, size=M_permutation),
                np.random.choice(list_1 + list_2, size=M_permutation))

    # then constitutes the pairs
    first_set = samples[u]
    second_set = samples[v]

    res = np.mean(first_set >= second_set, 0)
    res = np.log(res) - np.log(1 - res)
    return res
