import numpy as np
import torch
from torch.utils.data import DataLoader


# from torch.utils.data.sampler import SubsetRandomSampler

def get_statistics(vae, cortex_train_dataset, M_sampling=2, M_p=10000, permutation=False):
    """
    Output average over statistics in a symmetric way (a against b)
    forget the sets if permutation is True
    """
    labels_train = np.load('data/labels_train.npy')  # cell types

    # Compute sample rate for whole dataset
    data_loader_train = DataLoader(cortex_train_dataset, batch_size=128, shuffle=False, num_workers=1, pin_memory=True)
    px_scales = []
    for sample_batch, _, _, batch_index in data_loader_train:
        sample_batch = sample_batch.repeat(1, M_sampling).view(-1, sample_batch.size(
            1))  # sample_batch.repeat(1, sample_batch)
        batch_index = batch_index.repeat(1, M_sampling).view(-1, 1)
        if torch.cuda.is_available():
            sample_batch = sample_batch.cuda(async=True)
            batch_index = batch_index.cuda(async=True)
        px_scales += [vae.get_sample_rate(sample_batch, batch_index)]

    couple_celltypes = (0, 1)  # the couple types on which to study DE
    px_scale = torch.cat(px_scales).numpy()
    sample_rate_a = px_scale[labels_train == couple_celltypes[0]]
    sample_rate_b = px_scale[labels_train == couple_celltypes[1]]

    # agregate dataset
    samples = np.vstack((sample_rate_a, sample_rate_b))

    # prepare the pairs for sampling
    list_1 = list(np.arange(sample_rate_a.shape[0]))
    list_2 = list(sample_rate_a.shape[0] + np.arange(sample_rate_b.shape[0]))
    if not permutation:
        # case1: no permutation, sample from A and then from B
        u, v = np.random.choice(list_1, size=M_p), np.random.choice(list_2, size=M_p)
    else:
        # case2: permutation, sample from A+B twice
        u, v = (np.random.choice(list_1 + list_2, size=M_p), \
                np.random.choice(list_1 + list_2, size=M_p))

    # then constitutes the pairs
    first_set = samples[u]
    second_set = samples[v]

    res = np.mean(first_set >= second_set, 0)
    res = np.log(res) - np.log(1 - res)
    return res
