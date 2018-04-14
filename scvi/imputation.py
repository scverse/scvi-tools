import numpy as np

import torch


def imputation(vae, data_loader, rate=0.1):
    distance_list = torch.FloatTensor([])
    for sample_batch, local_l_mean, local_l_var, batch_index in data_loader:
        if torch.cuda.is_available():
            sample_batch = sample_batch.cuda()
            batch_index = batch_index.cuda()
        dropout_batch = sample_batch.clone()
        indices = torch.nonzero(dropout_batch)
        i, j = indices[:, 0], indices[:, 1]
        ix = torch.LongTensor(np.random.choice(range(len(i)), int(np.floor(rate * len(i))), replace=False))
        if torch.cuda.is_available:
            ix = ix.cuda()
        dropout_batch[i[ix], j[ix]] *= 0
        _, _, px_rate, _, _, _, _, _ = vae(dropout_batch, batch_index)
        distance_list = torch.cat([distance_list,
                                   torch.abs(px_rate[i[ix], j[ix]].data - sample_batch[i[ix], j[ix]]).cpu()])
    return torch.median(distance_list)
