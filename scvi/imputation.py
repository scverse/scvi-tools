import numpy as np

import torch


def imputation(vae, data_loader, rate=0.1):
    distance_list = torch.FloatTensor([])
    for sample_batch, local_l_mean, local_l_var, batch_index in data_loader:
        dropout_batch = sample_batch.clone()
        indices = torch.nonzero(dropout_batch)
        i, j = indices[:, 0], indices[:, 1]
        ix = torch.LongTensor(np.random.choice(range(len(i)), int(np.floor(rate * len(i))), replace=False))
        dropout_batch[i[ix], j[ix]] *= 0

        if torch.cuda.is_available():
            batch_index = batch_index.cuda()
            dropout_batch = dropout_batch.cuda()
            sample_batch = sample_batch.cuda()
            distance_list = distance_list.cuda()
            ix = ix.cuda()
            i = i.cuda()
            j = j.cuda()

        _, _, px_rate, _, _, _, _, _ = vae(dropout_batch, batch_index)
        distance_list = torch.cat([distance_list, torch.abs(px_rate[i[ix], j[ix]].data - sample_batch[i[ix], j[ix]])])
    return torch.median(distance_list)
