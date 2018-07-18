import numpy as np
import torch


def imputation(vae, data_loader, rate=0.1):
    distance_list = torch.FloatTensor([])
    for tensors in data_loader:
        sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors
        dropout_batch = sample_batch.clone()
        indices = torch.nonzero(dropout_batch)
        i, j = indices[:, 0], indices[:, 1]
        ix = torch.LongTensor(np.random.choice(range(len(i)), int(np.floor(rate * len(i))), replace=False))
        dropout_batch[i[ix], j[ix]] *= 0

        ix, i, j = [t.to(sample_batch.device) for t in [ix, i, j]]
        px_rate = vae.get_sample_rate(dropout_batch, batch_index=batch_index, y=labels)
        distance_list = torch.cat([distance_list, torch.abs(px_rate[i[ix], j[ix]] - sample_batch[i[ix], j[ix]]).cpu()])
    return distance_list
