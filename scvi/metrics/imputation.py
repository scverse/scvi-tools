import numpy as np
import torch

from scvi.utils import to_cuda, no_grad, eval_modules


@no_grad()
@eval_modules()
def imputation(vae, data_loader, rate=0.1):
    distance_list = torch.FloatTensor([])
    for tensorlist in data_loader:
        if vae.use_cuda:
            tensorlist = to_cuda(tensorlist)
        sample_batch, local_l_mean, local_l_var, batch_index, labels = tensorlist
        sample_batch = sample_batch.type(torch.float32)
        dropout_batch = sample_batch.clone()
        indices = torch.nonzero(dropout_batch)
        i, j = indices[:, 0], indices[:, 1]
        ix = torch.LongTensor(np.random.choice(range(len(i)), int(np.floor(rate * len(i))), replace=False))
        dropout_batch[i[ix], j[ix]] *= 0

        if vae.use_cuda:
            ix, i, j = to_cuda([ix, i, j], async=False)
        px_rate = vae.get_sample_rate(dropout_batch, labels, batch_index=batch_index)
        distance_list = torch.cat([distance_list, torch.abs(px_rate[i[ix], j[ix]] - sample_batch[i[ix], j[ix]]).cpu()])
    return torch.median(distance_list)
