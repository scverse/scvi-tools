import numpy as np
import torch

from scvi.utils import to_cuda, no_grad, eval_modules


@no_grad()
@eval_modules()
def imputation(model, data_loader, rate=0.1, use_cuda=True):
    distance_list = torch.FloatTensor([])
    for tensors in data_loader:
        tensors = to_cuda(tensors, use_cuda=use_cuda)
        sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors
        sample_batch = sample_batch.type(torch.float32)
        dropout_batch = sample_batch.clone()
        indices = torch.nonzero(dropout_batch)
        i, j = indices[:, 0], indices[:, 1]
        ix = torch.LongTensor(np.random.choice(range(len(i)), int(np.floor(rate * len(i))), replace=False))
        dropout_batch[i[ix], j[ix]] *= 0

        if use_cuda:
            ix, i, j = to_cuda([ix, i, j], async=False)
        px_rate = model.get_sample_rate(dropout_batch, batch_index=batch_index, y=labels)
        distance_list = torch.cat([distance_list, torch.abs(px_rate[i[ix], j[ix]] - sample_batch[i[ix], j[ix]]).cpu()])
    return distance_list
