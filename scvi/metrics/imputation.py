import copy

import matplotlib.pyplot as plt
import numpy as np
import torch
from scipy.stats import kde
from sklearn import neighbors

from scvi.dataset import GeneExpressionDataset
from scvi.dataset.data_loaders import DataLoaderWrapper


def imputation(infer, name, rate=0.1, n_samples=1, n_epochs=1, corruption="uniform"):
    corrupted_data = copy.deepcopy(infer.gene_dataset.X)

    if corruption == "uniform":  # multiply the entry n with a Ber(0.9) random variable.
        i, j = np.nonzero(corrupted_data)
        ix = np.random.choice(range(len(i)), int(np.floor(rate * len(i))), replace=False)
        i, j = i[ix], j[ix]
        corrupted_data[i, j] *= np.random.binomial(n=np.ones(len(ix), dtype=np.int64), p=0.9)
    elif corruption == "binomial":  # multiply the entry n with a Bin(n, 0.9) random variable.
        i, j = (k.ravel() for k in np.indices(corrupted_data.shape))
        ix = np.random.choice(range(len(i)), int(np.floor(rate * len(i))), replace=False)
        i, j = i[ix], j[ix]
        corrupted_data[i, j] = np.random.binomial(n=corrupted_data[i, j].astype(np.int64), p=0.2)

    infer.gene_dataset = gene_dataset = GeneExpressionDataset(
        *GeneExpressionDataset.get_attributes_from_matrix(
            corrupted_data,
            batch_indices=infer.gene_dataset.batch_indices,
            labels=infer.gene_dataset.labels
        )
    )

    original_data_loaders_loop = infer.data_loaders.loop
    infer.data_loaders.loop = ['corrupted_%s' % s for s in infer.data_loaders.loop]
    original_keys = list(infer.data_loaders.dict.keys())
    for key in original_keys:
        kwargs = copy.copy(infer.data_loaders.kwargs)
        kwargs['collate_fn'] = gene_dataset.collate_fn
        kwargs['sampler'] = copy.copy(infer.data_loaders[key].sampler)
        infer.data_loaders['corrupted_%s' % key] = DataLoaderWrapper(gene_dataset, use_cuda=infer.use_cuda, **kwargs)

    infer.train(n_epochs=n_epochs)
    infer.data_loaders.loop = original_data_loaders_loop

    original_list = []
    imputed_list = []
    batch_size = infer.data_loaders.kwargs["batch_size"] // n_samples
    for tensors, corrupted_tensors in \
        zip(infer.data_loaders[name].sequential(batch_size=batch_size),
            infer.data_loaders['corrupted_%s' % name].sequential(batch_size=batch_size)):
        batch = tensors[0]
        actual_batch_size = batch.size(0)
        dropout_batch, _, _, batch_index, labels = corrupted_tensors
        px_rate = infer.model.get_sample_rate(dropout_batch, batch_index=batch_index, y=labels, n_samples=n_samples)

        indices_dropout = torch.nonzero(batch - dropout_batch)
        i = indices_dropout[:, 0]
        j = indices_dropout[:, 1]

        batch = batch.unsqueeze(0).expand((n_samples, batch.size(0), batch.size(1)))
        original = np.array(batch[:, i, j].view(-1).cpu())
        imputed = np.array(px_rate[:, i, j].view(-1).cpu())

        cells_index = np.tile(np.array(i.cpu()), n_samples)

        original_list += [original[cells_index == i] for i in range(actual_batch_size)]
        imputed_list += [imputed[cells_index == i] for i in range(actual_batch_size)]

    return original_list, imputed_list


def plot_imputation(original, imputed, title="Imputation"):
    y = imputed
    x = original

    ymax = 10
    mask = x < ymax
    x = x[mask]
    y = y[mask]

    mask = y < ymax
    x = x[mask]
    y = y[mask]

    l_minimum = np.minimum(x.shape[0], y.shape[0])

    x = x[:l_minimum]
    y = y[:l_minimum]

    data = np.vstack([x, y])

    plt.figure(figsize=(5, 5))

    axes = plt.gca()
    axes.set_xlim([0, ymax])
    axes.set_ylim([0, ymax])

    nbins = 50

    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
    k = kde.gaussian_kde(data)
    xi, yi = np.mgrid[0:ymax:nbins * 1j, 0:ymax:nbins * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    plt.title(title, fontsize=12)
    plt.ylabel("Imputed counts")
    plt.xlabel('Original counts')

    plt.pcolormesh(yi, xi, zi.reshape(xi.shape), cmap="Reds")

    a, _, _, _ = np.linalg.lstsq(y[:, np.newaxis], x, rcond=-1)
    linspace = np.linspace(0, ymax)
    plt.plot(linspace, a * linspace, color='black')

    plt.plot(linspace, linspace, color='black', linestyle=":")
    plt.show()
    plt.savefig(title + '.png')


def proximity_imputation(real_latent1, normed_gene_exp_1, real_latent2, k=4):
    knn = neighbors.KNeighborsRegressor(k, weights='distance')
    y = knn.fit(real_latent1, normed_gene_exp_1).predict(real_latent2)
    return y
