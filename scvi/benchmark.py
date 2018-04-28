import numpy as np
import torch
from torch.utils.data import DataLoader
from torch.utils.data.sampler import SubsetRandomSampler

from scvi.clustering import entropy_batch_mixing
from scvi.dataset import CortexDataset
from scvi.differential_expression import get_statistics
from scvi.imputation import imputation
from scvi.models import VAE
from scvi.train import train
from scvi.utils import to_cuda
from scvi.visualization import show_t_sne

torch.set_grad_enabled(False)


def run_benchmarks(gene_dataset, model=VAE, n_epochs=1000, learning_rate=1e-3, use_batches=False, use_cuda=True,
                   show_batch_mixing=True):
    # options:
    # - gene_dataset: a GeneExpressionDataset object
    # call each of the 4 benchmarks:
    # - log-likelihood
    # - imputation
    # - batch mixing
    # - cluster scores
    example_indices = np.random.permutation(len(gene_dataset))
    tt_split = int(0.9 * len(gene_dataset))  # 10%/90% test/train split

    data_loader_train = DataLoader(gene_dataset, batch_size=128, pin_memory=use_cuda,
                                   sampler=SubsetRandomSampler(example_indices[:tt_split]))
    data_loader_test = DataLoader(gene_dataset, batch_size=128, pin_memory=use_cuda,
                                  sampler=SubsetRandomSampler(example_indices[tt_split:]))
    vae = model(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * use_batches, n_labels=gene_dataset.n_labels)

    if vae.using_cuda:
        vae.cuda()
    with torch.set_grad_enabled(True):
        stats = train(vae, data_loader_train, data_loader_test, n_epochs=n_epochs, learning_rate=learning_rate)

    if stats.n_epoch > 0:
        print(stats.history["LL_train"])
        print("Log-likelihood Train:", stats.history["LL_train"][-1])
        print("Log-likelihood Test:", stats.history["LL_test"][-1])

    # - imputation

    imputation_train = imputation(vae, data_loader_train)
    print("Imputation score on train (MAE) is:", imputation_train.item())

    # - batch mixing
    if gene_dataset.n_batches >= 2:
        latent = []
        batch_indices = []
        for tensor_list in data_loader_train:
            if vae.using_cuda:
                tensor_list = to_cuda(tensor_list)
            sample_batch, local_l_mean, local_l_var, batch_index, labels = tensor_list
            sample_batch = sample_batch.type(torch.float32)
            latent += [vae.sample_from_posterior_z(sample_batch, y=labels)]
            batch_indices += [batch_index]
        latent = torch.cat(latent)
        batch_indices = torch.cat(batch_indices)

    if gene_dataset.n_batches == 2:
        print("Entropy batch mixing :", entropy_batch_mixing(latent.cpu().numpy(), batch_indices.cpu().numpy()))
        if show_batch_mixing:
            show_t_sne(latent.cpu().numpy(), np.array([batch[0] for batch in batch_indices.cpu().numpy()]),
                       "Batch mixing t_SNE plot")

    # - differential expression
    #
    if type(gene_dataset) == CortexDataset:
        get_statistics(vae, data_loader_train, M_sampling=1, M_permutation=1)  # 200 - 100000
