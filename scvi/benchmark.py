import numpy as np
import torch
from torch.utils.data import DataLoader
from torch.utils.data.sampler import SubsetRandomSampler

from scvi.clustering import entropy_batch_mixing
from scvi.dataset import CortexDataset
from scvi.differential_expression import get_statistics
from scvi.imputation import imputation
from scvi.log_likelihood import compute_log_likelihood
from scvi.train import train
from scvi.vaec import VAEC, VAE
from scvi.visualization import show_t_sne


def run_benchmarks(gene_dataset, n_epochs=1000, learning_rate=1e-3, use_batches=False, use_cuda=True,
                   show_batch_mixing=True, semi_supervised=False):
    # options:
    # - gene_dataset: a GeneExpressionDataset object
    # call each of the 4 benchmarks:
    # - log-likelihood
    # - imputation
    # - batch mixing
    # - cluster scores
    example_indices = np.random.permutation(len(gene_dataset))
    tt_split = int(0.9 * len(gene_dataset))  # 10%/90% test/train split

    data_loader_train = DataLoader(gene_dataset, batch_size=128,pin_memory=use_cuda,
                                   sampler=SubsetRandomSampler(example_indices[:tt_split]))
    data_loader_test = DataLoader(gene_dataset, batch_size=128, pin_memory=use_cuda,
                                  sampler=SubsetRandomSampler(example_indices[tt_split:]))
    cls = VAEC if semi_supervised else VAE
    vae = cls(gene_dataset.nb_genes, batch=use_batches, n_batch=gene_dataset.n_batches,
              using_cuda=use_cuda, n_labels=gene_dataset.n_labels)

    if vae.using_cuda:
        vae.cuda()
    train(vae, data_loader_train, data_loader_test, n_epochs=n_epochs, learning_rate=learning_rate)

    # - log-likelihood
    vae.eval()  # Test mode - affecting dropout and batchnorm
    log_likelihood_train = compute_log_likelihood(vae, data_loader_train)
    log_likelihood_test = compute_log_likelihood(vae, data_loader_test)
    print("Log-likelihood Train:", log_likelihood_train)
    print("Log-likelihood Test:", log_likelihood_test)

    # - imputation

    imputation_train = imputation(vae, data_loader_train)
    print("Imputation score on train (MAE) is:", imputation_train.item())

    # - batch mixing
    if gene_dataset.n_batches >= 2:
        latent = []
        batch_indices = []
        for sample_batch, local_l_mean, local_l_var, batch_index, labels in data_loader_train:
            sample_batch = sample_batch.type(torch.FloatTensor)
            if vae.using_cuda:
                sample_batch = sample_batch.cuda(async=True)
            latent += [vae.sample_from_posterior_z(sample_batch, y=labels)]  # Just run a forward pass on all the data
            batch_indices += [batch_index]
        latent = torch.cat(latent)
        batch_indices = torch.cat(batch_indices)

    if gene_dataset.n_batches == 2:
        print("Entropy batch mixing :", entropy_batch_mixing(latent.data.cpu().numpy(), batch_indices.numpy()))
        if show_batch_mixing:
            show_t_sne(latent.data.cpu().numpy(), np.array([batch[0] for batch in batch_indices.numpy()]),
                       "Batch mixing t_SNE plot")

    # - differential expression
    #
    if type(gene_dataset) == CortexDataset:
        get_statistics(vae, data_loader_train, M_sampling=1, M_permutation=1)  # 200 - 100000
