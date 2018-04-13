import torch
from torch.utils.data import DataLoader

from scvi.clustering import entropy_batch_mixing
from scvi.imputation import imputation
from scvi.scvi import VAE
from scvi.train import train


def run_benchmarks(gene_dataset_train, gene_dataset_test, n_epochs=1000, learning_rate=1e-3, use_batches=False):
    # options:
    # - gene_dataset: a GeneExpressionDataset object
    # call each of the 4 benchmarks:
    # - log-likelihood
    # - imputation
    # - batch mixing
    # - cluster scores

    torch.backends.cudnn.benchmark = True

    data_loader_train = DataLoader(gene_dataset_train, batch_size=128, shuffle=True, num_workers=1)
    data_loader_test = DataLoader(gene_dataset_test, batch_size=128, shuffle=True, num_workers=1)
    vae = VAE(gene_dataset_train.nb_genes, batch=use_batches, n_batch=gene_dataset_train.n_batches)
    if torch.cuda.is_available():
        vae.cuda()
    train(vae, data_loader_train, data_loader_test, n_epochs=n_epochs, learning_rate=learning_rate)

    # - log-likelihood
    vae.eval()  # Test mode - affecting dropout and batchnorm
    log_likelihood_train = vae.compute_log_likelihood(data_loader_train)
    log_likelihood_test = vae.compute_log_likelihood(data_loader_test)
    print("Log-likelihood Train:", log_likelihood_train)
    print("Log-likelihood Test:", log_likelihood_test)

    # - imputation

    imputation_train = imputation(vae, data_loader_train)
    print("Imputation score on train (MAE) is:", imputation_train)

    # - batch mixing
    if gene_dataset_train.n_batches == 2:
        latent = []
        batch_indices = []
        for sample_batch, local_l_mean, local_l_var, batch_index in data_loader_train:
            latent += [vae.sample_from_posterior(sample_batch)]  # Just run a forward pass on all the data
            batch_indices += [batch_index]
        latent = torch.cat(latent)
        batch_indices = torch.cat(batch_indices)
        print("Entropy batch mixing :", entropy_batch_mixing(latent.data.numpy(), batch_indices.numpy()))
