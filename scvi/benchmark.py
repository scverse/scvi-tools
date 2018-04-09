import torch
from torch.autograd import Variable
from torch.utils.data import DataLoader

from scvi.clustering import entropy_batch_mixing
from scvi.imputation import imputation
from scvi.log_likelihood import log_zinb_positive
from scvi.scvi import VAE
from scvi.train import train

if torch.cuda.is_available():
    dtype = torch.cuda.FloatTensor
else:
    dtype = torch.FloatTensor


def compute_log_likelihood(vae, gene_dataset):
    data_loader_test = DataLoader(gene_dataset, batch_size=gene_dataset.total_size, shuffle=False, num_workers=1)
    for i_batch, (sample_batched, local_l_mean, local_l_var, batch_index) in enumerate(data_loader_test):
        sample_batched = Variable(sample_batched.type(dtype), requires_grad=False)
        px_scale, px_r, px_rate, px_dropout, qz_m, qz_v, ql_m, ql_v = vae(sample_batched)
        log_lkl = torch.mean(-log_zinb_positive(sample_batched, px_rate, torch.exp(px_r), px_dropout)).data[0]
    return log_lkl


def run_benchmarks(gene_dataset, n_epochs=1000, learning_rate=1e-3):
    # options:
    # - gene_dataset: a GeneExpressionDataset object
    # call each of the 4 benchmarks:
    # - log-likelihood
    # - imputation
    # - batch mixing
    # - cluster scores

    torch.backends.cudnn.benchmark = True

    data_loader = DataLoader(gene_dataset, batch_size=128, shuffle=True, num_workers=1)
    vae = VAE(gene_dataset.nb_genes)
    train(vae, data_loader, n_epochs=n_epochs, learning_rate=learning_rate)

    # - log-likelihood

    log_likelihood = compute_log_likelihood(vae, gene_dataset)
    print("Log-likelihood :", log_likelihood)

    # - imputation

    imputation_score = imputation(vae, gene_dataset)
    print("Imputation score (MAE) is:", imputation_score)

    # - batch mixing
    if gene_dataset.n_batches == 2:
        vae(gene_dataset.get_all())  # Just run a forward pass on all the data
        latent = vae.z.data.numpy()
        batches = gene_dataset.get_batches()
        print("Entropy batch mixing :", entropy_batch_mixing(latent, batches))
