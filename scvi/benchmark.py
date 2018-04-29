import matplotlib.pyplot as plt
import numpy as np
import torch
from torch.utils.data import DataLoader
from torch.utils.data.sampler import SubsetRandomSampler

from scvi.dataset import CortexDataset
from scvi.metrics.clustering import entropy_batch_mixing, get_latent
from scvi.metrics.differential_expression import get_statistics
from scvi.metrics.imputation import imputation
from scvi.metrics.visualization import show_t_sne
from scvi.models import VAE, SVAEC
from scvi.models.modules import Classifier
from scvi.train import train, train_classifier


def run_benchmarks(gene_dataset, model=VAE, n_epochs=1000, lr=1e-3, use_batches=False, use_cuda=True,
                   show_batch_mixing=True):
    # options:
    # - gene_dataset: a GeneExpressionDataset object
    # call each of the 4 benchmarks:
    # - log-likelihood
    # - imputation
    # - batch mixing
    # - cluster scores
    example_indices = np.random.permutation(len(gene_dataset))
    tt_split = int(0.5 * len(gene_dataset))  # 50%/50% train/test split

    data_loader_train = DataLoader(gene_dataset, batch_size=128, pin_memory=use_cuda,
                                   sampler=SubsetRandomSampler(example_indices[:tt_split]))
    data_loader_test = DataLoader(gene_dataset, batch_size=128, pin_memory=use_cuda,
                                  sampler=SubsetRandomSampler(example_indices[tt_split:]))
    vae = model(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * use_batches, n_labels=gene_dataset.n_labels,
                use_cuda=use_cuda)

    stats = train(vae, data_loader_train, data_loader_test, n_epochs=n_epochs, lr=lr)

    # - log-likelihood
    print("Log-likelihood Train:", stats.history["LL_train"][-1])
    print("Log-likelihood Test:", stats.history["LL_test"][-1])

    # - imputation
    imputation_train = imputation(vae, data_loader_train)
    print("Imputation score on train (MAE) is:", imputation_train.item())

    # - batch mixing
    if gene_dataset.n_batches == 2:
        latent, batch_indices = get_latent(vae, data_loader_train)
        print("Entropy batch mixing :", entropy_batch_mixing(latent.cpu().numpy(), batch_indices.cpu().numpy()))
        if show_batch_mixing:
            show_t_sne(latent.cpu().numpy(), np.array([batch[0] for batch in batch_indices.cpu().numpy()]),
                       "Batch mixing t_SNE plot")

    # - differential expression
    if type(gene_dataset) == CortexDataset:
        get_statistics(vae, data_loader_train, M_sampling=1, M_permutation=1)  # 200 - 100000


# Pipeline to compare different semi supervised models
def run_benchmarks_classification(gene_dataset, n_latent=10, n_epochs=10, n_epochs_classifier=10, lr=1e-2,
                                  use_batches=False, use_cuda=True, verbose=False, record_frequency=1):
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(12, 5))

    alpha = 100  # in Kingma, 0.1 * len(gene_dataset), but pb when : len(gene_dataset) >> 1

    # Create the dataset
    example_indices = np.random.permutation(len(gene_dataset))
    tt_split = int(0.5 * len(gene_dataset))  # 50%/50% train/test split

    data_loader_train = DataLoader(gene_dataset, batch_size=128, pin_memory=use_cuda,
                                   sampler=SubsetRandomSampler(example_indices[:tt_split]))
    data_loader_test = DataLoader(gene_dataset, batch_size=128, pin_memory=use_cuda,
                                  sampler=SubsetRandomSampler(example_indices[tt_split:]))

    # Now we try out the different models and compare the classification accuracy

    # ========== The M1 model ===========
    print("Trying M1 model")
    vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * use_batches, use_cuda=use_cuda,
              n_labels=gene_dataset.n_labels)
    train(vae, data_loader_train, data_loader_test, n_epochs=n_epochs,
          lr=lr, verbose=verbose, record_frequency=record_frequency)

    # Then we train a classifier on the latent space
    cls = Classifier(n_input=n_latent, n_labels=gene_dataset.n_labels, n_layers=3, use_cuda=use_cuda)
    for param in vae.z_encoder.parameters():
        param.requires_grad = False
    cls_stats = train_classifier(vae, cls, data_loader_train, data_loader_test, n_epochs=n_epochs_classifier,
                                 lr=lr, verbose=verbose, record_frequency=record_frequency)

    axes[0].plot(cls_stats.history["Accuracy_train"], label='classifier')
    axes[1].plot(cls_stats.history["Accuracy_test"])

    # ========== The M1+M2 model, first encoder frozen ===========
    print("Trying out standard M1+M2 model")
    prior = torch.FloatTensor(
        [(gene_dataset.labels == i).type(torch.float32).mean() for i in range(gene_dataset.n_labels)])

    svaec = SVAEC(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * use_batches, n_labels=gene_dataset.n_labels,
                  y_prior=prior, n_latent=n_latent, use_cuda=use_cuda)

    # Use a pretrained z encoder, freeze its weights
    svaec.z_encoder.load_state_dict(vae.z_encoder.state_dict())
    for param in svaec.z_encoder.parameters():
        param.requires_grad = False
    stats = train(svaec, data_loader_train, data_loader_test, n_epochs=n_epochs, lr=lr,
                  reconstruction_ratio=0, classification_ratio=alpha, verbose=verbose,
                  record_frequency=record_frequency)

    # We don't train the first z encoder in this procedure
    axes[0].plot(stats.history["Accuracy_train"], label='M1+M2 (frozen)')
    axes[1].plot(stats.history["Accuracy_test"])

    # ========== The M1+M2 model w/o reconstruction loss ===========
    print("Trying out M1+M2 w/o reconstruction loss")
    prior = torch.FloatTensor(
        [(gene_dataset.labels == i).type(torch.float32).mean() for i in range(gene_dataset.n_labels)])

    vaec = SVAEC(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * use_batches, n_labels=gene_dataset.n_labels,
                 y_prior=prior, use_cuda=use_cuda,
                 n_latent=n_latent)

    # Use a pretrained z encoder
    # vaec.z_encoder.load_state_dict(vae.z_encoder.state_dict())
    stats = train(vaec, data_loader_train, data_loader_test, n_epochs=n_epochs, lr=lr,
                  reconstruction_ratio=0, classification_ratio=alpha, verbose=verbose,
                  record_frequency=record_frequency)

    axes[0].plot(stats.history["Accuracy_train"], label='M1+M2 (no recons)')
    axes[1].plot(stats.history["Accuracy_test"])

    # ========== The M1+M2 model trained jointly ===========
    print("Trying out M1+M2 optimized jointly")
    prior = torch.FloatTensor(
        [(gene_dataset.labels == i).type(torch.float32).mean() for i in range(gene_dataset.n_labels)])

    vaec = SVAEC(gene_dataset.nb_genes, n_labels=gene_dataset.n_labels, y_prior=prior, n_latent=n_latent,
                 use_cuda=use_cuda)

    # vaec.z_encoder.load_state_dict(vae.z_encoder.state_dict())
    stats = train(vaec, data_loader_train, data_loader_test, n_epochs=n_epochs, lr=lr,
                  classification_ratio=alpha, verbose=verbose, record_frequency=record_frequency)

    axes[0].plot(stats.history["Accuracy_train"], label='M1+M2 (train all)')
    axes[1].plot(stats.history["Accuracy_test"])

    # ========== Classifier trained on the latent space of M1+M2 ===========
    print("Trying to classify on M1+M2's z1 latent space")
    cls = Classifier(n_input=n_latent, n_labels=gene_dataset.n_labels, n_layers=3, use_cuda=use_cuda)

    stats = train_classifier(vaec, cls, data_loader_train, data_loader_test, n_epochs=n_epochs_classifier,
                             lr=lr, verbose=verbose, record_frequency=record_frequency)  # alpha

    axes[0].plot(stats.history["Accuracy_train"], label='M1+M2+classifier')
    axes[1].plot(stats.history["Accuracy_test"])

    # Now plot the results
    axes[0].set_ylim(0, 1)
    axes[0].set_ylabel('accuracy')
    axes[0].set_xlabel('n_epochs')
    axes[0].set_title('acc. train')
    axes[0].legend()
    axes[1].set_xlabel('n_epochs')
    axes[1].set_title('acc. test')

    plt.tight_layout()
    plt.savefig("result_classification.png")
