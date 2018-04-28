import numpy as np
import torch
import matplotlib.pyplot as plt

from torch.utils.data import DataLoader
from torch.utils.data.sampler import SubsetRandomSampler

from scvi.models.svaec import SVAEC
from scvi.modules import Classifier
from scvi.train import train, train_classifier
from scvi.clustering import entropy_batch_mixing
from scvi.dataset import CortexDataset
from scvi.differential_expression import get_statistics
from scvi.imputation import imputation
from scvi.models import VAE
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
    vae = model(gene_dataset.nb_genes, batch=use_batches, n_batch=gene_dataset.n_batches,
                using_cuda=use_cuda, n_labels=gene_dataset.n_labels)

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


# Pipeline to compare different semi supervised models
def run_benchmarks_classification(gene_dataset, n_latent=10, n_epochs=10, n_epochs_classifier=10, learning_rate=1e-3,
                                  use_batches=False, use_cuda=True, verbose=False, record_frequency=1):
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(12, 5))

    lr_classifier = 1e-3 * len(gene_dataset)
    alpha = 0.1

    # Create the dataset
    example_indices = np.random.permutation(len(gene_dataset))
    tt_split = int(0.9 * len(gene_dataset))  # 10%/90% test/train split

    data_loader_train = DataLoader(gene_dataset, batch_size=128, pin_memory=use_cuda,
                                   sampler=SubsetRandomSampler(example_indices[:tt_split]))
    data_loader_test = DataLoader(gene_dataset, batch_size=128, pin_memory=use_cuda,
                                  sampler=SubsetRandomSampler(example_indices[tt_split:]))

    # Now we try out the different models and compare the classification accuracy

    # ========== The M1 model ===========
    print("Trying M1 model")
    vae = VAE(gene_dataset.nb_genes, batch=use_batches, n_batch=gene_dataset.n_batches,
              using_cuda=use_cuda, n_labels=gene_dataset.n_labels)
    # We train it
    if vae.using_cuda:
        vae.cuda()
    with torch.set_grad_enabled(True):
        train(vae, data_loader_train, data_loader_test, n_epochs=n_epochs,
              learning_rate=learning_rate, verbose=verbose, record_frequency=record_frequency)

    # Then we train a classifier on the latent space
    cls = Classifier(n_input=n_latent, n_labels=gene_dataset.n_labels, n_layers=3)
    with torch.set_grad_enabled(True):
        cls_stats = train_classifier(vae, cls, data_loader_train, data_loader_test,
                                     n_epochs=n_epochs_classifier, learning_rate=lr_classifier, verbose=verbose,
                                     record_frequency=record_frequency)
    accuracy_train = cls_stats.history["Accuracy_train"]
    accuracy_test = cls_stats.history["Accuracy_test"]

    axes[0].plot(accuracy_train, label='classifier')
    axes[1].plot(accuracy_test)

    # ========== The M1+M2 model, first encoder frozen ===========
    print("Trying out standard M1+M2 model")
    prior = torch.FloatTensor(
        [(gene_dataset.labels == i).type(torch.float32).mean() for i in range(gene_dataset.n_labels)])

    svaec = SVAEC(gene_dataset.nb_genes, n_labels=gene_dataset.n_labels, y_prior=prior, n_latent=n_latent)

    # Use a pretrained z encoder, freeze its weights
    # svaec.z_encoder.load_state_dict(vae.z_encoder.state_dict())
    for param in svaec.z_encoder.parameters():
        param.requires_grad = False
    with torch.set_grad_enabled(True):
        stats = train(svaec, data_loader_train, data_loader_test, n_epochs=n_epochs, learning_rate=learning_rate,
                      reconstruction_ratio=0, classification_ratio=alpha, verbose=verbose,
                      record_frequency=record_frequency)

    accuracy_train = stats.history["Accuracy_train"]
    accuracy_test = stats.history["Accuracy_test"]
    # We don't train the first z encoder in this procedure
    axes[0].plot(accuracy_train, label='M1+M2 (frozen)')
    axes[1].plot(accuracy_test)

    # ========== The M1+M2 model w/o reconstruction loss ===========
    print("Trying out M1+M2 w/o reconstruction loss")
    prior = torch.FloatTensor(
        [(gene_dataset.labels == i).type(torch.float32).mean() for i in range(gene_dataset.n_labels)])

    vaec = SVAEC(gene_dataset.nb_genes, n_labels=gene_dataset.n_labels, y_prior=prior, n_latent=n_latent)

    # Use a pretrained z encoder
    # vaec.z_encoder.load_state_dict(vae.z_encoder.state_dict())
    with torch.set_grad_enabled(True):
        stats = train(vaec, data_loader_train, data_loader_test, n_epochs=n_epochs, learning_rate=learning_rate,
                      reconstruction_ratio=0, classification_ratio=alpha, verbose=verbose,
                      record_frequency=record_frequency)

    accuracy_train = stats.history["Accuracy_train"]
    accuracy_test = stats.history["Accuracy_test"]

    axes[0].plot(accuracy_train, label='M1+M2 (no recons)')
    axes[1].plot(accuracy_test)

    # ========== The M1+M2 model trained jointly ===========
    print("Trying out M1+M2 optimized jointly")
    prior = torch.FloatTensor(
        [(gene_dataset.labels == i).type(torch.float32).mean() for i in range(gene_dataset.n_labels)])

    vaec = SVAEC(gene_dataset.nb_genes, n_labels=gene_dataset.n_labels, y_prior=prior, n_latent=n_latent)

    # vaec.z_encoder.load_state_dict(vae.z_encoder.state_dict())
    with torch.set_grad_enabled(True):
        stats = train(vaec, data_loader_train, data_loader_test,
                      n_epochs=n_epochs, learning_rate=learning_rate, classification_ratio=alpha, verbose=verbose,
                      record_frequency=record_frequency)
    accuracy_train = stats.history["Accuracy_train"]
    accuracy_test = stats.history["Accuracy_test"]

    axes[0].plot(accuracy_train, label='M1+M2 (train all)')
    axes[1].plot(accuracy_test)

    # ========== Classifier trained on the latent space of M1+M2 ===========
    print("Trying to classify on M1+M2's z1 latent space")
    cls = Classifier(n_input=n_latent, n_labels=gene_dataset.n_labels, n_layers=3)
    with torch.set_grad_enabled(True):
        stats = train_classifier(vaec, cls, data_loader_train, data_loader_test,
                                 n_epochs=n_epochs_classifier, learning_rate=lr_classifier, verbose=verbose,
                                 record_frequency=record_frequency)

    accuracy_train = stats.history["Accuracy_train"]
    accuracy_test = stats.history["Accuracy_test"]

    axes[0].plot(accuracy_train, label='M1+M2+classifier')
    axes[1].plot(accuracy_test)

    # Now plot the results
    axes[0].set_ylim(0, 1)
    axes[0].set_ylabel('accuracy')
    axes[0].set_xlabel('n_epochs')
    axes[0].set_title('acc. train')
    axes[0].legend()
    axes[1].set_xlabel('n_epochs')
    axes[1].set_title('acc. test')
    plt.tight_layout()
    plt.show()
