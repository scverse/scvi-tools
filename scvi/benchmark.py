import matplotlib.pyplot as plt
import numpy as np
import torch
from torch.utils.data import DataLoader
from torch.utils.data.sampler import SubsetRandomSampler

from scvi.dataset import CortexDataset
from scvi.dataset.utils import get_data_loaders, get_raw_data
from scvi.metrics.adapt_encoder import adapt_encoder
from scvi.metrics.classification import compute_accuracy_rf, compute_accuracy_svc
from scvi.metrics.clustering import entropy_batch_mixing, get_latent
from scvi.metrics.differential_expression import get_statistics
from scvi.metrics.imputation import imputation
from scvi.metrics.visualization import show_t_sne
from scvi.models import VAE, SVAEC
from scvi.models.classifier import Classifier
from scvi.train import train, train_classifier, train_semi_supervised_jointly, train_semi_supervised_alternately


def run_benchmarks(gene_dataset, model=VAE, n_epochs=1000, lr=1e-3, use_batches=False, use_cuda=True,
                   show_batch_mixing=True, benchmark=False, tt_split=0.9):
    # options:
    # - gene_dataset: a GeneExpressionDataset object
    # call each of the 4 benchmarks:
    # - log-likelihood
    # - imputation
    # - batch mixing
    # - cluster scores
    example_indices = np.random.permutation(len(gene_dataset))
    tt_split = int(tt_split * len(gene_dataset))  # 90%/10% train/test split

    data_loader_train = DataLoader(gene_dataset, batch_size=128, pin_memory=use_cuda,
                                   sampler=SubsetRandomSampler(example_indices[:tt_split]),
                                   collate_fn=gene_dataset.collate_fn)
    data_loader_test = DataLoader(gene_dataset, batch_size=128, pin_memory=use_cuda,
                                  sampler=SubsetRandomSampler(example_indices[tt_split:]),
                                  collate_fn=gene_dataset.collate_fn)
    vae = model(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * use_batches, n_labels=gene_dataset.n_labels,
                use_cuda=use_cuda)
    stats = train(vae, data_loader_train, data_loader_test, n_epochs=n_epochs, lr=lr, benchmark=benchmark)

    if isinstance(vae, VAE):
        best_ll = adapt_encoder(vae, data_loader_test, n_path=1, n_epochs=1, record_freq=1)
        print("Best ll was :", best_ll)

    # - log-likelihood
    print("Log-likelihood Train:", stats.history["LL_train"][stats.best_index])
    print("Log-likelihood Test:", stats.history["LL_test"][stats.best_index])

    # - imputation
    imputation_test = imputation(vae, data_loader_test)
    print("Imputation score on test (MAE) is:", imputation_test.item())

    # - batch mixing
    if gene_dataset.n_batches == 2:
        latent, batch_indices, labels = get_latent(vae, data_loader_train)
        print("Entropy batch mixing :", entropy_batch_mixing(latent, batch_indices))
        if show_batch_mixing:
            show_t_sne(latent, np.array([batch[0] for batch in batch_indices]))

    # - differential expression
    if type(gene_dataset) == CortexDataset:
        get_statistics(vae, data_loader_train, M_sampling=1, M_permutation=1)  # 200 - 100000


# Pipeline to compare different semi supervised models
def run_benchmarks_classification(gene_dataset, n_latent=10, n_epochs=10, n_epochs_classifier=10, lr=1e-2,
                                  use_batches=False, use_cuda=True, tt_split=0.9):
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(12, 5))
    alpha = 100  # in Kingma, 0.1 * len(gene_dataset), but pb when : len(gene_dataset) >> 1

    data_loader_all, data_loader_labelled, data_loader_unlabelled = get_data_loaders(gene_dataset, 10,
                                                                                     batch_size=128,
                                                                                     pin_memory=use_cuda)
    # Now we try out the different models and compare the classification accuracy

    (data_train, labels_train), (data_test, labels_test) = get_raw_data(data_loader_labelled, data_loader_unlabelled)
    accuracy_train_svc, accuracy_test_svc = compute_accuracy_svc(data_train, labels_train, data_test, labels_test,
                                                                 unit_test=True)
    accuracy_train_rf, accuracy_test_rf = compute_accuracy_rf(data_train, labels_train, data_test, labels_test,
                                                              unit_test=True)

    # ========== The M1 model ===========
    print("Trying M1 model")
    vae = VAE(gene_dataset.nb_genes, n_latent=n_latent,
              n_batch=gene_dataset.n_batches * use_batches, use_cuda=use_cuda,
              n_labels=gene_dataset.n_labels)
    train_semi_supervised_jointly(vae, data_loader_all, data_loader_labelled, data_loader_unlabelled,
                                  n_epochs=n_epochs, lr=lr)

    # Then we train a classifier on the latent space
    cls = Classifier(n_input=n_latent, n_labels=gene_dataset.n_labels, n_layers=3, use_cuda=use_cuda)
    for param in vae.z_encoder.parameters():
        param.requires_grad = False
    cls_stats = train_classifier(vae, cls, data_loader_labelled, data_loader_unlabelled, n_epochs=n_epochs_classifier,
                                 lr=lr)

    axes[0].plot(cls_stats.history["Accuracy_train"], label='classifier')
    axes[1].plot(cls_stats.history["Accuracy_test"])

    # ========== The M1+M2 model, first encoder frozen ===========
    print("Trying out standard M1+M2 model")
    prior = torch.FloatTensor([(gene_dataset.labels == i).mean() for i in range(gene_dataset.n_labels)])

    svaec = SVAEC(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * use_batches, n_labels=gene_dataset.n_labels,
                  y_prior=prior, n_latent=n_latent, use_cuda=use_cuda)

    # Use a pretrained z encoder, freeze its weights
    svaec.z_encoder.load_state_dict(vae.z_encoder.state_dict())
    for param in svaec.z_encoder.parameters():
        param.requires_grad = False
    stats = train_semi_supervised_jointly(svaec, data_loader_all, data_loader_labelled, data_loader_unlabelled,
                                          n_epochs=n_epochs, lr=lr, classification_ratio=alpha)

    # We don't train the first z encoder in this procedure
    axes[0].plot(stats.history["Accuracy_train"], label='M1+M2 (frozen)')
    axes[1].plot(stats.history["Accuracy_test"])

    # ========== The M1+M2 model trained jointly ===========
    print("Trying out M1+M2 optimized jointly")
    svaec = SVAEC(gene_dataset.nb_genes, n_labels=gene_dataset.n_labels, y_prior=prior, n_latent=n_latent,
                  use_cuda=use_cuda, logreg_classifier=False)

    stats = train_semi_supervised_jointly(svaec, data_loader_all, data_loader_labelled, data_loader_unlabelled,
                                          n_epochs=n_epochs, lr=lr, classification_ratio=100, record_freq=10)

    svaec = SVAEC(gene_dataset.nb_genes, n_labels=gene_dataset.n_labels, y_prior=prior, n_latent=n_latent,
                  use_cuda=use_cuda, logreg_classifier=True)

    stats = train_semi_supervised_alternately(svaec, data_loader_all, data_loader_labelled, data_loader_unlabelled,
                                              n_epochs=n_epochs, lr=lr, record_freq=10, lr_classification=0.05)

    axes[0].plot(stats.history["Accuracy_train"], label='M1+M2 (train all)')
    axes[1].plot(stats.history["Accuracy_test"])

    # ========== Classifier trained on the latent space of M1+M2 ===========
    print("Trying to classify on M1+M2's z1 latent space")
    cls = Classifier(n_input=n_latent, n_labels=gene_dataset.n_labels, n_layers=3, use_cuda=use_cuda)

    stats = train_classifier(svaec, cls, data_loader_labelled, data_loader_unlabelled, n_epochs=n_epochs_classifier,
                             lr=lr)

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
