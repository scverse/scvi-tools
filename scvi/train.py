from itertools import cycle

import numpy as np
import torch
from torch.nn import functional as F

from scvi.metrics.stats import Stats
from scvi.utils import to_cuda, enable_grad


@enable_grad()
def train(vae, data_loader_train, data_loader_test, n_epochs=20, lr=0.001, kl=None,
          early_stopping_criterion=(20, 0.01), verbose=True, record_frequency=10,
          reconstruction_ratio=1, classification_ratio=0):
    # Defining the optimizer
    optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, vae.parameters()), lr=lr, eps=0.01)

    # initialize
    (patience, threshold) = (early_stopping_criterion[0], early_stopping_criterion[1])
    current_performances_kl = np.ones((patience))
    current_performances_reconst = np.ones((patience))

    # Getting access to the stats during training
    stats = Stats(verbose, record_frequency, n_epochs=n_epochs)
    stats.callback(vae, data_loader_train, data_loader_test)

    # Training the model
    for epoch in range(n_epochs):
        # initialize kl, reconst
        total_current_kl = 0
        total_current_reconst = 0
        for i_batch, (tensors_train, tensors_test) in enumerate(zip(data_loader_train, cycle(data_loader_test))):
            if vae.use_cuda:
                tensors_train = to_cuda(tensors_train)
                tensors_test = to_cuda(tensors_test)
            sample_batch_train, local_l_mean_train, local_l_var_train, batch_index_train, labels_train = tensors_train
            sample_batch_test, local_l_mean_test, local_l_var_test, batch_index_test, _ = tensors_test
            sample_batch_train = sample_batch_train.type(torch.float32)
            sample_batch_test = sample_batch_test.type(torch.float32)

            if kl is None:
                kl_ponderation = min(1, epoch / 400.)
            else:
                kl_ponderation = kl

            reconst_loss_train, kl_divergence_train = vae(sample_batch_train, local_l_mean_train, local_l_var_train,
                                                          batch_index=batch_index_train, y=labels_train)
            reconst_loss_test, kl_divergence_test = vae(sample_batch_test, local_l_mean_test, local_l_var_test,
                                                        batch_index=batch_index_test, y=None)
            reconst_loss_mean_train = reconstruction_ratio * torch.mean(reconst_loss_train)
            kl_divergence_mean_train = torch.mean(kl_divergence_train)
            reconst_loss_mean_test = reconstruction_ratio * torch.mean(reconst_loss_test)
            kl_divergence_mean_test = torch.mean(kl_divergence_test)

            train_loss = (reconst_loss_mean_test + reconst_loss_mean_train +
                          kl_ponderation * (kl_divergence_mean_train + kl_divergence_mean_test))

            batch_size = sample_batch_train.size(0)
            total_current_kl += kl_divergence_mean_train.item() * batch_size
            total_current_reconst += reconst_loss_mean_train.item() * batch_size

            # Add a classification loss (most semi supervised VAE papers do)
            if hasattr(vae, "classify") and labels_train is not None:
                classification_loss = F.cross_entropy(vae.classify(sample_batch_train), labels_train.view(-1))
                train_loss += classification_loss * classification_ratio

            optimizer.zero_grad()
            train_loss.backward()
            optimizer.step()

        current_performances_kl[:-1] = current_performances_kl[1:]
        current_performances_reconst[:-1] = current_performances_reconst[1:]
        current_performances_kl[-1] = total_current_kl / len(data_loader_train.sampler.indices)
        current_performances_reconst[-1] = total_current_reconst / len(data_loader_train.sampler.indices)

        # Computing the relative improvment of kl and reconstruction loss
        # over the chosen number of epochs
        reconst_relative_improvement = ((current_performances_reconst[0] - current_performances_reconst[-1])
                                        / current_performances_reconst[0])
        kl_relative_improvement = ((current_performances_kl[0] - current_performances_kl[-1])
                                   / current_performances_kl[0])
        # Test whether stopping criterions are met
        if epoch > patience and kl_relative_improvement < threshold and reconst_relative_improvement < threshold:
            # We then stop the iterations
            print("Stopping the training after %d epochs: over the %d last epochs,"
                  " kl divergence improvement was %.4f, reconstruction loss improvement was %.4f"
                  % (epoch, patience, kl_relative_improvement, reconst_relative_improvement))
            return stats
            break

        stats.callback(vae, data_loader_train, data_loader_test)
    return stats


@enable_grad()
def train_classifier(vae, classifier, data_loader_train, data_loader_test, n_epochs=250, lr=0.1,
                     verbose=True, record_frequency=1):
    # Train classifier on the z latent space of a vae
    optimizer = torch.optim.Adam(classifier.parameters(), lr=lr, eps=0.01)

    # Getting access to the stats during training
    stats = Stats(verbose, record_frequency, n_epochs=n_epochs)
    stats.callback(vae, data_loader_train, data_loader_test, classifier)

    for epoch in range(1, n_epochs + 1):
        for i_batch, tensors in enumerate(data_loader_train):
            if vae.use_cuda:
                tensors = to_cuda(tensors)
            sample_batch_train, _, _, _, labels_train = tensors
            sample_batch_train = sample_batch_train.type(torch.float32)
            # Get the latent space
            z = vae.sample_from_posterior_z(sample_batch_train, labels_train)

            optimizer.zero_grad()
            loss = F.cross_entropy(classifier(z), labels_train.view(-1))
            loss.backward()
            optimizer.step()

        stats.callback(vae, data_loader_train, data_loader_test, classifier)

    return stats
