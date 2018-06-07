from itertools import cycle

import torch
from torch.nn import functional as F

from scvi.metrics.stats import Stats, EarlyStopping
from scvi.utils import to_cuda, enable_grad


@enable_grad()
def train(vae, data_loader_train, data_loader_test, n_epochs=20, lr=0.001, kl=None, benchmark=False):
    # Defining the optimizer
    optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, vae.parameters()), lr=lr)

    # Getting access to the stats during training
    stats = Stats(n_epochs=n_epochs, benchmark=benchmark)
    stats.callback(vae, data_loader_train, data_loader_test)
    early_stopping = EarlyStopping(benchmark=benchmark)

    # Training the model
    for epoch in range(n_epochs):
        total_train_loss = 0
        for i_batch, (tensors_train) in enumerate(data_loader_train):
            if vae.use_cuda:
                tensors_train = to_cuda(tensors_train)
            sample_batch_train, local_l_mean_train, local_l_var_train, batch_index_train, labels_train = tensors_train
            sample_batch_train = sample_batch_train.type(torch.float32)

            if kl is None:
                kl_ponderation = min(1, epoch / 400.)
            else:
                kl_ponderation = kl

            reconst_loss_train, kl_divergence_train = vae(sample_batch_train, local_l_mean_train, local_l_var_train,
                                                          batch_index=batch_index_train, y=labels_train)
            train_loss = torch.mean(reconst_loss_train + kl_ponderation * kl_divergence_train)

            batch_size = sample_batch_train.size(0)
            total_train_loss += train_loss.item() * batch_size
            optimizer.zero_grad()
            train_loss.backward()
            optimizer.step()

        if not early_stopping.update(total_train_loss):
            break
        stats.callback(vae, data_loader_train, data_loader_test)
    stats.display_time()
    stats.set_best_params(vae)
    return stats


@enable_grad()
def train_semi_supervised(vae, data_loader_train, data_loader_test, n_epochs=20, lr=0.001, kl=None, benchmark=False,
                          classification_ratio=0):
    # Defining the optimizer
    optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, vae.parameters()), lr=lr)

    # Getting access to the stats during training
    stats = Stats(n_epochs=n_epochs, benchmark=benchmark)
    stats.callback(vae, data_loader_train, data_loader_test)
    early_stopping = EarlyStopping(benchmark=benchmark)

    # Training the model
    for epoch in range(n_epochs):
        # initialize kl, reconst
        total_train_loss = 0
        for i_batch, (tensors_train, tensors_test) in enumerate(zip(data_loader_train, cycle(data_loader_test))):
            with torch.no_grad():
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

            train_loss = torch.mean(reconst_loss_train + kl_ponderation * kl_divergence_train)
            test_loss = torch.mean(reconst_loss_test + kl_ponderation * kl_divergence_test)
            train_test_loss = train_loss + test_loss

            # Add a classification loss (most semi supervised VAE papers do)
            if hasattr(vae, "classify") and labels_train is not None:
                classification_loss = F.cross_entropy(vae.classify(sample_batch_train), labels_train.view(-1))
                train_test_loss += classification_loss * classification_ratio

            batch_size = sample_batch_train.size(0)
            total_train_loss += train_loss.item() * batch_size
            optimizer.zero_grad()
            train_test_loss.backward()
            optimizer.step()

        if not early_stopping.update(total_train_loss):
            break
        stats.callback(vae, data_loader_train, data_loader_test)
    stats.display_time()
    return stats


@enable_grad()
def train_classifier(vae, classifier, data_loader_train, data_loader_test, n_epochs=250, lr=0.1):
    # Train classifier on the z latent space of a vae
    optimizer = torch.optim.Adam(classifier.parameters(), lr=lr, eps=0.01)

    # Getting access to the stats during training
    stats = Stats(n_epochs=n_epochs)
    stats.callback(vae, data_loader_train, data_loader_test, classifier=classifier)

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

        stats.callback(vae, data_loader_train, data_loader_test, classifier=classifier)

    return stats
