from itertools import cycle

import torch
from torch.nn import functional as F

from scvi.metrics.stats import Stats, EarlyStopping
from scvi.utils import to_cuda, enable_grad

# If we're running in Jupyter notebook, then we want to use
# `tqdm_notebook`. Otherwise, we want to use `tqdm`.
try:
    get_ipython
    from tqdm import tqdm_notebook as tqdm
except NameError:
    from tqdm import tqdm


@enable_grad()
def train(vae, data_loader_train, data_loader_test, n_epochs=20, lr=0.001, kl=None, benchmark=False, verbose=False):
    # Defining the optimizer
    optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, vae.parameters()), lr=lr)

    # Getting access to the stats during training
    stats = Stats(n_epochs=n_epochs, benchmark=benchmark, verbose=verbose)
    stats.callback(vae, data_loader_train, data_loader_test)
    early_stopping = EarlyStopping(benchmark=benchmark)

    # Training the model
    for epoch in tqdm(range(n_epochs), desc="training"):
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
def train_semi_supervised_jointly(vae, data_loader_all, data_loader_labelled, data_loader_unlabelled, n_epochs=20,
                                  lr=0.001,
                                  kl=None, benchmark=False, record_freq=1, classification_ratio=100):
    # Defining the optimizer
    optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, vae.parameters()), lr=lr)

    # Getting access to the stats during training
    stats = Stats(n_epochs=n_epochs, benchmark=benchmark, names=["labelled", "unlabelled"], record_freq=record_freq)
    stats.callback(vae, data_loader_labelled, data_loader_unlabelled)
    early_stopping = EarlyStopping(benchmark=benchmark)

    # Training the model
    for epoch in range(n_epochs):
        # initialize kl, reconst
        total_train_loss = 0
        for i_batch, (tensors, tensor_classification) in enumerate(zip(data_loader_all, cycle(data_loader_labelled))):
            with torch.no_grad():
                if vae.use_cuda:
                    tensors = to_cuda(tensors)
            sample_batch, local_l_mean, local_l_var, batch_index, _ = tensors
            sample_batch = sample_batch.type(torch.float32)

            if kl is None:
                kl_ponderation = min(1, epoch / 400.)
            else:
                kl_ponderation = kl

            reconst_loss, kl_divergence = vae(sample_batch, local_l_mean, local_l_var, batch_index=batch_index, y=None)
            loss = torch.mean(reconst_loss + kl_ponderation * kl_divergence)

            # Add a classification loss
            if hasattr(vae, "classify"):
                with torch.no_grad():
                    if vae.use_cuda:
                        tensor_classification = to_cuda(tensor_classification)
                sample_batch_cl, _, _, _, labels_cl = tensor_classification
                sample_batch_cl = sample_batch_cl.type(torch.float32)
                classification_loss = F.cross_entropy(vae.classify(sample_batch_cl), labels_cl.view(-1))
                loss += classification_loss * classification_ratio

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        if not early_stopping.update(total_train_loss):
            break
        stats.callback(vae, data_loader_labelled, data_loader_unlabelled)
    stats.display_time()
    return stats


@enable_grad()
def train_semi_supervised_alternately(vae, data_loader_all, data_loader_labelled, data_loader_unlabelled,
                                      n_epochs=20, lr=0.001, kl=None, benchmark=False, lr_classification=0.1,
                                      record_freq=1, n_epochs_classifier=1):
    # Defining the optimizer
    optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, vae.parameters()), lr=lr)

    # Getting access to the stats during training
    stats = Stats(n_epochs=n_epochs, benchmark=benchmark, names=["labelled", "unlabelled"], record_freq=record_freq)
    stats.callback(vae, data_loader_labelled, data_loader_unlabelled)
    early_stopping = EarlyStopping(benchmark=benchmark)

    # Training the model
    for epoch in range(n_epochs):
        # initialize kl, reconst
        total_train_loss = 0
        for i_batch, tensors in enumerate(data_loader_all):
            with torch.no_grad():
                if vae.use_cuda:
                    tensors = to_cuda(tensors)
            sample_batch, local_l_mean, local_l_var, batch_index, _ = tensors
            sample_batch = sample_batch.type(torch.float32)

            if kl is None:
                kl_ponderation = min(1, epoch / 400.)
            else:
                kl_ponderation = kl

            reconst_loss, kl_divergence = vae(sample_batch, local_l_mean, local_l_var, batch_index=batch_index, y=None)
            loss = torch.mean(reconst_loss + kl_ponderation * kl_divergence)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        # Optimize the classification loss at the end of each epoch
        if hasattr(vae, "classifier"):
            if hasattr(vae.classifier, "update_parameters"):
                vae.classifier.update_parameters(vae, data_loader_labelled)
            else:
                # Minimize the cross entropy over multiple batches of the training data.
                train_classifier(vae, vae.classifier, data_loader_labelled,
                                 n_epochs=n_epochs_classifier, lr=lr_classification, benchmark=True)

        if not early_stopping.update(total_train_loss):
            break
        stats.callback(vae, data_loader_labelled, data_loader_unlabelled)
    stats.display_time()
    return stats


@enable_grad()
def train_classifier(vae, classifier, *data_loaders, n_epochs=250, lr=0.1, benchmark=False):
    # Train classifier on the z latent space of a vae
    optimizer = torch.optim.Adam(classifier.parameters(), lr=lr)

    # Getting access to the stats during training
    stats = Stats(n_epochs=n_epochs, benchmark=benchmark)
    stats.callback(vae, *data_loaders, classifier=classifier)

    for epoch in range(n_epochs):
        for i_batch, tensors in enumerate(data_loaders[0]):
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

        stats.callback(vae, *data_loaders, classifier=classifier)

    return stats
