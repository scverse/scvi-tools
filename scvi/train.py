import sys
from itertools import cycle

import torch
from torch.nn import functional as F
from tqdm import trange

from scvi.metrics.stats import Stats, EarlyStopping
from scvi.utils import to_cuda, enable_grad


class Trainer:
    def on_epoch_begin(self):
        self.early_stopping_loss = 0
        self.kl_weight = self.kl if self.kl is not None else min(1, self.epoch / 400.)

    def loop(self, model, *tensors_list, classifier=None):
        reconst_loss, kl_divergence = model(*tensors_list[0])
        loss = torch.mean(reconst_loss + self.kl_weight * kl_divergence)
        self.early_stopping_loss += loss.item()
        return loss

    def on_epoch_end(self, model, *data_loaders, classifier=None):
        self.stats.callback(model, *data_loaders, classifier=classifier)
        return self.early_stopping.update(self.early_stopping_loss)

    @enable_grad()
    def train(self, model, *data_loaders, n_epochs=20, lr=0.001, kl=None, benchmark=False, verbose=False,
              record_freq=30, use_cuda=True, classifier=None):
        self.kl = kl
        self.use_cuda = use_cuda
        self.epoch = 0
        trainable_parameters = list(model.parameters())
        if classifier is not None:
            trainable_parameters = list(classifier.parameters())
        self.optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, trainable_parameters), lr=lr)
        self.stats = Stats(n_epochs=n_epochs, benchmark=benchmark, verbose=verbose, record_freq=record_freq,
                           use_cuda=use_cuda)
        self.stats.callback(model, *data_loaders, classifier=classifier, increment=False)
        self.early_stopping = EarlyStopping(benchmark=benchmark)

        # Training the model
        with trange(n_epochs, desc="training", file=sys.stdout, disable=verbose or benchmark) as pbar:
            # We have to use tqdm this way so it works in Jupyter notebook.
            # See https://stackoverflow.com/questions/42212810/tqdm-in-jupyter-notebook
            for epoch in pbar:
                pbar.update(1)
                self.on_epoch_begin()
                for tensors_list in zip(data_loaders[0], *[cycle(data_loader) for data_loader in data_loaders[1:]]):
                    updated_tensors_list = []
                    for tensors in tensors_list:
                        if self.use_cuda:
                            tensors = to_cuda(tensors)
                        tensors = (tensors[0].type(torch.float32),) + tensors[1:]
                        updated_tensors_list += [tensors]
                    loss = self.loop(model, *updated_tensors_list, classifier=classifier)
                    self.optimizer.zero_grad()
                    loss.backward()
                    self.optimizer.step()
                self.epoch += 1
                if not self.on_epoch_end(model, *data_loaders):
                    break
            self.stats.callback(model, *data_loaders, classifier=classifier, increment=False)
            self.stats.display_time()
            self.stats.set_best_params(model)
            return self.stats


class JointSemiSupervisedTrainer(Trainer):
    def __init__(self, classification_ratio=100):
        super(JointSemiSupervisedTrainer, self).__init__()
        self.classification_ratio = classification_ratio

    def loop(self, model, *tensors_list, classifier=None):
        all, labelled = tensors_list[:2]
        sample_batch, local_l_mean, local_l_var, batch_index, _ = all
        reconst_loss, kl_divergence = model(sample_batch, local_l_mean, local_l_var, batch_index=batch_index)
        loss = torch.mean(reconst_loss + self.kl_weight * kl_divergence)

        sample_batch, _, _, _, y = labelled
        classification_loss = F.cross_entropy(model.classify(sample_batch), y.view(-1))
        loss += classification_loss * self.classification_ratio
        return loss


class AlternateSemiSupervisedTrainer(Trainer):
    def __init__(self, n_epochs_classifier=1, lr_classification=0.1):
        super(AlternateSemiSupervisedTrainer, self).__init__()
        self.n_epochs_classifier = n_epochs_classifier
        self.lr_classification = lr_classification
        self.classifier_trainer = ClassifierTrainer()

    def loop(self, model, *tensors_list, classifier=None):
        reconst_loss, kl_divergence = model(*(tensors_list[0][:-1]))  # Not giving the label
        loss = torch.mean(reconst_loss + self.kl_weight * kl_divergence)
        self.early_stopping_loss += loss.item()
        return loss

    def on_epoch_end(self, model, *data_loaders):
        data_loader_labelled = data_loaders[-1]
        if hasattr(model.classifier, "update_parameters"):
            model.classifier.update_parameters(model, data_loader_labelled, use_cuda=self.use_cuda)
        else:
            self.classifier_trainer.train(
                model, data_loader_labelled, benchmark=True, verbose=False, n_epochs=self.n_epochs_classifier,
                classifier=model.classifier, use_cuda=self.use_cuda)
        return super(AlternateSemiSupervisedTrainer, self).on_epoch_end(model, *data_loaders)


class ClassifierTrainer(Trainer):
    def __init__(self, use_model=False):
        self.use_model = use_model

    def loop(self, model, *tensors_list, classifier=None):
        x, _, _, _, labels_train = tensors_list[0]
        x = model.sample_from_posterior_z(x) if self.use_model else x
        return F.cross_entropy(classifier(x), labels_train.view(-1))
