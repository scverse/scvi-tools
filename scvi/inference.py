import sys
from collections import defaultdict

import torch
from torch.nn import functional as F
# model, *data_loaders are attributes of this class
from tqdm import trange

from scvi.dataset.utils import TrainTestDataLoaders, AlternateSemiSupervisedDataLoaders
from scvi.metrics.classification import compute_accuracy
# from scvi.inference.inferencetask import EarlyStopping
from scvi.metrics.early_stopping import EarlyStopping
from scvi.metrics.imputation import imputation
from scvi.metrics.log_likelihood import compute_log_likelihood
from scvi.utils import to_cuda, enable_grad

metrics = {'ll': compute_log_likelihood,
           'accuracy': compute_accuracy,
           'imputation': imputation
           }


class Inference:
    """
    Inference class cares about the general training of a PyTorch model with its statistics
    """
    metrics = []

    def __init__(self, model, gene_dataset, data_loaders=None, metrics_to_monitor=[], use_cuda=True, benchmark=False,
                 record_freq=1, verbose=True, early_stopping_metric=None, save_best_state_metric=None):
        self.model = model
        self.gene_dataset = gene_dataset
        self.data_loaders = data_loaders
        self.benchmark = benchmark
        self.epoch = 0

        self.metrics_to_monitor = metrics_to_monitor
        self.early_stopping_metric = early_stopping_metric
        self.save_best_state_metric = save_best_state_metric

        self.early_stopping = EarlyStopping(benchmark=benchmark)
        self.use_cuda = use_cuda and torch.cuda.is_available()
        if self.use_cuda:
            self.model.cuda()
        self.set_metrics_methods()

        self.verbose = verbose if not benchmark else False
        self.record_freq = record_freq

        self.history = defaultdict(lambda: [])

    def set_metrics_methods(self):
        def metric_method_getter(metric):
            def metric_method(name, **kargs):
                return metrics[metric](self, self.data_loaders[name], use_cuda=self.use_cuda, **kargs)

            return metric_method

        for metric in self.metrics:
            setattr(self, metric, metric_method_getter(metric))

    def compute_metrics(self):
        if self.verbose and (self.epoch == 0 or self.epoch == self.n_epochs or (self.epoch % self.record_freq == 0)):
            print("EPOCH [%d/%d]: " % (self.epoch, self.n_epochs))
            for name in self.data_loaders.to_monitor:
                for metric in self.metrics_to_monitor:
                    result = getattr(self, metric)(name=name)
                    self.history[metric + '_' + name] += [result]
                    print("%s_%s : %.4f" % (metric, name, result))

    @enable_grad()
    def fit(self, n_epochs=20, lr=1e-3):
        optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, self.model.parameters()), lr=lr)
        self.n_epochs = n_epochs
        self.compute_metrics()
        # Training the model

        with trange(n_epochs, desc="training", file=sys.stdout, disable=self.benchmark or self.verbose) as pbar:
            # We have to use tqdm this way so it works in Jupyter notebook.
            # See https://stackoverflow.com/questions/42212810/tqdm-in-jupyter-notebook
            self.on_epoch_begin()
            for epoch in pbar:
                pbar.update(1)
                for tensors_list in self.data_loaders:
                    loss = self.loss(*[to_cuda(tensors, use_cuda=self.use_cuda) for tensors in tensors_list])
                    optimizer.zero_grad()
                    loss.backward()
                    optimizer.step()
                if not self.on_epoch_end():
                    break
        if self.save_best_state_metric is not None:
            self.model.load_state_dict(self.best_state_dict)
            #print(self.best_state_dict)
            print("Now re-computing metric")
            self.compute_metrics()

    def on_epoch_begin(self):
        pass

    def on_epoch_end(self):
        self.epoch += 1
        self.compute_metrics()
        if self.save_best_state_metric is not None:
            if self.early_stopping.update_state(self.history[self.save_best_state_metric][-1]):
                self.best_state_dict = self.model.state_dict()
                self.best_epoch = self.epoch

        continue_training = True
        if self.early_stopping_metric is not None:
            continue_training = self.early_stopping.update(self.history[self.early_stopping_metric][-1])
        return continue_training


class VariationalInference(Inference):
    metrics = ['ll', 'imputation']  # 'de','show_tsne'

    def __init__(self, model, gene_dataset, train_size=0.1, **kwargs):
        super(VariationalInference, self).__init__(model, gene_dataset, **kwargs)
        self.kl = None
        self.data_loaders = TrainTestDataLoaders(self.gene_dataset, train_size=train_size, pin_memory=self.use_cuda)

    def loss(self, tensors):
        sample_batch, local_l_mean, local_l_var, batch_index, _ = tensors
        reconst_loss, kl_divergence = self.model(sample_batch, local_l_mean, local_l_var, batch_index)
        loss = torch.mean(reconst_loss + self.kl_weight * kl_divergence)
        return loss

    def on_epoch_begin(self):
        self.kl_weight = self.kl if self.kl is not None else min(1, self.epoch / self.n_epochs)


class SemiSupervisedVariationalInference(VariationalInference):
    metrics = VariationalInference.metrics + ['accuracy']


class AlternateSemiSupervisedVariationalInference(SemiSupervisedVariationalInference):
    def __init__(self, model, gene_dataset, n_labelled_samples_per_class=50, n_epochs_classifier=1,
                 lr_classification=0.1, **kwargs):
        super(AlternateSemiSupervisedVariationalInference, self).__init__(model, gene_dataset, **kwargs)

        self.n_epochs_classifier = n_epochs_classifier
        self.lr_classification = lr_classification
        self.data_loaders = AlternateSemiSupervisedDataLoaders(gene_dataset, n_labelled_samples_per_class)

        self.classifier_inference = ClassifierInference(
            model.classifier, gene_dataset, metrics_to_monitor=['accuracy'],
            data_loaders=self.data_loaders.classifier_data_loaders(), sampling_model=self.model, verbose=True,
        )

    def on_epoch_end(self):
        self.classifier_inference.fit(self.n_epochs_classifier)
        return super(AlternateSemiSupervisedVariationalInference, self).on_epoch_end()


class JointSemiSupervisedVariationalInference(SemiSupervisedVariationalInference):
    def __init__(self, model, gene_dataset, n_labelled_samples_per_class=50, classification_ratio=100,
                 **kwargs):
        super(JointSemiSupervisedVariationalInference, self).__init__(
            model, gene_dataset, n_labelled_samples_per_class=n_labelled_samples_per_class, **kwargs
        )
        self.classification_ratio = classification_ratio

    def loss(self, tensors_unlabelled, tensors_labelled):
        loss = super(JointSemiSupervisedVariationalInference, self).loss(tensors_unlabelled)
        sample_batch, _, _, _, y = tensors_labelled
        classification_loss = F.cross_entropy(self.model.classify(sample_batch), y.view(-1))
        loss += classification_loss * self.classification_ratio
        return loss


class ClassifierInference(Inference):
    metrics = ['accuracy']

    def __init__(self, *args, sampling_model=None, **kwargs):
        self.sampling_model = sampling_model
        super(ClassifierInference, self).__init__(*args, **kwargs)
        if 'data_loaders' not in kwargs:
            self.data_loaders = TrainTestDataLoaders(self.gene_dataset, train_size=0.1, pin_memory=self.use_cuda)

    def fit(self, *args, **kargs):
        if hasattr(self.model.classifier, "update_parameters"):
            self.model.classifier.update_parameters(self.model, self.data_loaders['train'], use_cuda=self.use_cuda)
        else:
            super(ClassifierInference, self).fit(*args, **kargs)

    def loss(self, tensors_labelled):
        x, _, _, _, labels_train = tensors_labelled
        x = self.sampling_model.sample_from_posterior_z(x) if self.sampling_model is not None else x
        return F.cross_entropy(self.model(x), labels_train.view(-1))
