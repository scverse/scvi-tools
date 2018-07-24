import torch
from torch.nn import functional as F

from scvi.dataset.data_loaders import TrainTestDataLoaders
from scvi.metrics.classification import compute_accuracy
from . import Inference


class ClassifierInference(Inference):
    r"""The ClassifierInference class for training a classifier either on the raw data or on top of the latent
        space of another model (VAE, VAEC, SVAEC).

    Args:
        :model: A model instance from class ``VAE``, ``VAEC``, ``SVAEC``
        :gene_dataset: A gene_dataset instance like ``CortexDataset()``
        :train_size: The train size, either a float between 0 and 1 or and integer for the number of training samples
            to use Default: ``0.8``.
        :\**kwargs: Other keywords arguments from the general Inference class.

    infer_cls = ClassifierInference(cls, cortex_dataset)
    infer_cls.train(n_epochs=1)
    infer_cls.accuracy('train')

    Examples:
        >>> gene_dataset = CortexDataset()
        >>> vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * False,
        ... n_labels=gene_dataset.n_labels)

        >>> cls = Classifier(vae.n_latent, n_labels=cortex_dataset.n_labels)
        >>> infer = ClassifierInference(gene_dataset, sampling_model=vae, train_size=0.5)
        >>> infer.train(n_epochs=20, lr=1e-3)
        >>> infer.accuracy('test')

        >>> cls = Classifier(gene_dataset.nb_genes, n_labels=cortex_dataset.n_labels)
        >>> infer = ClassifierInference(gene_dataset, train_size=0.5)
        >>> infer.train(n_epochs=20, lr=1e-3)
        >>> infer.accuracy('test')

    """
    default_metrics_to_monitor = ['accuracy']

    def __init__(self, *args, sampling_model=None, use_cuda=True, **kwargs):
        self.sampling_model = sampling_model
        super(ClassifierInference, self).__init__(*args, use_cuda=use_cuda, **kwargs)
        if 'data_loaders' not in kwargs:
            self.data_loaders = TrainTestDataLoaders(self.gene_dataset, train_size=0.1)

    def train(self, *args, **kargs):
        if hasattr(self.model, "update_parameters"):
            with torch.no_grad():
                self.model.update_parameters(self.sampling_model, self.data_loaders['train'])
        else:
            super(ClassifierInference, self).train(*args, **kargs)

    def loss(self, tensors_labelled):
        x, _, _, _, labels_train = tensors_labelled
        x = self.sampling_model.sample_from_posterior_z(x) if self.sampling_model is not None else x
        return F.cross_entropy(self.model(x), labels_train.view(-1))

    def accuracy(self, name, verbose=False):
        model, cls = (self.sampling_model, self.model) if hasattr(self, 'sampling_model') else (self.model, None)
        acc = compute_accuracy(model, self.data_loaders[name], classifier=cls)
        if verbose:
            print("Acc for %s is : %.4f" % (name, acc))
        return acc

    accuracy.mode = 'max'
