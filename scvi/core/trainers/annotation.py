import logging

import numpy as np
import torch
from torch.nn import functional as F

from scvi import _CONSTANTS
from scvi._compat import Literal
from scvi.core.data_loaders import AnnotationDataLoader
from scvi.data._anndata import get_from_registry

from .inference import UnsupervisedTrainer
from .trainer import Trainer

logger = logging.getLogger(__name__)


class ClassifierTrainer(Trainer):
    """
    Class for training a classifier either on the raw data or on top of the latent space of another model.

    Parameters
    ----------
    model
        A model instance from class ``VAE``, ``VAEC``, ``SCANVI``
    gene_dataset
        A gene_dataset instance like ``CortexDataset()``
    train_size
        The train size, a float between 0 and 1 representing proportion of dataset to use for training
        to use Default: ``0.9``.
    test_size
        The test size, a float between 0 and 1 representing proportion of dataset to use for testing
        to use Default: ``None``.
    sampling_model
        Model with z_encoder with which to first transform data.
    sampling_zl
        Transform data with sampling_model z_encoder and l_encoder and concat.
    **kwargs
        Other keywords arguments from the general Trainer class.

    Examples
    --------
    >>> gene_dataset = CortexDataset()
    >>> vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * False,
    ... n_labels=gene_dataset.n_labels)

    >>> classifier = Classifier(vae.n_latent, n_labels=cortex_dataset.n_labels)
    >>> trainer = ClassifierTrainer(classifier, gene_dataset, sampling_model=vae, train_size=0.5)
    >>> trainer.train(n_epochs=20, lr=1e-3)
    >>> trainer.test_set.accuracy()

    """

    def __init__(
        self,
        *args,
        train_size=0.9,
        test_size=None,
        sampling_model=None,
        sampling_zl=False,
        use_cuda=True,
        **kwargs,
    ):
        train_size = float(train_size)
        if train_size > 1.0 or train_size <= 0.0:
            raise ValueError(
                "train_size needs to be greater than 0 and less than or equal to 1"
            )
        self.sampling_model = sampling_model
        self.sampling_zl = sampling_zl
        super().__init__(*args, use_cuda=use_cuda, **kwargs)
        self.train_set, self.test_set, self.validation_set = self.train_test_validation(
            self.model,
            self.adata,
            train_size=train_size,
            test_size=test_size,
            type_class=AnnotationDataLoader,
        )
        self.train_set.to_monitor = ["accuracy"]
        self.test_set.to_monitor = ["accuracy"]
        self.validation_set.to_monitor = ["accuracy"]
        self.train_set.model_zl = sampling_zl
        self.test_set.model_zl = sampling_zl
        self.validation_set.model_zl = sampling_zl

    @property
    def scvi_data_loaders_loop(self):
        return ["train_set"]

    def __setattr__(self, key, value):
        if key in ["train_set", "test_set"]:
            value.sampling_model = self.sampling_model
        super().__setattr__(key, value)

    def loss(self, tensors_labelled):
        x = tensors_labelled[_CONSTANTS.X_KEY]
        b = tensors_labelled[_CONSTANTS.BATCH_KEY]
        labels_train = tensors_labelled[_CONSTANTS.LABELS_KEY]
        if self.sampling_model:
            if hasattr(self.sampling_model, "classify"):
                return F.cross_entropy(
                    self.sampling_model.classify(x, b), labels_train.view(-1)
                )
            else:
                if self.sampling_model.log_variational:
                    x = torch.log(1 + x)
                if self.sampling_zl:
                    x_z = self.sampling_model.z_encoder(x, b)[0]
                    x_l = self.sampling_model.l_encoder(x, b)[0]
                    x = torch.cat((x_z, x_l), dim=-1)
                else:
                    x = self.sampling_model.z_encoder(x, b)[0]
        return F.cross_entropy(self.model(x), labels_train.view(-1))

    def create_scvi_dl(
        self,
        model=None,
        adata=None,
        shuffle=False,
        indices=None,
        type_class=AnnotationDataLoader,
    ):
        return super().create_scvi_dl(model, adata, shuffle, indices, type_class)


class SemiSupervisedTrainer(UnsupervisedTrainer):
    """
    Class for the semi-supervised training of an autoencoder.

    This parent class can be inherited to specify the different training schemes for semi-supervised learning

    Parameters
    ----------
    n_labelled_samples_per_class
        number of labelled samples per class

    """

    def __init__(
        self,
        model,
        adata,
        n_labelled_samples_per_class=50,
        indices_labelled=None,
        indices_unlabelled=None,
        n_epochs_classifier=1,
        lr_classification=5 * 1e-3,
        classification_ratio=50,
        seed=0,
        scheme: Literal["joint", "alternate", "both"] = "both",
        **kwargs,
    ):
        super().__init__(model, adata, **kwargs)
        self.model = model
        self.adata = adata
        self.n_epochs_classifier = n_epochs_classifier
        self.lr_classification = lr_classification
        self.classification_ratio = classification_ratio
        self.scheme = scheme

        if scheme == "joint":
            self.n_epochs_classifier = 0

        if indices_labelled is None and indices_unlabelled is None:
            n_labelled_samples_per_class_array = [
                n_labelled_samples_per_class
            ] * self.adata.uns["_scvi"]["summary_stats"]["n_labels"]
            labels = np.array(
                get_from_registry(self.adata, _CONSTANTS.LABELS_KEY)
            ).ravel()
            np.random.seed(seed=seed)
            permutation_idx = np.random.permutation(len(labels))
            labels = labels[permutation_idx]
            indices = []
            current_nbrs = np.zeros(len(n_labelled_samples_per_class_array))
            for idx, (label) in enumerate(labels):
                label = int(label)
                if current_nbrs[label] < n_labelled_samples_per_class_array[label]:
                    indices.insert(0, idx)
                    current_nbrs[label] += 1
                else:
                    indices.append(idx)
            indices = np.array(indices)
            total_labelled = sum(n_labelled_samples_per_class_array)
            indices_labelled = permutation_idx[indices[:total_labelled]]
            indices_unlabelled = permutation_idx[indices[total_labelled:]]

        self.classifier_trainer = ClassifierTrainer(
            model.classifier,
            self.adata,
            metrics_to_monitor=[],
            silent=True,
            frequency=0,
            sampling_model=self.model,
        )
        self.full_dataset = self.create_scvi_dl(shuffle=True)
        self.labelled_set = self.create_scvi_dl(indices=indices_labelled)
        self.unlabelled_set = self.create_scvi_dl(indices=indices_unlabelled)

        for scdl in [self.labelled_set, self.unlabelled_set]:
            scdl.to_monitor = ["elbo", "reconstruction_error", "accuracy"]
        # allow to track ELBO
        self.unlabelled_set.unlabeled = True
        self.full_dataset.unlabeled = True

    @property
    def scvi_data_loaders_loop(self):

        joint = ["full_dataset", "labelled_set"]
        full = ["full_dataset"]
        if len(self.labelled_set.indices) == 0 or self.scheme == "alternate":
            return full
        else:
            return joint

    def __setattr__(self, key, value):
        if key == "labelled_set":
            self.classifier_trainer.train_set = value
        super().__setattr__(key, value)

    def loss(self, tensors_all, tensors_labelled=None):
        loss = super().loss(tensors_all, feed_labels=False)
        if tensors_labelled is not None:
            sample_batch = tensors_labelled[_CONSTANTS.X_KEY]
            batch_index = tensors_labelled[_CONSTANTS.BATCH_KEY]
            y = tensors_labelled[_CONSTANTS.LABELS_KEY]
            classification_loss = F.cross_entropy(
                self.model.classify(sample_batch, batch_index), y.view(-1)
            )
            loss += classification_loss * self.classification_ratio
        return loss

    def on_epoch_end(self):
        self.model.eval()
        self.classifier_trainer.train(
            self.n_epochs_classifier, lr=self.lr_classification
        )
        self.model.train()
        return super().on_epoch_end()

    def create_scvi_dl(
        self,
        model=None,
        adata=None,
        shuffle=False,
        indices=None,
        type_class=AnnotationDataLoader,
    ):
        return super().create_scvi_dl(model, adata, shuffle, indices, type_class)
