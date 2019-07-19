from collections import namedtuple

import numpy as np
import logging

from sklearn import neighbors
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC

import torch
from torch.nn import functional as F

from scvi.inference import Posterior
from scvi.inference import Trainer
from scvi.inference.inference import UnsupervisedTrainer
from scvi.inference.posterior import unsupervised_clustering_accuracy

logger = logging.getLogger(__name__)


class AnnotationPosterior(Posterior):
    def __init__(self, *args, model_zl=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.model_zl = model_zl

    def accuracy(self):
        model, cls = (
            (self.sampling_model, self.model)
            if hasattr(self, "sampling_model")
            else (self.model, None)
        )
        acc = compute_accuracy(model, self, classifier=cls, model_zl=self.model_zl)
        logger.debug("Acc: %.4f" % (acc))
        return acc

    accuracy.mode = "max"

    @torch.no_grad()
    def hierarchical_accuracy(self):
        all_y, all_y_pred = self.compute_predictions()
        acc = np.mean(all_y == all_y_pred)

        all_y_groups = np.array([self.model.labels_groups[y] for y in all_y])
        all_y_pred_groups = np.array([self.model.labels_groups[y] for y in all_y_pred])
        h_acc = np.mean(all_y_groups == all_y_pred_groups)

        logger.debug("Hierarchical Acc : %.4f\n" % h_acc)
        return acc

    accuracy.mode = "max"

    @torch.no_grad()
    def compute_predictions(self, soft=False):
        """
        :return: the true labels and the predicted labels
        :rtype: 2-tuple of :py:class:`numpy.int32`
        """
        model, cls = (
            (self.sampling_model, self.model)
            if hasattr(self, "sampling_model")
            else (self.model, None)
        )
        return compute_predictions(
            model, self, classifier=cls, soft=soft, model_zl=self.model_zl
        )

    @torch.no_grad()
    def unsupervised_classification_accuracy(self):
        all_y, all_y_pred = self.compute_predictions()
        uca = unsupervised_clustering_accuracy(all_y, all_y_pred)[0]
        logger.debug("UCA : %.4f" % (uca))
        return uca

    unsupervised_classification_accuracy.mode = "max"

    @torch.no_grad()
    def nn_latentspace(self, posterior):
        data_train, _, labels_train = self.get_latent()
        data_test, _, labels_test = posterior.get_latent()
        nn = KNeighborsClassifier()
        nn.fit(data_train, labels_train)
        score = nn.score(data_test, labels_test)
        return score


class ClassifierTrainer(Trainer):
    r"""The ClassifierInference class for training a classifier either on the raw data or on top of the latent
        space of another model (VAE, VAEC, SCANVI).

    Args:
        :model: A model instance from class ``VAE``, ``VAEC``, ``SCANVI``
        :gene_dataset: A gene_dataset instance like ``CortexDataset()``
        :train_size: The train size, either a float between 0 and 1 or and integer for the number of training samples
            to use Default: ``0.8``.
        :test_size: The test size, either a float between 0 and 1 or and integer for the number of test samples
            to use Default: ``None``.
        :sampling_model: Model with z_encoder with which to first transform data.
        :sampling_zl: Transform data with sampling_model z_encoder and l_encoder and concat.
        :\**kwargs: Other keywords arguments from the general Trainer class.


    Examples:
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
        train_size=0.8,
        test_size=None,
        sampling_model=None,
        sampling_zl=False,
        use_cuda=True,
        **kwargs
    ):
        self.sampling_model = sampling_model
        self.sampling_zl = sampling_zl
        super().__init__(*args, use_cuda=use_cuda, **kwargs)
        self.train_set, self.test_set, self.validation_set = self.train_test_validation(
            self.model,
            self.gene_dataset,
            train_size=train_size,
            test_size=test_size,
            type_class=AnnotationPosterior,
        )
        self.train_set.to_monitor = ["accuracy"]
        self.test_set.to_monitor = ["accuracy"]
        self.validation_set.to_monitor = ["accuracy"]
        self.train_set.model_zl = sampling_zl
        self.test_set.model_zl = sampling_zl
        self.validation_set.model_zl = sampling_zl

    @property
    def posteriors_loop(self):
        return ["train_set"]

    def __setattr__(self, key, value):
        if key in ["train_set", "test_set"]:
            value.sampling_model = self.sampling_model
        super().__setattr__(key, value)

    def loss(self, tensors_labelled):
        x, _, _, _, labels_train = tensors_labelled
        if self.sampling_model:
            if hasattr(self.sampling_model, "classify"):
                return F.cross_entropy(
                    self.sampling_model.classify(x), labels_train.view(-1)
                )
            else:
                if self.sampling_model.log_variational:
                    x = torch.log(1 + x)
                if self.sampling_zl:
                    x_z = self.sampling_model.z_encoder(x)[0]
                    x_l = self.sampling_model.l_encoder(x)[0]
                    x = torch.cat((x_z, x_l), dim=-1)
                else:
                    x = self.sampling_model.z_encoder(x)[0]
        return F.cross_entropy(self.model(x), labels_train.view(-1))

    @torch.no_grad()
    def compute_predictions(self, soft=False):
        """
        :return: the true labels and the predicted labels
        :rtype: 2-tuple of :py:class:`numpy.int32`
        """
        model, cls = (
            (self.sampling_model, self.model)
            if hasattr(self, "sampling_model")
            else (self.model, None)
        )
        full_set = self.create_posterior(type_class=AnnotationPosterior)
        return compute_predictions(
            model, full_set, classifier=cls, soft=soft, model_zl=self.sampling_zl
        )


class SemiSupervisedTrainer(UnsupervisedTrainer):
    r"""The SemiSupervisedTrainer class for the semi-supervised training of an autoencoder.
    This parent class can be inherited to specify the different training schemes for semi-supervised learning
    """

    def __init__(
        self,
        model,
        gene_dataset,
        n_labelled_samples_per_class=50,
        n_epochs_classifier=1,
        lr_classification=5 * 1e-3,
        classification_ratio=50,
        seed=0,
        **kwargs
    ):
        """
        :param n_labelled_samples_per_class: number of labelled samples per class
        """
        super().__init__(model, gene_dataset, **kwargs)
        self.model = model
        self.gene_dataset = gene_dataset

        self.n_epochs_classifier = n_epochs_classifier
        self.lr_classification = lr_classification
        self.classification_ratio = classification_ratio
        n_labelled_samples_per_class_array = [
            n_labelled_samples_per_class
        ] * self.gene_dataset.n_labels
        labels = np.array(self.gene_dataset.labels).ravel()
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
            gene_dataset,
            metrics_to_monitor=[],
            show_progbar=False,
            frequency=0,
            sampling_model=self.model,
        )
        self.full_dataset = self.create_posterior(shuffle=True)
        self.labelled_set = self.create_posterior(indices=indices_labelled)
        self.unlabelled_set = self.create_posterior(indices=indices_unlabelled)

        for posterior in [self.labelled_set, self.unlabelled_set]:
            posterior.to_monitor = ["reconstruction_error", "accuracy"]

    @property
    def posteriors_loop(self):
        return ["full_dataset", "labelled_set"]

    def __setattr__(self, key, value):
        if key == "labelled_set":
            self.classifier_trainer.train_set = value
        super().__setattr__(key, value)

    def loss(self, tensors_all, tensors_labelled):
        loss = super().loss(tensors_all)
        sample_batch, _, _, _, y = tensors_labelled
        classification_loss = F.cross_entropy(
            self.model.classify(sample_batch), y.view(-1)
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

    def create_posterior(
        self,
        model=None,
        gene_dataset=None,
        shuffle=False,
        indices=None,
        type_class=AnnotationPosterior,
    ):
        return super().create_posterior(
            model, gene_dataset, shuffle, indices, type_class
        )


class JointSemiSupervisedTrainer(SemiSupervisedTrainer):
    def __init__(self, model, gene_dataset, **kwargs):
        kwargs.update({"n_epochs_classifier": 0})
        super().__init__(model, gene_dataset, **kwargs)


class AlternateSemiSupervisedTrainer(SemiSupervisedTrainer):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def loss(self, all_tensor):
        return UnsupervisedTrainer.loss(self, all_tensor)

    @property
    def posteriors_loop(self):
        return ["full_dataset"]


@torch.no_grad()
def compute_predictions(
    model, data_loader, classifier=None, soft=False, model_zl=False
):
    all_y_pred = []
    all_y = []

    for i_batch, tensors in enumerate(data_loader):
        sample_batch, _, _, _, labels = tensors
        all_y += [labels.view(-1).cpu()]

        if hasattr(model, "classify"):
            y_pred = model.classify(sample_batch)
        elif classifier is not None:
            # Then we use the specified classifier
            if model is not None:
                if model.log_variational:
                    sample_batch = torch.log(1 + sample_batch)
                if model_zl:
                    sample_z = model.z_encoder(sample_batch)[0]
                    sample_l = model.l_encoder(sample_batch)[0]
                    sample_batch = torch.cat((sample_z, sample_l), dim=-1)
                else:
                    sample_batch, _, _ = model.z_encoder(sample_batch)
            y_pred = classifier(sample_batch)
        else:  # The model is the raw classifier
            y_pred = model(sample_batch)

        if not soft:
            y_pred = y_pred.argmax(dim=-1)

        all_y_pred += [y_pred.cpu()]

    all_y_pred = np.array(torch.cat(all_y_pred))
    all_y = np.array(torch.cat(all_y))

    return all_y, all_y_pred


@torch.no_grad()
def compute_accuracy(vae, data_loader, classifier=None, model_zl=False):
    all_y, all_y_pred = compute_predictions(
        vae, data_loader, classifier=classifier, model_zl=model_zl
    )
    return np.mean(all_y == all_y_pred)


Accuracy = namedtuple(
    "Accuracy", ["unweighted", "weighted", "worst", "accuracy_classes"]
)


@torch.no_grad()
def compute_accuracy_tuple(y, y_pred):
    y = y.ravel()
    n_labels = len(np.unique(y))
    classes_probabilities = []
    accuracy_classes = []
    for cl in range(n_labels):
        idx = y == cl
        classes_probabilities += [np.mean(idx)]
        accuracy_classes += [
            np.mean((y[idx] == y_pred[idx])) if classes_probabilities[-1] else 0
        ]
        # This is also referred to as the "recall": p = n_true_positive / (n_false_negative + n_true_positive)
        # ( We could also compute the "precision": p = n_true_positive / (n_false_positive + n_true_positive) )
        accuracy_named_tuple = Accuracy(
            unweighted=np.dot(accuracy_classes, classes_probabilities),
            weighted=np.mean(accuracy_classes),
            worst=np.min(accuracy_classes),
            accuracy_classes=accuracy_classes,
        )
    return accuracy_named_tuple


@torch.no_grad()
def compute_accuracy_nn(data_train, labels_train, data_test, labels_test, k=5):
    clf = neighbors.KNeighborsClassifier(k, weights="distance")
    return compute_accuracy_classifier(
        clf, data_train, labels_train, data_test, labels_test
    )


@torch.no_grad()
def compute_accuracy_classifier(clf, data_train, labels_train, data_test, labels_test):
    clf.fit(data_train, labels_train)
    # Predicting the labels
    y_pred_test = clf.predict(data_test)
    y_pred_train = clf.predict(data_train)

    return (
        (
            compute_accuracy_tuple(labels_train, y_pred_train),
            compute_accuracy_tuple(labels_test, y_pred_test),
        ),
        y_pred_test,
    )


@torch.no_grad()
def compute_accuracy_svc(
    data_train,
    labels_train,
    data_test,
    labels_test,
    param_grid=None,
    verbose=0,
    max_iter=-1,
):
    if param_grid is None:
        param_grid = [
            {"C": [1, 10, 100, 1000], "kernel": ["linear"]},
            {"C": [1, 10, 100, 1000], "gamma": [0.001, 0.0001], "kernel": ["rbf"]},
        ]
    svc = SVC(max_iter=max_iter)
    clf = GridSearchCV(svc, param_grid, verbose=verbose)
    return compute_accuracy_classifier(
        clf, data_train, labels_train, data_test, labels_test
    )


@torch.no_grad()
def compute_accuracy_rf(
    data_train, labels_train, data_test, labels_test, param_grid=None, verbose=0
):
    if param_grid is None:
        param_grid = {"max_depth": np.arange(3, 10), "n_estimators": [10, 50, 100, 200]}
    rf = RandomForestClassifier(max_depth=2, random_state=0)
    clf = GridSearchCV(rf, param_grid, verbose=verbose)
    return compute_accuracy_classifier(
        clf, data_train, labels_train, data_test, labels_test
    )
