from collections import namedtuple

import numpy as np
import torch
from sklearn import neighbors
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from torch.nn import functional as F

from scvi.inference import Trainer, Posterior
from scvi.inference.posterior import unsupervised_clustering_accuracy
from scvi.inference.trainers import UnsupervisedTrainer

Accuracy = namedtuple('Accuracy', ['unweighted', 'weighted', 'worst', 'accuracy_classes'])


class SemiSupervisedPosterior(Posterior):
    def __init__(self, *args, **kwargs):
        super(SemiSupervisedPosterior, self).__init__(*args, **kwargs)
        self.to_monitor = ["ll", "accuracy"]

    def accuracy(self, verbose=False):
        model, cls = (self.sampling_model, self.model) if hasattr(self, 'sampling_model') else (self.model, None)
        acc = compute_accuracy(model, self.data_loader, classifier=cls)
        if verbose:
            print("Acc for %s is : %.4f" % (self.name, acc))
        return acc

    accuracy.mode = 'max'

    def hierarchical_accuracy(self, name, verbose=False):

        all_y, all_y_pred = compute_predictions(self.model, self.data_loaders[name])
        acc = np.mean(all_y == all_y_pred)

        all_y_groups = np.array([self.model.labels_groups[y] for y in all_y])
        all_y_pred_groups = np.array([self.model.labels_groups[y] for y in all_y_pred])
        h_acc = np.mean(all_y_groups == all_y_pred_groups)

        if verbose:
            print("Acc for %s is : %.4f\nHierarchical Acc for %s is : %.4f\n" % (name, acc, name, h_acc))
        return acc

    accuracy.mode = 'max'

    def unsupervised_accuracy(self, verbose=False):
        uca = unsupervised_classification_accuracy(self.model, self.data_loader)[0]
        if verbose:
            print("UCA for %s is : %.4f" % (self.name, uca))
        return uca

    unsupervised_accuracy.mode = 'max'


def compute_accuracy_tuple(y, y_pred):
    y = y.ravel()
    n_labels = len(np.unique(y))
    classes_probabilities = []
    accuracy_classes = []
    for cl in range(n_labels):
        idx = y == cl
        classes_probabilities += [np.mean(idx)]
        accuracy_classes += [np.mean((y[idx] == y_pred[idx])) if classes_probabilities[-1] else 0]
        # This is also referred to as the "recall": p = n_true_positive / (n_false_negative + n_true_positive)
        # ( We could also compute the "precision": p = n_true_positive / (n_false_positive + n_true_positive) )
        accuracy_named_tuple = Accuracy(
            unweighted=np.dot(accuracy_classes, classes_probabilities),
            weighted=np.mean(accuracy_classes),
            worst=np.min(accuracy_classes),
            accuracy_classes=accuracy_classes)
    return accuracy_named_tuple


def compute_predictions(model, data_loader, classifier=None):
    all_y_pred = []
    all_y = []

    for i_batch, tensors in enumerate(data_loader):
        sample_batch, _, _, _, labels = tensors
        all_y += [labels.view(-1)]

        if hasattr(model, 'classify'):
            y_pred = model.classify(sample_batch).argmax(dim=-1)
        elif classifier is not None:
            # Then we use the specified classifier
            if model is not None:
                sample_batch, _, _ = model.z_encoder(sample_batch)
            y_pred = classifier(sample_batch).argmax(dim=-1)
        else:  # The model is the raw classifier
            y_pred = model(sample_batch).argmax(dim=-1)
        all_y_pred += [y_pred]

    all_y_pred = np.array(torch.cat(all_y_pred))
    all_y = np.array(torch.cat(all_y))
    return all_y, all_y_pred


def compute_accuracy(vae, data_loader, classifier=None):
    all_y, all_y_pred = compute_predictions(vae, data_loader, classifier=classifier)
    return np.mean(all_y == all_y_pred)


def unsupervised_classification_accuracy(vae, data_loader, classifier=None):
    all_y, all_y_pred = compute_predictions(vae, data_loader, classifier=classifier)
    return unsupervised_clustering_accuracy(all_y, all_y_pred)


def compute_accuracy_svc(data_train, labels_train, data_test, labels_test, unit_test=False, verbose=0,
                         max_iter=-1):
    param_grid = [
        {'C': [1, 10, 100, 1000], 'kernel': ['linear']},
        {'C': [1, 10, 100, 1000], 'gamma': [0.001, 0.0001], 'kernel': ['rbf']}]
    if unit_test:
        param_grid = [{'C': [1], 'kernel': ['linear']}]
    svc = SVC(max_iter=max_iter)

    clf = GridSearchCV(svc, param_grid, verbose=verbose)
    return compute_accuracy_classifier(clf, data_train, labels_train, data_test, labels_test)


def compute_accuracy_rf(data_train, labels_train, data_test, labels_test, unit_test=False,
                        verbose=0):
    param_grid = {'max_depth': np.arange(3, 10), 'n_estimators': [10, 50, 100, 200]}
    if unit_test:
        param_grid = [{'max_depth': [3], 'n_estimators': [10]}]

    rf = RandomForestClassifier(max_depth=2, random_state=0)

    clf = GridSearchCV(rf, param_grid, verbose=verbose)
    return compute_accuracy_classifier(clf, data_train, labels_train, data_test, labels_test)


def compute_accuracy_nn(data_train, labels_train, data_test, labels_test, k=5):
    clf = neighbors.KNeighborsClassifier(k, weights='distance')
    return compute_accuracy_classifier(clf, data_train, labels_train, data_test, labels_test)


def compute_accuracy_classifier(clf, data_train, labels_train, data_test, labels_test):
    clf.fit(data_train, labels_train)
    # Predicting the labels
    y_pred_test = clf.predict(data_test)
    y_pred_train = clf.predict(data_train)

    return (compute_accuracy_tuple(labels_train, y_pred_train),
            compute_accuracy_tuple(labels_test, y_pred_test)), y_pred_test


class ClassifierTrainer(Trainer):
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
        super(ClassifierTrainer, self).__init__(*args, use_cuda=use_cuda, **kwargs)
        self.train_set, self.test_set = self.train_test(self.model, self.gene_dataset, cls=SemiSupervisedPosterior)
        self.posteriors_dict = {"train_set": self.train_set, "test_set": self.test_set}
        self.posteriors_loop = ["train_set"]
        if sampling_model is not None:
            self.train_set.sampling_model = sampling_model
            self.test_set.sampling_model = sampling_model

    def loss(self, tensors_labelled):
        x, _, _, _, labels_train = tensors_labelled
        x = self.sampling_model.sample_from_posterior_z(x) if self.sampling_model is not None else x
        return F.cross_entropy(self.model(x), labels_train.view(-1))


class SemiSupervisedTrainer(UnsupervisedTrainer):
    r"""The SemiSupervisedTrainer class for the semi-supervised training of an autoencoder.
    This parent class can be inherited to specify the different training schemes for semi-supervised learning
    """
    default_metrics_to_monitor = Trainer.default_metrics_to_monitor + ['accuracy']

    def __init__(self, model, gene_dataset, n_labelled_samples_per_class=50, n_epochs_classifier=1,
                 lr_classification=0.1, classification_ratio=1, seed=0, **kwargs):
        """
        :param n_labelled_samples_per_class: number of labelled samples per class
        """
        super(SemiSupervisedTrainer, self).__init__(model, gene_dataset, **kwargs)
        self.model = model
        self.gene_dataset = gene_dataset

        self.n_epochs_classifier = n_epochs_classifier
        self.lr_classification = lr_classification
        self.classification_ratio = classification_ratio
        n_labelled_samples_per_class_array = [n_labelled_samples_per_class] * self.gene_dataset.n_labels
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

        self.full_dataset = self.create_posterior(shuffle=True, cls=SemiSupervisedPosterior,
                                                  name='full_dataset')
        self.labelled_set = self.create_posterior(indices=indices_labelled, cls=SemiSupervisedPosterior,
                                                  name='labelled_set')
        self.unlabelled_set = self.create_posterior(indices=indices_unlabelled, cls=SemiSupervisedPosterior,
                                                    name='unlabelled_set')

        self.classifier_trainer = ClassifierTrainer(
            model.classifier, gene_dataset, metrics_to_monitor=[], verbose=True, frequency=0,
            sampling_model=self.model
        )
        self.classifier_trainer.posteriors_dict['train_set'] = self.labelled_set
        self.classifier_trainer.posteriors_loop = ['train_set']

        self.posteriors_dict = {'labelled_set': self.labelled_set,
                                'unlabelled_set': self.unlabelled_set,
                                'full_dataset': self.full_dataset}

        self.posteriors_loop = ['full_dataset', 'labelled_set']
        for posterior in [self.labelled_set, self.unlabelled_set]:
            posterior.to_monitor = ['ll', 'accuracy']

    def loss(self, tensors_all, tensors_labelled):
        loss = super(SemiSupervisedTrainer, self).loss(tensors_all)
        sample_batch, _, _, _, y = tensors_labelled
        classification_loss = F.cross_entropy(self.model.classify(sample_batch), y.view(-1))
        loss += classification_loss * self.classification_ratio
        return loss

    def on_epoch_end(self):
        self.classifier_trainer.train(self.n_epochs_classifier, lr=self.lr_classification)
        return super(SemiSupervisedTrainer, self).on_epoch_end()

    def svc_rf(self, **kwargs):
        raw_data = Posterior.raw_data(self.labelled_set, self.unlabelled_set)
        (data_train, labels_train), (data_test, labels_test) = raw_data
        svc_scores, _ = compute_accuracy_svc(data_train, labels_train, data_test, labels_test, **kwargs)
        rf_scores, _ = compute_accuracy_rf(data_train, labels_train, data_test, labels_test, **kwargs)
        return svc_scores, rf_scores


class JointSemiSupervisedTrainer(SemiSupervisedTrainer):
    def __init__(self, model, gene_dataset, **kwargs):
        kwargs.update({'n_epochs_classifier': 0})
        super(JointSemiSupervisedTrainer, self).__init__(model, gene_dataset, **kwargs)


class AlternateSemiSupervisedTrainer(SemiSupervisedTrainer):
    def __init__(self, *args, **kwargs):
        super(AlternateSemiSupervisedTrainer, self).__init__(*args, **kwargs)
        self.posteriors_loop = ['full_dataset']

    def loss(self, all_tensor):
        return UnsupervisedTrainer.loss(self, all_tensor)
