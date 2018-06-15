from collections import namedtuple

import numpy as np
import torch
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC

from scvi.utils import no_grad, eval_modules, to_cuda

Accuracy = namedtuple('Accuracy', ['unweighted', 'weighted', 'worst', 'accuracy_classes'])


def compute_accuracy_tuple(y, labels):
    labels = labels.ravel()
    n_labels = len(np.unique(labels))
    classes_probabilities = []
    accuracy_classes = []
    for cl in range(n_labels):
        idx = labels == cl
        classes_probabilities += [np.mean(idx)]
        accuracy_classes += [np.mean((labels[idx] == y[idx])) if classes_probabilities[-1] else 0]
        # This is also referred to as the "recall": p = n_true_positive / (n_false_negative + n_true_positive)
        # ( We could also compute the "precision": p = n_true_positive / (n_false_positive + n_true_positive) )
        accuracy_named_tuple = Accuracy(
            unweighted=np.dot(accuracy_classes, classes_probabilities),
            weighted=np.mean(accuracy_classes),
            worst=np.min(accuracy_classes),
            accuracy_classes=accuracy_classes)
    return accuracy_named_tuple


@no_grad()
@eval_modules()
def compute_accuracy(vae, data_loader, classifier=None):
    all_y_pred = []
    all_labels = []

    for i_batch, tensors in enumerate(data_loader):
        if vae.use_cuda:
            tensors = to_cuda(tensors)
        sample_batch, _, _, _, labels = tensors
        sample_batch = sample_batch.type(torch.float32)
        all_labels += [labels.view(-1)]

        if classifier is not None:
            # Then we use the specified classifier
            mu_z, _, _ = vae.z_encoder(sample_batch)
            y_pred = classifier(mu_z).argmax(dim=-1)
        else:
            # Then the vae must implement a classify function
            y_pred = vae.classify(sample_batch).argmax(dim=-1)
        all_y_pred += [y_pred]

    accuracy = (torch.cat(all_y_pred) == torch.cat(all_labels)).type(torch.float32).mean().item()

    return accuracy


def compute_accuracy_svc(data_train, labels_train, data_test, labels_test, unit_test=False, verbose=0):
    param_grid = [
        {'C': [1, 10, 100, 1000], 'kernel': ['linear']},
        {'C': [1, 10, 100, 1000], 'gamma': [0.001, 0.0001], 'kernel': ['rbf']}]
    if unit_test:
        param_grid = [{'C': [1], 'kernel': ['linear']}]
    svc = SVC()

    clf = GridSearchCV(svc, param_grid, verbose=verbose)
    clf.fit(data_train, labels_train)

    # Predicting the labels
    y_pred_test = clf.predict(data_test)
    y_pred_train = clf.predict(data_train)

    return (compute_accuracy_tuple(y_pred_train, labels_train),
            compute_accuracy_tuple(y_pred_test, labels_test))


def compute_accuracy_rf(data_train, labels_train, data_test, labels_test, unit_test=False, verbose=0):
    param_grid = {'max_depth': np.arange(3, 10), 'n_estimators': [10, 50, 100, 200]}
    if unit_test:
        param_grid = [{'max_depth': [3], 'n_estimators': [10]}]

    rf = RandomForestClassifier(max_depth=2, random_state=0)

    clf = GridSearchCV(rf, param_grid, verbose=verbose)
    clf.fit(data_train, labels_train)

    # Predicting the labels
    y_pred_test = clf.predict(data_test)
    y_pred_train = clf.predict(data_train)

    return (compute_accuracy_tuple(y_pred_train, labels_train),
            compute_accuracy_tuple(y_pred_test, labels_test))
