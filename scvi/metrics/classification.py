from collections import namedtuple

import numpy as np
import torch
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.utils.linear_assignment_ import linear_assignment
from sklearn import neighbors

Accuracy = namedtuple('Accuracy', ['unweighted', 'weighted', 'worst', 'accuracy_classes'])


def compute_accuracy_tuple(y, y_pred):
    y = y.ravel()
    classes_probabilities = []
    accuracy_classes = []
    for cl in np.unique(y):
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


def compute_predictions(vae, data_loader, classifier=None):
    all_y_pred = []
    all_y = []

    for i_batch, tensors in enumerate(data_loader):
        sample_batch, _, _, _, labels = tensors
        all_y += [labels.view(-1)]

        if classifier is not None:
            # Then we use the specified classifier
            if vae is not None:
                sample_batch, _, _ = vae.z_encoder(torch.log(1+sample_batch))
            y_pred = classifier(sample_batch).argmax(dim=-1)
        elif hasattr(vae, 'classify'):
            y_pred = vae.classify(sample_batch).argmax(dim=-1)
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


def unsupervised_clustering_accuracy(y, y_pred):
    """
    Unsupervised Clustering Accuracy
    """
    assert len(y_pred) == len(y)
    n_clusters = len(np.unique(y))
    reward_matrix = np.zeros((n_clusters, n_clusters), dtype=np.int64)
    for y_pred_, y_ in zip(y_pred, y):
        reward_matrix[y_pred_, y_] += 1
    cost_matrix = reward_matrix.max() - reward_matrix
    ind = linear_assignment(cost_matrix)
    return sum([reward_matrix[i, j] for i, j in ind]) * 1.0 / y_pred.size, ind


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
