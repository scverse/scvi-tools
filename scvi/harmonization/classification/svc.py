from scvi.metrics.classification import compute_accuracy_tuple

from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC as SVC_

class SVC:
    def fit(self, data_train, labels_train, verbose=False, unit_test=False, max_iter=-1):
        param_grid = [
            {'C': [1, 10, 100, 1000], 'kernel': ['linear']},
            {'C': [1, 10, 100, 1000], 'gamma': [0.001, 0.0001], 'kernel': ['rbf']}]
        if unit_test:
            param_grid = [{'C': [1], 'kernel': ['linear']}]
        svc = SVC_(max_iter=max_iter)

        self.clf = GridSearchCV(svc, param_grid, verbose=verbose)
        self.clf.fit(data_train, labels_train)
        self.accuracy_tuple_train = compute_accuracy_tuple(labels_train, self.clf.predict(data_train))
        return self.clf.score(data_train, labels_train)

    def score(self, data_test, labels_test):
        self.accuracy_tuple_test = compute_accuracy_tuple(labels_test, self.clf.predict(data_test))
        return self.clf.score(data_test, labels_test)
