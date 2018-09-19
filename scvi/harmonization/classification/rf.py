from scvi.metrics.classification import compute_accuracy_tuple
import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier

class RF:
    def fit(self, data_train, labels_train, verbose=False, unit_test=False, max_iter=-1):
        param_grid = {'max_depth': np.arange(3, 10), 'n_estimators': [10, 50, 100, 200]}
        if unit_test:
            param_grid = [{'max_depth': [3], 'n_estimators': [10]}]

        rf = RandomForestClassifier(max_depth=2, random_state=0)

        self.clf = GridSearchCV(rf, param_grid, verbose=verbose)
        self.clf.fit(data_train, labels_train)

        self.accuracy_tuple_train = compute_accuracy_tuple(labels_train, self.clf.predict(data_train))

    def score(self, data_test, labels_test):
        y_pred_test = self.clf.predict(data_test)
        self.accuracy_tuple_test = self.clf.fit(data_test, y_pred_test)

        return self.clf.score(data_test, labels_test)

