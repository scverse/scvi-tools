import numpy as np


def train_test_split(*Xs, train_size=0.75):
    """
    A substitute for the sklearn function to avoid the dependency
    """
    Xs = [np.array(X) for X in Xs]
    split_idx = int(train_size * len(Xs[0]))

    split_list = [[X[:split_idx], X[split_idx:]] for X in Xs]
    return [X for Xs in split_list for X in Xs]
