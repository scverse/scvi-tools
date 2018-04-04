import numpy as np


def dropout(X, rate=0.1):
    """
    X: original testing set
    ========
    returns:
    X_zero: copy of X with zeros
    i, j, ix: indices of where dropout is applied
    """
    X_zero = np.copy(X)
    # select non-zero subset
    i, j = np.nonzero(X_zero)

    # choice number 1 : select 10 percent of the non zero values (so that distributions overlap enough)
    ix = np.random.choice(range(len(i)), int(np.floor(rate * len(i))), replace=False)
    X_zero[i[ix], j[ix]] *= 0  # *np.random.binomial(1, rate)

    return X_zero, i, j, ix


def imputation_error(X_pred, X, i, j, ix):
    """
    X_pred: imputed dataset
    X: original dataset
    X_zero: zeros dataset
    i, j, ix: indices of where dropout was applied
    ========
    returns:
    median L1 distance between datasets at indices given
    """
    all_index = i[ix], j[ix]
    x, y = X_pred[all_index], X[all_index]
    return np.median(np.abs(x - y))
