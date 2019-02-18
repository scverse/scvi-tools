#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Test more memory efficient method to remove empty cells has same output."""


import numpy as np


def test_remove_empty_cells():
    """Test that our new method to remove empty cells is the same as the old."""

    X = np.random.poisson(lam=1, size=(1111, 3))

    def existing_method(X):
        to_keep = np.array((X.sum(axis=1) > 0)).ravel()
        removed_idx = []
        if X.shape != X[to_keep].shape:
            for i in range(len(to_keep)):
                if not to_keep[i]:
                    removed_idx.append(i)
        return X[to_keep], removed_idx

    def new_method(X):
        ne_cells = X.sum(axis=1) > 0
        if not ne_cells.all():
            X = X[np.where(ne_cells)]
        return X, np.where(~ne_cells)[0]

    #
    # Run existing method
    Xexisting, removed_idx = existing_method(X)
    if len(removed_idx):
        print("Cells with zero expression in all genes considered were removed, the indices of the removed cells "
              "in the expression matrix were:")
        print(removed_idx)
    #
    # Run new method
    Xnew, removed_idx = new_method(X)
    if len(removed_idx):
        print("Cells with zero expression in all genes considered were removed, the indices of the removed cells "
              "in the expression matrix were:")
        print(removed_idx)
    #
    # Check the same
    assert Xexisting.shape == Xnew.shape
    assert (Xexisting == Xnew).all()
