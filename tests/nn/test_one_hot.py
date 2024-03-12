import pytest
import torch

from scvi.nn import one_hot


def test_one_hot_deprecation_warning():
    """Public API should give deprecation warning"""
    x = torch.tensor([[1], [1], [2]])
    with pytest.warns(DeprecationWarning):
        one_hot(x, 3)
