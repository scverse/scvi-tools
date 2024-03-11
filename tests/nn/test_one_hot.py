import pytest
import torch

from scvi.nn import one_hot
from scvi.nn._utils import _one_hot


def test_one_hot_deprecation_warning():
    """Public API should give deprecation warning"""
    x = torch.tensor([[1], [1], [2]])
    with pytest.warns(DeprecationWarning):
        one_hot(x, 3)


def test_one_hot():
    # N x 1
    x = torch.tensor([[1], [1], [2]])
    expected = torch.tensor([[0, 1, 0], [0, 1, 0], [0, 0, 1]])
    assert torch.equal(_one_hot(x, 3), expected)

    expected = torch.tensor([[0, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
    assert torch.equal(_one_hot(x, 4), expected)

    # N x M x 1
    x = torch.tensor([[[1], [2]], [[0], [1]]])
    expected = torch.tensor(
        [
            [[0, 1, 0], [0, 0, 1]],
            [[1, 0, 0], [0, 1, 0]],
        ]
    )
    assert torch.equal(_one_hot(x, 3), expected)

    # N x 2 – error
    x = torch.tensor([[1, 2], [0, 1]])
    with pytest.raises(ValueError):
        _one_hot(x, 3)
