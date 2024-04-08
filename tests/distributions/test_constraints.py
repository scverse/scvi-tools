import torch

from scvi.distributions._constraints import open_interval


def test_open_interval():
    lower_bound = 0
    upper_bound = 1.0

    constraint = open_interval(lower_bound=lower_bound, upper_bound=upper_bound)

    # Within the interval is ok
    assert constraint.check(torch.tensor(0.5))

    # Values outside interval are not ok
    assert not constraint.check(torch.tensor(1.1))
    assert not constraint.check(torch.tensor(-0.1))

    # Endpoints are also not ok
    assert not constraint.check(torch.tensor(0))
    assert not constraint.check(torch.tensor(1.0))
