import warnings

import pytest

from scvi.model.base._training_mixin import _check_warmup


@pytest.mark.parametrize(
    "plan_kwargs,max_epochs,n_cells,batch_size",
    [
        ({}, 399, 1, 1),
        ({"n_epochs_kl_warmup": 200}, 100, 1, 1),
        ({"n_epochs_kl_warmup": None, "n_steps_kl_warmup": 100}, 1, 100, 10),
        ({"n_epochs_kl_warmup": None, "n_steps_kl_warmup": 12}, 1, 103, 10),
    ],
)
def test_unsupervised_training_mixin_warmup_warning(
    plan_kwargs, max_epochs, n_cells, batch_size
):
    with pytest.warns(UserWarning):
        _check_warmup(plan_kwargs, max_epochs, n_cells, batch_size)


@pytest.mark.parametrize(
    "plan_kwargs,max_epochs,n_cells,batch_size",
    [
        ({}, 400, 1, 1),
        ({"n_epochs_kl_warmup": 200}, 200, 1, 1),
        ({"n_epochs_kl_warmup": 200, "n_steps_kl_warmup": None}, 200, 1, 1),
        ({"n_epochs_kl_warmup": 200, "n_steps_kl_warmup": 100}, 200, 100, 10),
        ({"n_epochs_kl_warmup": None, "n_steps_kl_warmup": None}, 200, 1, 1),
        ({"n_epochs_kl_warmup": None}, 200, 1, 1),
        ({"n_epochs_kl_warmup": None, "n_steps_kl_warmup": 10}, 1, 102, 10),
    ],
)
def test_unsupervised_training_mixin_warmup_no_warning(
    plan_kwargs, max_epochs, n_cells, batch_size
):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _check_warmup(plan_kwargs, max_epochs, n_cells, batch_size)
