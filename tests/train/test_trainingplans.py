import pytest

from scvi.model import SCVI, JaxSCVI
from scvi.train import JaxTrainingPlan, TrainingPlan
from scvi.train._trainingplans import _compute_kl_weight


@pytest.mark.parametrize(
    "current,n_warm_up,min_kl_weight,max_kl_weight,expected",
    [
        (0, 400, 0.0, 1.0, 0.0),
        (200, 400, 0.0, 1.0, 0.5),
        (400, 400, 0.0, 1.0, 1.0),
        (0, 400, 0.5, 1.0, 0.5),
        (200, 400, 0.5, 1.0, 0.75),
        (400, 400, 0.0, 1.0, 1.0),
        (400, 400, 0.0, 2.0, 2.0),
    ],
)
def test_compute_kl_weight_linear_annealing(
    current, n_warm_up, min_kl_weight, max_kl_weight, expected
):
    kl_weight = _compute_kl_weight(
        current, 1, n_warm_up, None, max_kl_weight, min_kl_weight
    )
    assert kl_weight == pytest.approx(expected)
    kl_weight = _compute_kl_weight(
        1, current, None, n_warm_up, max_kl_weight, min_kl_weight
    )
    assert kl_weight == pytest.approx(expected)


@pytest.mark.parametrize("max_kl_weight", [1.0, 2.0])
def test_compute_kl_weight_no_annealing(max_kl_weight):
    assert _compute_kl_weight(1, 1, None, None, max_kl_weight, 0.0) == max_kl_weight


def test_compute_kl_weight_min_greater_max():
    with pytest.raises(ValueError):
        _compute_kl_weight(1, 1, 400, None, 0.5, 1.0)


@pytest.mark.parametrize(
    "epoch,step,n_epochs_kl_warmup,n_steps_kl_warmup,expected",
    [
        (0, 100, 100, 100, 0.0),
        (50, 200, 100, 1000, 0.5),
        (100, 200, 100, 1000, 1.0),
    ],
)
def test_compute_kl_precedence(
    epoch, step, n_epochs_kl_warmup, n_steps_kl_warmup, expected
):
    kl_weight = _compute_kl_weight(
        epoch, step, n_epochs_kl_warmup, n_steps_kl_warmup, 1.0, 0.0
    )
    assert kl_weight == expected


def test_loss_args(synthetic_adata):
    """Test that self._loss_args is set correctly."""
    SCVI.setup_anndata(synthetic_adata)
    JaxSCVI.setup_anndata(synthetic_adata)
    vae = SCVI(synthetic_adata)
    jax_vae = JaxSCVI(synthetic_adata)
    tp = TrainingPlan(vae.module)
    jax_tp = JaxTrainingPlan(jax_vae.module)

    loss_args = [
        "tensors",
        "inference_outputs",
        "generative_outputs",
        "kl_weight",
    ]
    assert len(tp._loss_args) == len(loss_args)
    assert len(jax_tp._loss_args) == len(loss_args)
    for arg in loss_args:
        assert arg in tp._loss_args
        assert arg in jax_tp._loss_args
