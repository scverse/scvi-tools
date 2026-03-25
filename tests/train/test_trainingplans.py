import pytest

import scvi
from scvi.data import synthetic_iid
from scvi.model import SCVI
from scvi.train import TrainingPlan, _trainingplans
from scvi.train._constants import METRIC_KEYS
from scvi.train._trainingplans import _compute_kl_weight


@pytest.mark.parametrize(
    ("current", "n_warm_up", "min_kl_weight", "max_kl_weight", "expected"),
    [
        (0, 400, 0.0, 1.0, 0.0),
        (200, 400, 0.0, 1.0, 0.5),
        (400, 400, 0.0, 1.0, 1.0),
        (0, 400, 0.5, 1.0, 0.5),
        (200, 400, 0.5, 1.0, 0.75),
        (400, 400, 0.0, 2.0, 2.0),
    ],
)
def test_compute_kl_weight_linear_annealing(
    current, n_warm_up, min_kl_weight, max_kl_weight, expected
):
    kl_weight = _compute_kl_weight(current, 1, n_warm_up, None, max_kl_weight, min_kl_weight)
    assert kl_weight == pytest.approx(expected)
    kl_weight = _compute_kl_weight(1, current, None, n_warm_up, max_kl_weight, min_kl_weight)
    assert kl_weight == pytest.approx(expected)


@pytest.mark.parametrize("max_kl_weight", [1.0, 2.0])
def test_compute_kl_weight_no_annealing(max_kl_weight):
    assert _compute_kl_weight(1, 1, None, None, max_kl_weight, 0.0) == max_kl_weight


def test_compute_kl_weight_min_greater_max():
    with pytest.raises(ValueError):
        _compute_kl_weight(1, 1, 400, None, 0.5, 1.0)


@pytest.mark.parametrize(
    ("epoch", "step", "n_epochs_kl_warmup", "n_steps_kl_warmup", "expected"),
    [
        (0, 100, 100, 100, 0.0),
        (50, 200, 100, 1000, 0.5),
        (100, 200, 100, 1000, 1.0),
    ],
)
def test_compute_kl_precedence(epoch, step, n_epochs_kl_warmup, n_steps_kl_warmup, expected):
    kl_weight = _compute_kl_weight(epoch, step, n_epochs_kl_warmup, n_steps_kl_warmup, 1.0, 0.0)
    assert kl_weight == expected


def test_jax_import_guard_includes_flax():
    """Regression test: the Jax import guard in _trainingplans must check for flax
    in addition to jax and optax, because JaxBaseModuleClass is only defined when
    flax is installed (GH-3762)."""
    import ast
    import inspect

    source = inspect.getsource(_trainingplans)
    tree = ast.parse(source)

    # Find the top-level if-block that imports JaxBaseModuleClass
    for node in ast.walk(tree):
        if not isinstance(node, ast.If):
            continue
        # Check if the body contains an import of JaxBaseModuleClass
        body_source = ast.dump(node)
        if "JaxBaseModuleClass" not in body_source:
            continue
        # Verify the guard condition includes "flax"
        guard_source = ast.dump(node.test)
        assert "flax" in guard_source, (
            "The Jax import guard in _trainingplans.py must check for 'flax' "
            "to match _base_module.py (see GH-3762)"
        )
        break
    else:
        pytest.skip("JaxBaseModuleClass import guard not found in _trainingplans.py")


def test_loss_args():
    """Test that self._loss_args is set correctly."""
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)
    vae = SCVI(adata)
    tp = TrainingPlan(vae.module)

    loss_args = [
        "tensors",
        "inference_outputs",
        "generative_outputs",
        "kl_weight",
    ]
    assert len(tp._loss_args) == len(loss_args)
    for arg in loss_args:
        assert arg in tp._loss_args


def test_semisupervisedtrainingplan_metrics():
    adata = scvi.data.synthetic_iid(n_labels=3)
    scvi.model.SCANVI.setup_anndata(
        adata,
        labels_key="labels",
        unlabeled_category="label_0",
        batch_key="batch",
    )
    model = scvi.model.SCANVI(adata)
    model.train(max_epochs=1, check_val_every_n_epoch=1)

    for mode in ["train", "validation"]:
        for metric in [
            METRIC_KEYS.ACCURACY_KEY,
            METRIC_KEYS.F1_SCORE_KEY,
            METRIC_KEYS.CLASSIFICATION_LOSS_KEY,
        ]:
            assert f"{mode}_{metric}" in model.history_
