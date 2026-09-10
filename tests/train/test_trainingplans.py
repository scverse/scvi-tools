import pytest
import torch

import scvi
from scvi.data import synthetic_iid
from scvi.model import SCVI
from scvi.train import TrainingPlan
from scvi.train._constants import METRIC_KEYS
from scvi.train._trainingplans import (
    _compilation_fell_back,
    _compute_kl_weight,
    _dynamo_frame_counts,
)


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


@pytest.mark.parametrize(
    ("before", "after", "expected"),
    [
        ((0, 0), (0, 0), False),  # dynamo never ran
        ((0, 0), (3, 3), False),  # everything compiled
        ((0, 0), (3, 1), False),  # partial success is still a compiled model
        ((0, 0), (3, 0), True),  # every frame fell back to eager
        ((5, 4), (8, 4), True),  # counters are process-wide, so compare deltas
        ((5, 4), (8, 5), False),
    ],
)
def test_compilation_fell_back_only_when_every_attempted_frame_failed(before, after, expected):
    assert _compilation_fell_back(before, after) is expected


@pytest.mark.optional
def test_compile_restores_the_global_dynamo_error_suppression():
    # `compile=True` turns dynamo's error suppression on so that a failed compilation
    # degrades to eager instead of raising. That is process-wide state, so it has to be
    # put back once training is over rather than leaking into unrelated later code.
    # Actually runs torch.compile on CPU, which is slow and not the intended use case
    # (compile is meant for mps) — kept out of the base suite, covered in --optional.
    import torch

    previous = torch._dynamo.config.suppress_errors
    torch._dynamo.config.suppress_errors = False
    try:
        adata = synthetic_iid(batch_size=32, n_genes=16)
        SCVI.setup_anndata(adata)
        model = SCVI(adata, n_latent=2, n_hidden=8, n_layers=1)
        model.train(
            max_epochs=1,
            accelerator="cpu",
            devices=1,
            enable_progress_bar=False,
            plan_kwargs={"compile": True},
        )
        assert torch._dynamo.config.suppress_errors is False
    finally:
        torch._dynamo.config.suppress_errors = previous


def test_compile_warns_when_every_frame_fell_back_to_eager(monkeypatch):
    # A silent fallback is the dangerous case: the user asked for a compiled model,
    # got an eager one, and nothing said so.
    from scvi.train import _trainingplans

    adata = synthetic_iid(batch_size=32, n_genes=16)
    SCVI.setup_anndata(adata)
    model = SCVI(adata, n_latent=2, n_hidden=8, n_layers=1)
    plan = TrainingPlan(model.module, compile=True)

    monkeypatch.setattr(_trainingplans, "_dynamo_frame_counts", lambda: (99, 0))
    plan._compile_frames_before = (0, 0)
    with pytest.warns(UserWarning, match="fell back to eager"):
        plan.on_train_end()


@pytest.mark.skipif(not torch.backends.mps.is_available(), reason="requires an MPS device")
def test_compile_actually_compiles_frames_on_mps():
    # The point of the warning above is that a silent eager fallback is indistinguishable
    # from a compiled run without checking. Assert the positive case on mps.
    adata = synthetic_iid(batch_size=32, n_genes=16)
    SCVI.setup_anndata(adata)
    model = SCVI(adata, n_latent=2, n_hidden=8, n_layers=1)
    before = _dynamo_frame_counts()
    model.train(
        max_epochs=1,
        accelerator="mps",
        devices=1,
        enable_progress_bar=False,
        plan_kwargs={"compile": True},
    )
    after = _dynamo_frame_counts()
    assert after[1] > before[1]
    assert not _compilation_fell_back(before, after)


def _make_plan(**kwargs):
    adata = synthetic_iid(batch_size=32, n_genes=16)
    SCVI.setup_anndata(adata)
    model = SCVI(adata, n_latent=2, n_hidden=8, n_layers=1)
    return TrainingPlan(model.module, **kwargs)


def test_fused_optimizer_is_used_when_the_backend_supports_it():
    plan = _make_plan()
    optimizer = plan.get_optimizer_creator()(
        p for p in plan.module.parameters() if p.requires_grad
    )
    assert optimizer.defaults.get("fused") is True


def test_fused_optimizer_can_be_turned_off():
    plan = _make_plan(fused_optimizer=False)
    optimizer = plan.get_optimizer_creator()(
        p for p in plan.module.parameters() if p.requires_grad
    )
    assert not optimizer.defaults.get("fused")


def test_fused_optimizer_falls_back_when_the_backend_rejects_it():
    # torch raises at construction time when fused is not available for these params,
    # and the message differs across versions and devices, so the fallback is driven by
    # the exception rather than by a device allowlist of our own.
    class RejectsFused(torch.optim.Adam):
        def __init__(self, params, **kwargs):
            if kwargs.pop("fused", False):
                raise RuntimeError("fused is not supported on this device")
            super().__init__(params, **kwargs)

    plan = _make_plan()
    params = [p for p in plan.module.parameters() if p.requires_grad]
    optimizer = plan._optimizer_creator_fn(RejectsFused)(params)
    assert not optimizer.defaults.get("fused")
    assert isinstance(optimizer, RejectsFused)
