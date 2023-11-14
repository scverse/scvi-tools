import os

import pytest

import scvi
from scvi.train._callbacks import MetricsCallback, SaveCheckpoint


def test_modelcheckpoint_callback(save_path: str):
    def check_checkpoint_logging(model, adata):
        assert any(isinstance(c, SaveCheckpoint) for c in model.trainer.callbacks)
        callback = [c for c in model.trainer.callbacks if isinstance(c, SaveCheckpoint)]
        assert len(callback) == 1
        callback = callback[0]
        assert callback.best_model_path is not None
        assert callback.best_model_score is not None
        assert os.path.exists(callback.best_model_path)

        log_dirs = os.listdir(scvi.settings.logging_dir)
        assert len(log_dirs) >= 1
        checkpoints_dir = os.path.join(scvi.settings.logging_dir, log_dirs[0])
        checkpoint_dirs = os.listdir(checkpoints_dir)
        assert len(checkpoint_dirs) >= 1
        checkpoint_dir = os.path.join(checkpoints_dir, checkpoint_dirs[0])
        checkpoint = model.__class__.load(checkpoint_dir, adata=adata)
        assert checkpoint.is_trained_

    def test_model_cls(model_cls, adata):
        scvi.settings.logging_dir = os.path.join(save_path, model_cls.__name__)

        # enable_checkpointing=True, default monitor
        model = model_cls(adata)
        model.train(max_epochs=5, enable_checkpointing=True)
        check_checkpoint_logging(model, adata)

        # enable_checkpointing=True, custom monitor
        model = model_cls(adata)
        model.train(
            max_epochs=5,
            enable_checkpointing=True,
            checkpointing_monitor="elbo_validation",
        )
        check_checkpoint_logging(model, adata)

        # manual callback
        model = model_cls(adata)
        model.train(
            max_epochs=5,
            callbacks=[SaveCheckpoint(monitor="elbo_validation")],
        )
        check_checkpoint_logging(model, adata)

    old_logging_dir = scvi.settings.logging_dir
    adata = scvi.data.synthetic_iid()

    scvi.model.SCVI.setup_anndata(adata)
    test_model_cls(scvi.model.SCVI, adata)

    scvi.model.SCANVI.setup_anndata(adata, "labels", "label_0")
    test_model_cls(scvi.model.SCANVI, adata)

    scvi.settings.logging_dir = old_logging_dir


def test_metricscallback_init():
    def dummy_metric(model) -> float:
        return 0.0

    callback = MetricsCallback(dummy_metric)
    assert callback.metric_fns == {"dummy_metric": dummy_metric}

    metrics = [dummy_metric]
    callback = MetricsCallback(metrics)
    assert callback.metric_fns == {"dummy_metric": dummy_metric}

    metrics = {"dummy_metric": dummy_metric, "dummy_metric2": dummy_metric}
    callback = MetricsCallback(metrics)
    assert len(callback.metric_fns) == 2

    with pytest.raises(TypeError):
        MetricsCallback(0)

    with pytest.raises(TypeError):
        MetricsCallback([0])


def test_metricscallback_with_scvi():
    adata = scvi.data.synthetic_iid()
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata)

    def dummy_metric(model: scvi.model.SCVI) -> float:
        _ = model.get_latent_representation()
        return 0.0

    model.train(
        max_epochs=1,
        check_val_every_n_epoch=1,
        additional_val_metrics={"dummy_metric": dummy_metric},
    )
    assert "dummy_metric" in model.history

    model = scvi.model.SCVI(adata)
    model.train(
        max_epochs=1,
        check_val_every_n_epoch=1,
        additional_val_metrics=[dummy_metric],
    )
    assert "dummy_metric" in model.history

    with pytest.warns(UserWarning):
        # no validation step
        model.train(
            max_epochs=1,
            check_val_every_n_epoch=None,
            additional_val_metrics=[dummy_metric],
        )
