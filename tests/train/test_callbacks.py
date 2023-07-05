import pytest

import scvi
from scvi.train._callbacks import MetricsCallback


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
