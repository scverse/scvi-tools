import scvi


def test_metricscallback_with_scvi():
    adata = scvi.data.synthetic_iid()
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata)

    def dummy_metric(model: scvi.model.SCVI) -> float:
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
