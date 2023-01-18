import scvi


def test_model_tuner_init():
    model_cls = scvi.model.SCVI
    scvi.autotune.ModelTuner(model_cls)


def test_model_tuner_fit():
    model_cls = scvi.model.SCVI
    tuner = scvi.autotune.ModelTuner(model_cls)

    adata = scvi.data.synthetic_iid()
    model_cls.setup_anndata(adata)
    results = tuner.fit(adata, num_samples=1, max_epochs=1)
    assert results is not None


def test_model_tuner_info():
    model_cls = scvi.model.SCVI
    tuner = scvi.autotune.ModelTuner(model_cls)

    tuner.info()
