import pytest

import scvi


@pytest.mark.optional
def test_model_tuner_init():
    model_cls = scvi.model.SCVI
    scvi.autotune.ModelTuner(model_cls)


@pytest.mark.optional
def test_model_tuner_fit(save_path):
    model_cls = scvi.model.SCVI
    tuner = scvi.autotune.ModelTuner(model_cls)

    adata = scvi.data.synthetic_iid()
    model_cls.setup_anndata(adata)
    results = tuner.fit(
        adata,
        use_defaults=True,
        num_samples=1,
        max_epochs=1,
        logging_dir=save_path,
    )
    assert results is not None


@pytest.mark.optional
def test_model_tuner_info():
    model_cls = scvi.model.SCVI
    tuner = scvi.autotune.ModelTuner(model_cls)

    tuner.info()
