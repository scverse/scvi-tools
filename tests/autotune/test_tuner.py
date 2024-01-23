from scvi.autotune import ModelTuner
from scvi.data import synthetic_iid
from scvi.model import SCVI


def test_model_tuner_init():
    ModelTuner(SCVI)


def test_model_tuner_fit(save_path: str):
    tuner = ModelTuner(SCVI)

    adata = synthetic_iid()
    SCVI.setup_anndata(adata)
    results = tuner.fit(
        adata,
        use_defaults=True,
        num_samples=1,
        max_epochs=1,
        logging_dir=save_path,
    )
    assert results is not None


def test_model_tuner_info():
    tuner = ModelTuner(SCVI)

    tuner.info()
