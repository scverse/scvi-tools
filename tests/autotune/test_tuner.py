import os
import shutil

import scvi


def test_model_tuner_init():
    model_cls = scvi.model.SCVI
    scvi.autotune.ModelTuner(model_cls)


def test_model_tuner_fit():
    model_cls = scvi.model.SCVI
    tuner = scvi.autotune.ModelTuner(model_cls)

    adata = scvi.data.synthetic_iid()
    model_cls.setup_anndata(adata)
    logging_dir = os.path.join(os.getcwd(), "tests", "data", "autotune_logs")
    results = tuner.fit(
        adata,
        use_defaults=True,
        num_samples=1,
        max_epochs=1,
        logging_dir=logging_dir,
    )
    assert results is not None
    shutil.rmtree(logging_dir)


def test_model_tuner_info():
    model_cls = scvi.model.SCVI
    tuner = scvi.autotune.ModelTuner(model_cls)

    tuner.info()
