import scvi
from scvi.inference import auto_tune_scvi_model
from hyperopt import hp
import logging

logger = logging.getLogger("scvi.inference.autotune")
logger.setLevel(logging.WARNING)


def test_autotune_cortex(save_path):
    n_epochs = 1
    max_evals = 1
    reserve_timeout = 5
    fmin_timeout = 10

    cortex_dataset = scvi.dataset.cortex(save_path=save_path)

    best_trainer, trials = auto_tune_scvi_model(
        gene_dataset=cortex_dataset,
        parallel=False,
        exp_key="cortex_dataset",
        train_func_specific_kwargs={"n_epochs": n_epochs},
        max_evals=max_evals,
        reserve_timeout=reserve_timeout,
        fmin_timeout=fmin_timeout,
        save_path=save_path,  # temp dir, see conftest.py
    )

    space = {
        "model_tunable_kwargs": {"dropout_rate": hp.uniform("dropout_rate", 0.1, 0.3)},
        "train_func_tunable_kwargs": {"lr": hp.loguniform("lr", -4.0, -3.0)},
    }

    best_trainer, trials = auto_tune_scvi_model(
        gene_dataset=cortex_dataset,
        space=space,
        parallel=False,
        exp_key="cortex_dataset_custom_space",
        train_func_specific_kwargs={"n_epochs": n_epochs},
        max_evals=max_evals,
        reserve_timeout=reserve_timeout,
        fmin_timeout=fmin_timeout,
        save_path=save_path,
    )
