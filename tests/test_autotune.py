import os

from hyperopt import hp

from scvi.dataset import CortexDataset, LoomDataset, PbmcDataset
from scvi.inference.autotune import auto_tune_scvi_model
from notebooks.utils.autotune_advanced_notebook import custom_objective_hyperopt

n_epochs = 1
n_epochs_brain_large = 1
max_evals = 1
reserve_timeout = 5
fmin_timeout = 10


def test_default_usages(save_path):
    # Default
    cortex_dataset = CortexDataset(save_path=save_path)
    best_trainer, trials = auto_tune_scvi_model(
        gene_dataset=cortex_dataset,
        parallel=False,
        exp_key="cortex_dataset",
        train_func_specific_kwargs={"n_epochs": n_epochs},
        max_evals=max_evals,
        reserve_timeout=reserve_timeout,
        fmin_timeout=fmin_timeout,
    )

    # Custom space
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
    )

    # Custom objective metric
    pbmc_dataset = PbmcDataset(
        save_path=save_path, save_path_10X=os.path.join(save_path, "10X")
    )
    best_trainer, trials = auto_tune_scvi_model(
        gene_dataset=pbmc_dataset,
        metric_name="entropy_batch_mixing",
        posterior_name="train_set",
        parallel=False,
        exp_key="pbmc_entropy_batch_mixing",
        train_func_specific_kwargs={"n_epochs": n_epochs},
        max_evals=max_evals,
        reserve_timeout=reserve_timeout,
        fmin_timeout=fmin_timeout,
    )


def test_custom_fn(save_path):
    synthetic_dataset = LoomDataset(
        filename="simulation_1.loom",
        save_path=os.path.join(save_path, "simulation/"),
        url="https://github.com/YosefLab/scVI-data/raw/master/simulation/simulation_1.loom",
    )

    objective_kwargs = dict(dataset=synthetic_dataset, n_epochs=n_epochs)
    best_trainer, trials = auto_tune_scvi_model(
        custom_objective_hyperopt=custom_objective_hyperopt,
        objective_kwargs=objective_kwargs,
        parallel=False,
        exp_key="synthetic_dataset_scanvi",
        max_evals=max_evals,
        reserve_timeout=reserve_timeout,
        fmin_timeout=fmin_timeout,
    )
