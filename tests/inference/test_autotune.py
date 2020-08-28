import scvi
from scvi.inference import auto_tune_scvi_model
from scvi import _CONSTANTS
import os
from hyperopt import hp
from hyperopt import STATUS_OK
import torch
from collections import defaultdict
import numpy as np

from scvi.core.trainers import SemiSupervisedTrainer
from scvi.core.models import SCANVAE
from sklearn.model_selection import train_test_split
from scvi.dataset._anndata import get_from_registry


import logging

logger = logging.getLogger("scvi.inference.autotune")
logger.setLevel(logging.WARNING)

n_epochs = 1
max_evals = 1
reserve_timeout = 5
fmin_timeout = 10


def test_autotune_cortex(save_path):

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


def custom_objective_hyperopt(
    space, is_best_training=False, dataset=None, n_epochs=None
):
    """Custom objective function for advanced autotune tutorial."""
    space = defaultdict(dict, space)
    model_tunable_kwargs = space["model_tunable_kwargs"]
    trainer_tunable_kwargs = space["trainer_tunable_kwargs"]
    train_func_tunable_kwargs = space["train_func_tunable_kwargs"]

    trainer_specific_kwargs = {}
    model_specific_kwargs = {}
    train_func_specific_kwargs = {}
    trainer_specific_kwargs["use_cuda"] = bool(torch.cuda.device_count())
    train_func_specific_kwargs["n_epochs"] = n_epochs

    # add hardcoded parameters
    # disable scVI progbar
    trainer_specific_kwargs["show_progbar"] = False
    trainer_specific_kwargs["frequency"] = 1

    # merge params with fixed param precedence
    model_tunable_kwargs.update(model_specific_kwargs)
    trainer_tunable_kwargs.update(trainer_specific_kwargs)
    train_func_tunable_kwargs.update(train_func_specific_kwargs)

    scanvi = SCANVAE(
        dataset.uns["_scvi"]["summary_stats"]["n_genes"],
        dataset.uns["_scvi"]["summary_stats"]["n_batch"],
        dataset.uns["_scvi"]["summary_stats"]["n_labels"],
        **model_tunable_kwargs
    )
    trainer_scanvi = SemiSupervisedTrainer(scanvi, dataset, **trainer_tunable_kwargs)
    batch_indices = get_from_registry(dataset, _CONSTANTS.BATCH_KEY)
    trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(
        indices=(batch_indices == 1)
    )
    trainer_scanvi.unlabelled_set.to_monitor = ["reconstruction_error", "accuracy"]
    indices_labelled = batch_indices == 0

    if not is_best_training:
        # compute k-fold accuracy on a 20% validation set
        k = 5
        accuracies = np.zeros(k)
        indices_labelled = batch_indices == 0
        for i in range(k):
            indices_labelled_train, indices_labelled_val = train_test_split(
                indices_labelled.nonzero()[0], test_size=0.2
            )
            trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(
                indices=indices_labelled_train
            )
            trainer_scanvi.labelled_set.to_monitor = [
                "reconstruction_error",
                "accuracy",
            ]
            trainer_scanvi.validation_set = trainer_scanvi.create_posterior(
                indices=indices_labelled_val
            )
            trainer_scanvi.validation_set.to_monitor = ["accuracy"]
            trainer_scanvi.train(**train_func_tunable_kwargs)
            accuracies[i] = trainer_scanvi.history["accuracy_unlabelled_set"][-1]
        return {"loss": -accuracies.mean(), "space": space, "status": STATUS_OK}
    else:
        trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(
            indices=indices_labelled
        )
        trainer_scanvi.labelled_set.to_monitor = ["reconstruction_error", "accuracy"]
        trainer_scanvi.train(**train_func_tunable_kwargs)
        return trainer_scanvi


def test_autotune_simulated_data(save_path):

    synthetic_dataset = scvi.dataset.annotation_simulation(
        1, save_path=os.path.join(save_path, "simulation/")
    )
    objective_kwargs = dict(dataset=synthetic_dataset, n_epochs=n_epochs)
    best_trainer, trials = auto_tune_scvi_model(
        custom_objective_hyperopt=custom_objective_hyperopt,
        objective_kwargs=objective_kwargs,
        parallel=True,
        exp_key="synthetic_dataset_scanvi",
        max_evals=max_evals,
        reserve_timeout=reserve_timeout,
        fmin_timeout=fmin_timeout,
        save_path=save_path,
    )
