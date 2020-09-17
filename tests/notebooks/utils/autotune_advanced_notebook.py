import numpy as np
import torch

from collections import defaultdict
from hyperopt import STATUS_OK
from scvi.core.modules import SCANVAE
from scvi.core.trainers import SemiSupervisedTrainer
from sklearn.model_selection import train_test_split
from scvi.dataset._anndata import get_from_registry
from scvi import _CONSTANTS


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
    trainer_specific_kwargs["silent"] = True
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
    trainer_scanvi.unlabelled_set = trainer_scanvi.create_scvi_dl(
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
            trainer_scanvi.labelled_set = trainer_scanvi.create_scvi_dl(
                indices=indices_labelled_train
            )
            trainer_scanvi.labelled_set.to_monitor = [
                "reconstruction_error",
                "accuracy",
            ]
            trainer_scanvi.validation_set = trainer_scanvi.create_scvi_dl(
                indices=indices_labelled_val
            )
            trainer_scanvi.validation_set.to_monitor = ["accuracy"]
            trainer_scanvi.train(**train_func_tunable_kwargs)
            accuracies[i] = trainer_scanvi.history["accuracy_unlabelled_set"][-1]
        return {"loss": -accuracies.mean(), "space": space, "status": STATUS_OK}
    else:
        trainer_scanvi.labelled_set = trainer_scanvi.create_scvi_dl(
            indices=indices_labelled
        )
        trainer_scanvi.labelled_set.to_monitor = ["reconstruction_error", "accuracy"]
        trainer_scanvi.train(**train_func_tunable_kwargs)
        return trainer_scanvi
