import atexit
import logging
import os
import pickle

from collections import defaultdict
from functools import partial
from subprocess import Popen
from typing import Any, Dict, Type, Union

from hyperopt import fmin, tpe, Trials, hp, STATUS_OK
from hyperopt.mongoexp import MongoTrials

import torch

from . import Trainer
from .inference import UnsupervisedTrainer
from ..dataset import GeneExpressionDataset
from ..models import VAE


RUNNING_PROCESSES = []

# initialize scVI logger
logging.basicConfig()
logger = logging.getLogger(__name__)


# Kill all subprocesses if parent dies
@atexit.register
def cleanup():
    map(lambda p: p.terminate(), RUNNING_PROCESSES)


def auto_tuned_scvi_model(
    gene_dataset: GeneExpressionDataset,
    model_class: VAE = VAE,
    trainer_class: Trainer = UnsupervisedTrainer,
    model_specific_kwargs: dict = None,
    trainer_specific_kwargs: dict = None,
    train_func_specific_kwargs: dict = None,
    space: dict = None,
    use_batches: bool = False,
    max_evals: int = 10,
    parallel: bool = False,
    exp_key: str = "exp1",
    verbose: bool = True,
    save_path: str = ".",
    use_cpu: bool = False,
    n_worker_per_gpu: int = 1,
) -> (Type[Trainer], Trials):
    """Perform automatic hyperparameter optimization of an scVI model
    and return best model and hyperopt Trials object.
    Trials object contains hyperparameter space and loss history for each trial.
    We provide a default hyperparameter search space (see source code),
    but we recommend the user to build a custom one for each application.
    Convention: fixed parameters (no default) have precedence over tunable parameters (default).

    :param gene_dataset: scVI gene dataset
    :param model_class: scVI model class (e.g ``VAE``, ``VAEC``, ``SCANVI``)
    :param trainer_class: Trainer class (e.g ``UnsupervisedTrainer``)
    :param model_specific_kwargs: dict of fixed parameters which will be passed to the model.
    :param trainer_specific_kwargs: dict of fixed parameters which will be passed to the trainer.
    :param train_func_specific_kwargs: dict of fixed parameters which will be passed to the train method.
    :param space: dict containing up to three sub-dicts with keys "model_tunable_kwargs",
    "trainer_tunable_kwargs" or "train_func_tunable_kwargs".
    Each of those dict contains hyperopt defined parameter spaces (e.g. ``hp.choice(..)``)
    which will be passed to the corresponding object : model, trainer or train method
    when performing hyperoptimization. Default: mutable, see source code.
    :param use_batches: If False, pass n_batch=0 to model else pass gene_dataset.n_batches
    :param max_evals: Maximum number of trainings.
    :param parallel: If True, use MongoTrials object to run trainings in parallel.
    :param exp_key: Name of the experiment in MongoDb.
    If already exists in db, hyperopt will run a numebr of trainings equal to
    the difference between current and previous max_evals.
    :param verbose: if True, output parameters used for each training to stdout.
    :param use_cpu: if True, also launch cpu-only workers
    :return: Trainer object for the best model and Trials object containing logs for the different runs.

    Examples:
        >>> gene_dataset = CortexDataset()
        >>> best_trainer, trials = auto_tuned_scvi_parameters(gene_dataset)
        >>> posterior = best_trainer.create_posterior()
    """
    # default specific kwargs
    model_specific_kwargs = model_specific_kwargs if model_specific_kwargs else {}
    trainer_specific_kwargs = trainer_specific_kwargs if trainer_specific_kwargs else {}
    train_func_specific_kwargs = (
        train_func_specific_kwargs if train_func_specific_kwargs else {}
    )

    # default early stopping
    if "early_stopping_kwargs" not in trainer_specific_kwargs:
        early_stopping_kwargs = {
            "early_stopping_metric": "ll",
            "save_best_state_metric": "ll",
            "patience": 50,
            "threshold": 3,
        }
        trainer_specific_kwargs["early_stopping_kwargs"] = early_stopping_kwargs

    # default search space
    space = (
        space
        if space
        else {
            "model_tunable_kwargs": {
                "n_latent": 5 + hp.randint("n_latent", 11),  # [5, 15]
                "n_hidden": hp.choice("n_hidden", [64, 128, 256]),
                "n_layers": 1 + hp.randint("n_layers", 5),
                "dropout_rate": hp.uniform("dropout_rate", 0.1, 0.9),
                "reconstruction_loss": hp.choice("reconstruction_loss", ["zinb", "nb"]),
            },
            "train_func_tunable_kwargs": {"lr": hp.choice("lr", [0.01, 0.001, 0.0001])},
        }
    )

    # logging
    if verbose:
        logger.setLevel(logging.INFO)

    logger.info(
        "Fixed parameters: \n"
        "model: \n"
        + str(model_specific_kwargs)
        + "\n"
        + "trainer: \n"
        + str(trainer_specific_kwargs)
        + "\n"
        + "train method: \n"
        + str(train_func_specific_kwargs)
        + "\n"
    )

    # build a partial objective function restricted to the search space
    objective_hyperopt = partial(
        _objective_function,
        **{
            "gene_dataset": gene_dataset,
            "model_class": model_class,
            "trainer_class": trainer_class,
            "model_specific_kwargs": model_specific_kwargs,
            "trainer_specific_kwargs": trainer_specific_kwargs,
            "train_func_specific_kwargs": train_func_specific_kwargs,
            "use_batches": use_batches,
        },
    )
    if parallel:
        n_gpus = torch.cuda.device_count()
        if not n_gpus and not use_cpu:
            raise ValueError("No GPUs detected by PyTorch and use_cpu is set to False")

        # filter out logs for clarity
        hp_logger = logging.getLogger("hyperopt.mongoexp")
        hp_logger.addFilter(logging.Filter("scvi"))

        # run mongo daemon
        mongo_path = os.path.join(os.path.abspath("."), "mongo")
        if not os.path.exists(mongo_path):
            os.makedirs(mongo_path)
        mongo_logfile = open(os.path.join(mongo_path, "mongo_logfile"), "w")
        RUNNING_PROCESSES.append(
            Popen(
                [
                    "mongod",
                    "--quiet",
                    "--dbpath={path}".format(path=mongo_path),
                    "--port=1234",
                ],
                stdout=mongo_logfile,
            )
        )

        # run one hyperopt worker per gpu available
        for gpu in range(n_gpus):
            for i in range(n_worker_per_gpu):
                sub_env = {
                    "PYTHONPATH": ".",
                    "CUDA_VISIBLE_DEVICES": str(gpu),
                    "WORKER_NAME": "worker_gpu_{i}".format(i=str(gpu) + ":" + str(i)),
                }
                sub_env = {**os.environ, **sub_env}
                RUNNING_PROCESSES.append(
                    Popen(
                        [
                            "hyperopt-mongo-worker",
                            "--mongo=localhost:1234/scvi_db",
                            "--poll-interval=0.1",
                            "--max-consecutive-failures=1",
                            "--reserve-timeout=10.0",
                        ],
                        env=sub_env,
                    )
                )

        # run one hyperopt worker per cpu (though not specifically assigned)
        # minus two to prevent overloading loss
        # FIXME set cpu affinity and use all existing CPUs minus one
        for cpu in range(max(0, os.cpu_count() * use_cpu - 1)):
            sub_env = {
                "PYTHONPATH": ".",
                "CUDA_VISIBLE_DEVICES": "",
                "WORKER_NAME": "worker_cpu_{i}".format(i=cpu),
            }
            sub_env = {**os.environ, **sub_env}
            RUNNING_PROCESSES.append(
                Popen(
                    [
                        "hyperopt-mongo-worker",
                        "--mongo=localhost:1234/scvi_db",
                        "--poll-interval=0.1",
                        "--max-consecutive-failures=1",
                        "--reserve-timeout=10.0",
                    ],
                    env=sub_env,
                )
            )

    # instantiate Trials object
    trials = (
        MongoTrials("mongo://localhost:1234/scvi_db/jobs", exp_key=exp_key)
        if parallel
        else Trials()
    )

    # run hyperoptimization
    _ = fmin(
        objective_hyperopt,
        space=space,
        algo=tpe.suggest,
        max_evals=max_evals,
        trials=trials,
        show_progressbar=not parallel,  # progbar useless in parallel mode
    )

    # kill all subprocesses and close logfile
    map(lambda p: p.terminate(), RUNNING_PROCESSES)
    mongo_logfile.close()

    # return best model, trained
    best_space = trials.best_trial["result"]["space"]
    best_trainer = objective_hyperopt(best_space, is_best_training=True)

    # pickle trainer and save model (overkill?)
    with open(os.path.join(save_path, "best_trainer_{key}".format(key=exp_key)), "wb") as f:
        pickle.dump(best_trainer, f)
    torch.save(
        best_trainer.model.state_dict(),
        os.path.join(save_path, "best_model_{key}".format(key=exp_key)),
    )

    # remove object containing thread.lock (otherwise pickle.dump throws)
    if hasattr(trials, "handle"):
        del trials.handle
    with open(os.path.join(save_path, "trials_{key}".format(key=exp_key)), "wb") as f:
        pickle.dump(trials, f)

    return best_trainer, trials


def _objective_function(
    space: dict,
    gene_dataset: GeneExpressionDataset,
    model_class: Type[VAE] = VAE,
    trainer_class: Type[Trainer] = UnsupervisedTrainer,
    model_specific_kwargs: dict = None,
    trainer_specific_kwargs: dict = None,
    train_func_specific_kwargs: dict = None,
    use_batches: bool = False,
    is_best_training: bool = False,
) -> Union[Dict[str, Any], Type[VAE]]:
    """Objective function for automatic hyperparameter optimization.
    Train a scVI model and return the best value of the early-stopping metric (e.g, log-likelihood).
    Convention: fixed parameters (no default) have precedence over tunable parameters (default).

    :param space: dict containing up to three sub-dicts with keys "model_tunable_kwargs",
    "trainer_tunable_kwargs" or "train_func_tunable_kwargs".
    Each of those dict contains hyperopt defined parameter spaces (e.g. ``hp.choice(..)``)
    which will be passed to the corresponding object : model, trainer or train method
    when performing hyperoptimization.
    :param gene_dataset: scVI gene dataset
    :param model_class: scVI model class (e.g ``VAE``, ``VAEC``, ``SCANVI``)
    :param trainer_class: Trainer class (e.g ``UnsupervisedTrainer``)
    :param model_specific_kwargs: dict of fixed parameters which will be passed to the model.
    :param trainer_specific_kwargs: dict of fixed parameters which will be passed to the trainer.
    :param train_func_specific_kwargs: dict of fixed parameters which will be passed to the train method.
    :param use_batches: If False, pass n_batch=0 to model else pass gene_dataset.n_batches
    :param is_best_training: True if training the model with the best hyperparameters
    :return: best value of the early stopping metric, and best model if is_best_training
    """
    # hyperopt params
    space = defaultdict(dict, space)
    model_tunable_kwargs = space["model_tunable_kwargs"]
    trainer_tunable_kwargs = space["trainer_tunable_kwargs"]
    train_func_tunable_kwargs = space["train_func_tunable_kwargs"]

    # use_cuda default
    if "use_cuda" not in trainer_specific_kwargs:
        trainer_specific_kwargs["use_cuda"] = bool(torch.cuda.device_count())
    if "n_epochs" not in {**train_func_specific_kwargs, **train_func_tunable_kwargs}:
        train_func_specific_kwargs["n_epochs"] = 1000

    # add hardcoded parameters
    # disable scVI progbar
    trainer_specific_kwargs["show_progbar"] = False
    if is_best_training:
        trainer_specific_kwargs["train_size"] = 1.0
        # no monitoring, will crash otherwise
        trainer_specific_kwargs["frequency"] = None
        trainer_specific_kwargs["early_stopping_kwargs"] = {}
    else:
        # evaluate at each epoch
        trainer_specific_kwargs["frequency"] = 1

    # merge params with fixed param precedence
    model_tunable_kwargs.update(model_specific_kwargs)
    trainer_tunable_kwargs.update(trainer_specific_kwargs)
    train_func_tunable_kwargs.update(train_func_specific_kwargs)

    if not is_best_training:
        # FIXME had to do info level because hyperopt basicConfig
        # is info level in MongoWorker
        logger.info(
            "Worker : {name}".format(name=os.environ["WORKER_NAME"]) + "\n"
            "Parameters being tested: \n"
            "model: \n"
            + str(model_tunable_kwargs)
            + "\n"
            + "trainer: \n"
            + str(trainer_tunable_kwargs)
            + "\n"
            + "train method: \n"
            + str(train_func_tunable_kwargs)
            + "\n"
        )

    # define model
    model = model_class(
        n_input=gene_dataset.nb_genes,
        n_batch=gene_dataset.n_batches * use_batches,
        **model_tunable_kwargs,
    )

    # define trainer
    trainer = trainer_class(model, gene_dataset, **trainer_tunable_kwargs)

    # train model
    trainer.train(**train_func_tunable_kwargs)

    # if training the best model, return model else return criterion
    if is_best_training:
        return trainer
    else:
        # select metric from early stopping kwargs or default to "ll_test_set"
        metric = None
        if "early_stopping_kwargs" in trainer_specific_kwargs:
            early_stopping_kwargs = trainer_specific_kwargs["early_stopping_kwargs"]
            if "early_stopping_metric" in early_stopping_kwargs:
                metric = early_stopping_kwargs["early_stopping_metric"]
                metric += "_" + trainer.early_stopping.on
        metric = metric if metric else "ll_test_set"

        return {
            "loss": trainer.history[metric][-1],
            "status": STATUS_OK,
            "history": trainer.history,
            "space": space,
        }
