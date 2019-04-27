import atexit
import logging
import os
import pickle
import time

from collections import defaultdict
from functools import partial
from multiprocessing import Process, Queue
from subprocess import Popen
from typing import Any, Callable, Dict, Type, Union

from hyperopt import fmin, tpe, Trials, hp, STATUS_OK
from hyperopt.mongoexp import MongoTrials

import torch

from . import Trainer
from .inference import UnsupervisedTrainer
from ..dataset import GeneExpressionDataset
from ..models import VAE


# register running process to terminate at exit
RUNNING_PROCESSES = []
OPEN_FILES = []

# initialize scVI logger
logging.basicConfig()
logger = logging.getLogger(__name__)


# Kill all subprocesses if parent dies
@atexit.register
def cleanup_processes_files():
    for p in RUNNING_PROCESSES:
        p.terminate()
    for f in OPEN_FILES:
        f.close()


def auto_tune_scvi_model(
    exp_key: str,
    gene_dataset: GeneExpressionDataset,
    model_class: VAE = VAE,
    trainer_class: Trainer = UnsupervisedTrainer,
    model_specific_kwargs: dict = None,
    trainer_specific_kwargs: dict = None,
    train_func_specific_kwargs: dict = None,
    space: dict = None,
    use_batches: bool = False,
    max_evals: int = 100,
    parallel: bool = False,
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

    :param exp_key: Name of the experiment in MongoDb.
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
    If already exists in db, hyperopt will run a numebr of trainings equal to
    the difference between current and previous max_evals.
    :param verbose: if True, output parameters used for each training to stdout.
    :param use_cpu: if True, also launch cpu-only workers
    :return: Trainer object for the best model and Trials object containing logs for the different runs.

    Examples:
        >>> gene_dataset = CortexDataset()
        >>> best_trainer, trials = auto_tuned_scvi_parameters(gene_dataset)
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
        trials = auto_tune_parallel(
            objective_hyperopt=objective_hyperopt,
            exp_key=exp_key,
            space=space,
            max_evals=max_evals,
            save_path=save_path,
            use_cpu=use_cpu,
            n_worker_per_gpu=n_worker_per_gpu,
        )

    else:
        trials = Trials()

        # run hyperoptimization
        _ = fmin(
            fn=objective_hyperopt,
            space=space,
            algo=tpe.suggest,
            max_evals=max_evals,
            trials=trials,
        )

    # return best model, trained
    best_space = trials.best_trial["result"]["space"]
    best_trainer = objective_hyperopt(best_space, is_best_training=True)

    # pickle trainer and save model (overkill?)
    with open(
        os.path.join(save_path, "best_trainer_{key}".format(key=exp_key)), "wb"
    ) as f:
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


def auto_tune_parallel(
    objective_hyperopt: Callable,
    exp_key: str,
    space: dict = None,
    max_evals: int = 100,
    save_path: str = ".",
    use_cpu: bool = False,
    n_worker_per_gpu: int = 1,
) -> MongoTrials:
    running_workers = []
    n_gpus = torch.cuda.device_count()
    if not n_gpus and not use_cpu:
        raise ValueError("No GPUs detected by PyTorch and use_cpu is set to False")

    # run mongo daemon
    mongo_path = os.path.join(save_path, "mongo")
    if not os.path.exists(mongo_path):
        os.makedirs(mongo_path)
    mongo_logfile = open(os.path.join(mongo_path, "mongo_logfile"), "w")
    OPEN_FILES.append(mongo_logfile)
    mongod_process = Popen(
        ["mongod", "--quiet", "--dbpath={path}".format(path=mongo_path), "--port=1234"],
        stdout=mongo_logfile,
    )
    RUNNING_PROCESSES.append(mongod_process)

    # run one hyperopt worker per gpu available
    print(str(time.time()) + "launching workers")
    for gpu_id in range(n_gpus):
        for sub_id in range(n_worker_per_gpu):
            p = launch_hyperopt_worker(gpu=True, id=str(gpu_id), sub_id=str(sub_id))
            running_workers.append(p)

    # run one hyperopt worker per cpu (though not specifically assigned)
    # minus two to prevent overloading loss
    # FIXME set cpu affinity and use all AVAILABLE CPUs minus one
    for cpu_id in range(max(0, os.cpu_count() * use_cpu - 1)):
        p = launch_hyperopt_worker(gpu=False, id=str(cpu_id))
        running_workers.append(p)
    RUNNING_PROCESSES.extend(running_workers)

    # run hyperoptimization, in a forked process
    # this allows to terminate if the workers crash
    # since mongo is not thread-safe, trials must be instantiated in fork
    queue = Queue()
    fmin_kwargs = {
        "queue": queue,
        "fn": objective_hyperopt,
        "exp_key": exp_key,
        "space": space,
        "algo": tpe.suggest,
        "max_evals": max_evals,
        "show_progressbar": False,  # progbar useless in parallel mode
    }
    fmin_process = Process(target=_fmin_parallel, kwargs=fmin_kwargs)
    fmin_process.start()
    RUNNING_PROCESSES.append(fmin_process)

    # wait and raise if all workers have died to prevent hanging
    while fmin_process.is_alive():

        def is_dead(p):
            return p.poll() is not None

        is_dead_statuses = list(map(is_dead, running_workers))
        print(is_dead_statuses)
        if all(is_dead_statuses):
            raise RuntimeError(
                "All workers have died, check stdout/stderr for tracebacks"
            )
        time.sleep(10)

    # get result, terminate all subprocesses and close logfile
    # fmin_process is done so no need to wait but queue might fail (see link below)
    # https://docs.python.org/3.6/library/multiprocessing.html#multiprocessing.Queue
    time.sleep(1)
    trials = queue.get_nowait()
    cleanup_processes_files()

    return trials


def launch_hyperopt_worker(
    gpu: bool = True, id: str = None, sub_id: str = None
) -> Popen:
    device_type = "gpu" if gpu else "cpu"
    sub_id = ":" + sub_id if sub_id else str()
    worker_name = "worker_" + device_type + "_" + id + sub_id
    # worker_log_file = open("log" + worker_name, "w")
    cuda_visible_devices = id if gpu else str()
    sub_env = {
        # FIXME pythonpath unnecessary if scVI installed
        "PYTHONPATH": ".",
        "CUDA_VISIBLE_DEVICES": cuda_visible_devices,
        "WORKER_NAME": worker_name,
    }
    sub_env = {**os.environ, **sub_env}
    p = Popen(
        [
            "hyperopt-mongo-worker",
            "--mongo=localhost:1234/scvi_db",
            "--poll-interval=0.1",
            "--max-consecutive-failures=1",
            "--reserve-timeout=10.0",
        ],
        env=sub_env,
        # stdout=worker_log_file,
        # stderr=worker_log_file,
    )
    return p


def _fmin_parallel(
    queue: Queue,
    fn: Callable,
    exp_key: str,
    space: dict,
    algo: Callable,
    max_evals: int,
    show_progressbar: bool,
):
    # instantiate Trials object
    trials = MongoTrials("mongo://localhost:1234/scvi_db/jobs", exp_key=exp_key)

    # run hyperoptimization
    _ = fmin(
        fn=fn,
        space=space,
        algo=algo,
        max_evals=max_evals,
        trials=trials,
        show_progressbar=show_progressbar,
    )
    # queue uses pickle so remove attribute containing thread.lock
    if hasattr(trials, "handle"):
        del trials.handle
    queue.put(trials)


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
            "Worker : {name}".format(name=os.environ.get("WORKER_NAME", "0")) + "\n"
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
