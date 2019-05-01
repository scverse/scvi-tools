import atexit
import datetime
import logging
import os
import pickle
import time

from collections import defaultdict
from functools import partial
from multiprocessing import Process, Queue
from subprocess import Popen
from typing import Any, Callable, Dict, List, Type, Union

from hyperopt import fmin, tpe, Trials, hp, STATUS_OK
from hyperopt.mongoexp import MongoTrials

import torch

from . import Trainer
from .inference import UnsupervisedTrainer
from ..dataset import GeneExpressionDataset
from ..models import VAE

# TODO: add database watcher and inform user about progress

# register running process and open files to terminate/close at exit
RUNNING_PROCESSES = []
OPEN_FILES = []

# initialize scVI logger
# if handler already added (by user) don't add default handler
# FIXME we force a unique handler..
logger = logging.getLogger(__name__)
if not len(logger.handlers):
    ch = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)


# Kill all subprocesses if parent dies
@atexit.register
def _cleanup_processes_files():
    """cleanup function"""
    logger.info("Cleaning up: terminating processes")
    for p in RUNNING_PROCESSES:
        p.terminate()
    logger.info("Cleaning up: closing files")
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
    save_path: str = ".",
    use_cpu: bool = False,
    n_workers_per_gpu: int = 1,
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
    :param save_path
    :param use_cpu: if True, also launch cpu-only workers
    :param n_workers_per_gpu: Number of workers ton launch per gpu found by torch
    :return: Trainer object for the best model and Trials object containing logs for the different runs.

    Examples:
        >>> gene_dataset = CortexDataset()
        >>> best_trainer, trials = auto_tuned_scvi_parameters(gene_dataset)
    """
    logger.info("Starting experiment: {exp_key}".format(exp_key=exp_key))
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
                "dropout_rate": hp.choice("dropout_rate", [0.1, 0.3, 0.5, 0.7, 0.9]),
                "reconstruction_loss": hp.choice("reconstruction_loss", ["zinb", "nb"]),
            },
            "train_func_tunable_kwargs": {"lr": hp.choice("lr", [0.01, 0.005, 0.001, 0.0005, 0.0001])},
        }
    )

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
        trials = _auto_tune_parallel(
            objective_hyperopt=objective_hyperopt,
            exp_key=exp_key,
            space=space,
            max_evals=max_evals,
            save_path=save_path,
            use_cpu=use_cpu,
            n_workers_per_gpu=n_workers_per_gpu,
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


def _auto_tune_parallel(
    objective_hyperopt: Callable,
    exp_key: str,
    space: dict = None,
    max_evals: int = 100,
    save_path: str = ".",
    use_cpu: bool = False,
    n_workers_per_gpu: int = 1,
) -> MongoTrials:
    # start by running fmin process so that workers don't timeout
    # run hyperoptimization, in a forked process
    # this allows to warn if the workers crash
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

    running_workers = []
    n_gpus = torch.cuda.device_count()
    if not n_gpus and not use_cpu:
        raise ValueError("No GPUs detected by PyTorch and use_cpu is set to False")

    # run mongo
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

    # run n hyperopt worker per gpu available
    logger.info(
        "Launching {n_workers_per_gpu} worker.s for each of the {n_gpus} gpu.s found by torch"
        "found".format(n_workers_per_gpu=n_workers_per_gpu, n_gpus=n_gpus)
    )
    for gpu_id in range(n_gpus):
        for sub_id in range(n_workers_per_gpu):
            p = _launch_hyperopt_worker(
                exp_key=exp_key, gpu=True, id=str(gpu_id), sub_id=str(sub_id)
            )
            running_workers.append(p)

    # run one hyperopt worker per cpu (though not specifically assigned)
    # minus one to prevent overloading the cores
    # FIXME set cpu affinity and use AVAILABLE CPUs not EXISTING
    logger.info(
        "Launching {n_cpu} cpu worker.s".format(
            n_cpu=max(0, os.cpu_count() * use_cpu - 1)
        )
    )
    for cpu_id in range(max(0, os.cpu_count() * use_cpu - 1)):
        p = _launch_hyperopt_worker(exp_key=exp_key, gpu=False, id=str(cpu_id))
        running_workers.append(p)
    RUNNING_PROCESSES.extend(running_workers)

    # wait and warn if all workers have died in case call to fmin hangs
    check_on_workers_process = Process(
        target=check_on_workers,
        kwargs={"running_workers": running_workers},
    )
    check_on_workers_process.start()
    RUNNING_PROCESSES.append(check_on_workers_process)

    # wait for fmin then close worker checking process
    fmin_process.join()
    check_on_workers_process.terminate()

    # get result, terminate all subprocesses and close logfile
    # fmin_process is done but queue might fail if no wait (see link to standard doc)
    # https://docs.python.org/3.6/library/multiprocessing.html#multiprocessing.Queue
    time.sleep(1)
    trials = queue.get_nowait()
    logger.info("Done with experiment {exp_key}, cleaning up".format(exp_key=exp_key))
    _cleanup_processes_files()

    return trials


def _launch_hyperopt_worker(
    exp_key:str, gpu: bool = True, id: str = None, sub_id: str = None
) -> Popen:
    device_type = "gpu" if gpu else "cpu"
    sub_id = ":" + sub_id if sub_id else str()
    worker_name = "worker_" + device_type + "_" + id + sub_id
    cuda_visible_devices = id if gpu else str()
    sub_env = {
        # FIXME pythonpath unnecessary if scVI installed
        "PYTHONPATH": ".",
        "CUDA_VISIBLE_DEVICES": cuda_visible_devices,
        "WORKER_NAME": worker_name,
    }
    sub_env = {**os.environ, **sub_env}
    # FIXME hard-coded worker params, not universal?
    p = Popen(
        [
            "hyperopt-mongo-worker",
            "--mongo=localhost:1234/scvi_db",
            "--exp-key={exp_key}".format(exp_key=exp_key),
            "--poll-interval=0.1",
            "--max-consecutive-failures=1",
            "--reserve-timeout=30.0",
        ],
        env=sub_env,
    )
    return p


def check_on_workers(running_workers: List[Popen]):
        for worker in running_workers:
            worker.communicate()
        logger.info(
            "All workers have died, check stdout/stderr for error tracebacks."
            "If ReserveTimeout was raised wait for fmin to finish. \n"
            "If fmin doesn't finish, terminate the process and check that "
            "the exp_key hasn't been used before or with a lesser max_evals."
        )


def _fmin_parallel(
    queue: Queue,
    fn: Callable,
    exp_key: str,
    space: dict,
    algo: Callable,
    max_evals: int,
    show_progressbar: bool,
):
    logger.info("Minimization process launched.")
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
    # queue.put uses pickle so remove attribute containing thread.lock
    if hasattr(trials, "handle"):
        logger.info("Deleting Trial handle for pickling")
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
    start_time = time.monotonic()
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
    elapsed_time = time.monotonic() - start_time
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
        metric_history = trainer.history[metric]
        logger.info(
            "[Worker : {name}] Training of {n_epochs} epochs finished in {time}".format(
                name=os.environ.get("WORKER_NAME", "0"),
                n_epochs=len(metric_history),
                time=str(datetime.timedelta(seconds=elapsed_time)),
            )
        )

        return {
            "loss": metric_history[-1],
            "elapsed_time": elapsed_time,
            "status": STATUS_OK,
            "history": trainer.history,
            "space": space,
        }
