import datetime
import logging
import os
import pickle
import sys
import time

from collections import defaultdict
from functools import partial, wraps
from io import IOBase
from multiprocessing import Process, Queue
from subprocess import Popen
from typing import Any, Callable, Dict, List, Type, Union

from hyperopt import fmin, tpe, Trials, hp, STATUS_OK
from hyperopt.mongoexp import as_mongo_str, MongoJobs, MongoTrials, MongoWorker, ReserveTimeout, Shutdown, WaitQuit

import torch

from . import Trainer
from .inference import UnsupervisedTrainer
from ..dataset import GeneExpressionDataset
from ..models import VAE

# TODO: add database watcher and inform user about progress, choose adequate logging levels

# register running process and open files to terminate/close at exit
RUNNING_PROCESSES: List[Union[Process, Popen]] = []
OPEN_FILES: List[IOBase] = []

# initialize scVI logger
# if handler already added (by user) don't add default handler
# FIXME we force a unique handler..
logger = logging.getLogger(__name__)
if not len(logger.handlers):
    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    ch.setFormatter(formatter)
    logger.addHandler(ch)


# cleanup helper functions
def _cleanup_processes_files():
    """cleanup function, starts with latest processes/files"""
    logger.info("Cleaning up: terminating processes")
    for p in RUNNING_PROCESSES[::-1]:
        if isinstance(p, Popen):
            if p.poll() is not None:
                p.terminate()
        if isinstance(p, Process):
            if p.is_alive():
                p.terminate()
    logger.info("Cleaning up: closing files")
    for f in OPEN_FILES[::-1]:
        if not f.closed:
            f.close()


def _cleanup_decorator(func):
    @wraps(func)
    def decorated(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception:
            logger.info("Caught Exception, starting cleanup")
            _cleanup_processes_files()
            raise
    return decorated


@_cleanup_decorator
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
    reserve_timeout: float = 30.0,
    mongo_port: str = "1234",
    mongo_host: str = "localhost",
    db_name:str = "scvi_db",

) -> (Type[Trainer], Trials):
    """Perform automatic hyperparameter optimization of an scVI model
    and return best model and hyperopt Trials object.
    Trials object contains hyperparameter space and loss history for each trial.
    We provide a default hyperparameter search space (see source code),
    but we recommend the user to build a custom one for each application.
    Convention: fixed parameters (no default) have precedence over tunable parameters (default).

    :param exp_key: Name of the experiment in MongoDb.
    :param gene_dataset: scVI gene dataset.
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
    :param max_evals: Maximum number of evaluations of the objective
    :param parallel: If True, use MongoTrials object to run trainings in parallel.
        If already exists in db, hyperopt will run a numebr of trainings equal to
        the difference between current and previous max_evals.
    :param save_path: Path where to save best model, trainer, trials and mongo files.
    :param use_cpu: if True, also launch cpu-only workers.
    :param n_workers_per_gpu: Number of workers ton launch per gpu found by torch.
    :return: Trainer object for the best model and Trials object containing logs for the different runs.

    Examples:
        >>> from scvi.dataset import CortexDataset
        >>> gene_dataset = CortexDataset()
        >>> best_trainer, trials = auto_tune_scvi_model(gene_dataset)
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
            "train_func_tunable_kwargs": {
                "lr": hp.choice("lr", [0.01, 0.005, 0.001, 0.0005, 0.0001])
            },
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
        logger.info("Starting parallel hyperoptimization")
        trials = _auto_tune_parallel(
            objective_hyperopt=objective_hyperopt,
            exp_key=exp_key,
            space=space,
            max_evals=max_evals,
            save_path=save_path,
            use_cpu=use_cpu,
            n_workers_per_gpu=n_workers_per_gpu,
            reserve_timeout=reserve_timeout,
            mongo_port=mongo_port,
            mongo_host=mongo_host,
            db_name=db_name,
        )

    else:
        logger.info("Starting sequential hyperoptimization")
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
    logger.info("Training best model with full training set")
    best_space = trials.best_trial["result"]["space"]
    best_trainer = objective_hyperopt(best_space, is_best_training=True)

    # pickle trainer and save model (overkill?)
    logger.info("Pickling best model, trainer as well as Trials object")
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
    reserve_timeout: float = 30.0,
    mongo_port: str = "1234",
    mongo_host: str = "localhost",
    db_name: str = "scvi_db",
) -> MongoTrials:

    """Parallel version of the hyperoptimization procedure.
    Called by auto_tune_scvi_model when parallel=True

    :param objective_hyperopt: Callable, the objective function to minimize
    :param exp_key: Name of the experiment in MongoDb.
    :param space: dict containing up to three sub-dicts with keys "model_tunable_kwargs",
        "trainer_tunable_kwargs" or "train_func_tunable_kwargs".
        Each of those dict contains hyperopt defined parameter spaces (e.g. ``hp.choice(..)``)
        which will be passed to the corresponding object : model, trainer or train method
        when performing hyperoptimization. Default: mutable, see source code.
    :param max_evals: Maximum number of evaluations of the objective.
    :param save_path: Path where to save best model, trainer, trials and mongo files.
    :param use_cpu: if True, also launch cpu-only workers.
    :param n_workers_per_gpu: Number of workers ton launch per gpu found by torch.
    :return: MongoTrials object containing the results of the procedure.
    """
    # run mongo
    logger.info("Starting MongoDb process")
    mongo_path = os.path.join(save_path, "mongo")
    if not os.path.exists(mongo_path):
        os.makedirs(mongo_path)
    mongo_logfile = open(os.path.join(mongo_path, "mongo_logfile"), "w")
    OPEN_FILES.append(mongo_logfile)
    mongod_process = Popen(
        ["mongod", "--quiet", "--dbpath={path}".format(path=mongo_path), "--port={port}".format(port=mongo_port)],
        stdout=mongo_logfile,
    )
    mongo_port_address = os.path.join(mongo_host + ":" + mongo_port, db_name)
    RUNNING_PROCESSES.append(mongod_process)

    # start by running fmin process so that workers don't timeout
    # run hyperoptimization, in a forked process
    # this allows to warn if the workers crash
    # since mongo is not thread-safe, trials must be instantiated in fork
    logger.info("Starting minimization procedure")
    queue = Queue()
    fmin_kwargs = {
        "queue": queue,
        "fn": objective_hyperopt,
        "exp_key": exp_key,
        "space": space,
        "algo": tpe.suggest,
        "max_evals": max_evals,
        "show_progressbar": False,  # progbar useless in parallel mode
        "mongo_port_address": mongo_port_address,
    }
    fmin_process = Process(target=_fmin_parallel, kwargs=fmin_kwargs)
    fmin_process.start()
    RUNNING_PROCESSES.append(fmin_process)

    # start worker launcher
    logger.info("Starting worker launcher")
    launcher_kwargs = {
        "exp_key": exp_key,
        "use_cpu": use_cpu,
        "n_workers_per_gpu": n_workers_per_gpu,
        "reserve_timeout": reserve_timeout,
        "workdir": mongo_path,
        "mongo_port_address": mongo_port_address,
    }
    workers_process = Process(
        target=launch_workers,
        kwargs=launcher_kwargs,
    )
    workers_process.start()
    RUNNING_PROCESSES.append(workers_process)

    # wait for fmin then close worker process
    trials = queue.get(block=True)
    fmin_process.join()
    logger.info("Finished minimization procedure for experiment {exp_key}".format(exp_key=exp_key))
    RUNNING_PROCESSES.remove(fmin_process)
    logger.info("Terminating worker launcher")
    workers_process.terminate()
    logger.info("worker launcher terminated")
    RUNNING_PROCESSES.remove(workers_process)

    # get result, terminate all forked processed and close files
    _cleanup_processes_files()

    return trials


@_cleanup_decorator
def _fmin_parallel(
    queue: Queue,
    fn: Callable,
    exp_key: str,
    space: dict,
    algo: Callable = tpe.suggest,
    max_evals: int = 100,
    show_progressbar: bool = False,
    mongo_port_address: str = "localhost:1234/scvi_db",
):
    """Launches a minimization procedure
    """
    # instantiate Trials object
    trials = MongoTrials(as_mongo_str(os.path.join(mongo_port_address, "jobs")), exp_key=exp_key)

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
    logger.info("Putting Trials in Queue")
    queue.put(trials)


@_cleanup_decorator
def launch_workers(
    exp_key: str,
    use_cpu: bool = False,
    n_workers_per_gpu: int = 1,
    reserve_timeout: float = 30.0,
    workdir: str = ".",
    mongo_port_address: str = "localhost:1234/scvi_db",
):
    n_gpus = torch.cuda.device_count()
    if not n_gpus and not use_cpu:
        raise ValueError("No GPUs detected by PyTorch and use_cpu is set to False")

    running_workers = []
    # run n hyperopt worker per gpu available
    logger.info(
        "Starting {n_workers_per_gpu} worker.s for each of the {n_gpus} gpu.s found by torch"
        "found".format(n_workers_per_gpu=n_workers_per_gpu, n_gpus=n_gpus)
    )
    for gpu_id in range(n_gpus):
        for sub_id in range(n_workers_per_gpu):
            worker_kwargs = {
                "exp_key": exp_key,
                "workdir": workdir,
                "gpu": True,
                "hw_id": str(gpu_id),
                "sub_id": str(sub_id),
                "reserve_timeout": reserve_timeout,
                "mongo_port_address": mongo_port_address,
            }
            p = Process(
                target=_hyperopt_worker,
                kwargs=worker_kwargs,
            )
            # FIXME won't terminate if parent is killed (SIGKILL)
            p.start()
            running_workers.append(p)

    # run one hyperopt worker per cpu (though not specifically assigned)
    # minus one to prevent overloading the cores
    # FIXME set cpu affinity and use AVAILABLE CPUs not EXISTING
    logger.info(
        "Starting {n_cpu} cpu worker.s".format(
            n_cpu=max(0, os.cpu_count() * use_cpu - 1)
        )
    )
    for cpu_id in range(max(0, os.cpu_count() * use_cpu - 1)):
        worker_kwargs = {
            "exp_key": exp_key,
            "workdir": workdir,
            "gpu": False,
            "hw_id": str(cpu_id),
            "reserve_timeout": reserve_timeout,
            "mongo_port_address": mongo_port_address,
        }
        p = Process(
            target=_hyperopt_worker,
            kwargs=worker_kwargs,
        )
        # FIXME won't terminate if parent is killed (SIGKILL)
        p.start()
        running_workers.append(p)
    RUNNING_PROCESSES.extend(running_workers)

    # wait and warn if all workers have died in case call to fmin hangs
    workers_watchdog(running_workers=running_workers)


def _hyperopt_worker(
    exp_key: str,
    workdir: str = ".",
    gpu: bool = True,
    hw_id: str = None,
    sub_id: str = None,
    poll_interval: float = 1.0,
    reserve_timeout: float = 30.0,
    mongo_port_address: str = "localhost:1234/scvi_db",
):
    """Launches a hyperopt MongoWorker which runs jobs until ReserveTimeout is raised.
    """
    # worker naming convention
    device_type = "gpu" if gpu else "cpu"
    sub_id = ":" + sub_id if sub_id else str()
    os.environ["WORKER_NAME"] = "worker_" + device_type + "_" + hw_id + sub_id

    os.environ["CUDA_VISIBLE_DEVICES"] = hw_id if gpu else str()
    # FIXME is this stil necessary?
    sys.path.append(".")

    mjobs = MongoJobs.new_from_connection_str(
        os.path.join(as_mongo_str(mongo_port_address), 'jobs')
    )
    mworker = MongoWorker(mjobs,
                          float(poll_interval),
                          workdir=workdir,
                          exp_key=exp_key)

    while True:
        # FIXME we don't protect ourselves from memory leaks, bad cleanup, etc.
        # The tradeoff is that a large dataset must be reloaded once for
        # each subprocess.
        try:
            mworker.run_one(reserve_timeout=float(reserve_timeout))
        except ReserveTimeout:
            logger.info(
                "Caught ReserveTimeout. "
                "Exiting after failing to reserve job for {time} seconds".format(time=reserve_timeout)
            )
            break


def workers_watchdog(running_workers: List[Process]):
    """Checks that worker in running_workers are stil running, if not inform user.
    """
    while True:
        one_alive = False
        for worker in running_workers:
            one_alive = one_alive or worker.is_alive()
        # if all workers are dead, inform user
        if not one_alive:
            logger.info(
                "All workers have died, check stdout/stderr for error tracebacks."
                "If ReserveTimeout was raised wait for fmin to finish. \n"
                "If fmin doesn't finish, terminate the process and check that "
                "the exp_key hasn't been used before or is used with a higher max_evals."
            )
        time.sleep(5)


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
) -> Union[Dict[str, Any], Trainer]:
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
        logger.info(
            "[{name}]".format(name=os.environ.get("WORKER_NAME", "0"))
            + "Parameters being tested: \n"
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
    logger.info("Instantiating model")
    model = model_class(
        n_input=gene_dataset.nb_genes,
        n_batch=gene_dataset.n_batches * use_batches,
        **model_tunable_kwargs,
    )

    # define trainer
    logger.info("Instantiating trainer")
    trainer = trainer_class(model, gene_dataset, **trainer_tunable_kwargs)

    # train model
    logger.info("Starting training")
    trainer.train(**train_func_tunable_kwargs)
    logger.info("Finished training")
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
        worker_name = os.environ.get("WORKER_NAME", "0")
        logger.info(
            "[{name}] Training of {n_epochs} epochs finished in {time}".format(
                name=worker_name,
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
            "worker_name": worker_name,
        }
