import datetime
import logging
import multiprocessing
import os
import pickle
import threading
import time
from collections import defaultdict
from functools import partial, wraps
from logging.handlers import QueueListener, QueueHandler
from queue import Empty
from subprocess import Popen
from typing import Any, Callable, Dict, List, TextIO, Type, Union

import numpy as np
import pymongo
import torch
import tqdm
from hyperopt import fmin, tpe, Trials, hp, STATUS_OK, STATUS_FAIL
from hyperopt.mongoexp import (
    as_mongo_str,
    MongoJobs,
    MongoTrials,
    MongoWorker,
    ReserveTimeout,
)

from scvi._settings import autotune_formatter
from scvi.dataset import DownloadableDataset, GeneExpressionDataset
from scvi.models import VAE
from . import Trainer, UnsupervisedTrainer

# spawning is required for processes relying on cuda, and for windows
multiprocessing.set_start_method("spawn", force=True)

# instantiate logger, handler and formatter
# logger_all is used to send *all* autotune logs to a logfile
logger_all = logging.getLogger(__name__ + ".all")
logger_all.setLevel(logging.DEBUG)
logger = logging.getLogger(__name__)
# instantiate hyperopt and autotune file handlers as global variables for clean up
fh_autotune = None
fh_hyperopt = None


class FminTimeoutError(Exception):
    """Thrown if fmin process hasn't finished in the allotted
    time after all workers have died.
    """


class DispatchHandler(logging.Handler):
    """A simple dispatcher for logging events.

    It dispatches events to loggers based on the name in the received record,
    which then get dispatched, by the logging system, to the handlers, configured for those loggers.
    """

    def emit(self, record: logging.LogRecord):
        record_logger = logging.getLogger(record.name)
        if record.levelno >= record_logger.level:
            record_logger.handle(record)


class StoppableThread(threading.Thread):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.stop_event = threading.Event()

    def stop(self):
        self.stop_event.set()


# register running process and open files to terminate/close at exit
started_processes: List[Union[multiprocessing.Process, Popen, QueueListener]] = []
started_threads: List[StoppableThread] = []
started_queues: List[multiprocessing.Queue] = []
open_files: List[TextIO] = []


# cleanup helpers
def _cleanup_processes_files():
    """Cleanup function, starts with latest processes/files.

    Terminates processes, sets stop events to stop threads, closes open files.
    """
    logger_all.info("Cleaning up")
    logger_all.debug("Cleaning up: closing files.")
    for f in open_files[::-1]:
        if not f.closed:
            f.close()
    logger_all.debug("Cleaning up: closing queues.")
    for q in started_queues:
        q.close()
    logger_all.debug("Cleaning up: setting cleanup_event and joining threads.")
    for t in started_threads[::-1]:
        if t.is_alive():
            logger_all.debug("Closing Thread {}.".format(t.name))
            t.stop_event.set()
            t.join()
        else:
            logger_all.debug("Thread {} already done.".format(t.name))
    logger_all.debug("Cleaning up: terminating processes.")
    for p in started_processes[::-1]:
        if isinstance(p, Popen):
            if p.poll() is not None:
                logger_all.debug("Terminating mongod process.")
                p.terminate()
                p.wait()
            else:
                logger_all.debug("mongodd process already done.")
        if isinstance(p, multiprocessing.Process):
            if p.is_alive():
                logger_all.debug("Terminating Process {}.".format(p.name))
                p.terminate()
            else:
                logger_all.debug("Process {} already done.".format(p.name))
        if isinstance(p, QueueListener):
            if p._thread is not None and not p.queue._closed:
                p.stop()


def _cleanup_logger():
    """Removes added handlers."""
    logger_all.debug("Cleaning up: removing added logging handler.")
    hp_logger = logging.getLogger("hyperopt")
    for handler in hp_logger.handlers:
        if handler == fh_hyperopt:
            logger_all.debug("Cleaning up: removing hyperopt FileHandler.")
            hp_logger.removeHandler(fh_hyperopt)
            break
    for handler in logger_all.handlers:
        if handler == fh_autotune:
            logger_all.debug("Cleaning up: removing autotune FileHandler.")
            logger_all.removeHandler(fh_autotune)


def _cleanup_decorator(func: Callable):
    """Decorates top-level calls in order to launch cleanup when an Exception is caught."""

    @wraps(func)
    def decorated(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logger_all.exception(
                "Caught {exception} in {func}, starting cleanup.".format(
                    exception=e.args, func=func.__name__
                ),
                exc_info=True,
            )
            _cleanup_processes_files()
            _cleanup_logger()
            raise

    return decorated


def _error_logger_decorator(func: Callable):
    """Decorates top-level calls in order to launch cleanup when an Exception is caught."""

    @wraps(func)
    def decorated(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logger_all.exception(
                "Caught {exception} in {func}, starting cleanup.".format(
                    exception=e.args, func=func.__name__
                ),
                exc_info=True,
            )
            raise

    return decorated


def configure_asynchronous_logging(logging_queue: multiprocessing.Queue):
    """Helper for asynchronous logging - Writes all logs to a queue."""
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    queue_handler = QueueHandler(logging_queue)
    queue_handler.setLevel(logging.DEBUG)
    root_logger.addHandler(queue_handler)
    logger_all.debug("Asynchronous logging has been set.")


def _asynchronous_logging_method_decorator(func: Callable):
    """Decorates top-level calls in order to launch cleanup when an Exception is caught."""

    @wraps(func)
    def decorated(self, *args, **kwargs):
        configure_asynchronous_logging(self.logging_queue)
        return func(self, *args, **kwargs)

    return decorated


@_cleanup_decorator
def auto_tune_scvi_model(
    exp_key: str,
    gene_dataset: GeneExpressionDataset = None,
    delayed_populating: bool = False,
    custom_objective_hyperopt: Callable = None,
    objective_kwargs: Dict[str, Any] = None,
    model_class: VAE = VAE,
    trainer_class: Trainer = UnsupervisedTrainer,
    metric_name: str = None,
    metric_kwargs: Dict[str, Any] = None,
    posterior_name: str = "test_set",
    model_specific_kwargs: dict = None,
    trainer_specific_kwargs: dict = None,
    train_func_specific_kwargs: dict = None,
    space: dict = None,
    max_evals: int = 100,
    train_best: bool = True,
    pickle_result: bool = True,
    save_path: str = ".",
    use_batches: bool = False,
    parallel: bool = True,
    n_cpu_workers: int = None,
    gpu_ids: List[int] = None,
    n_workers_per_gpu: int = 1,
    reserve_timeout: float = 180.0,
    fmin_timeout: float = 300.0,
    fmin_timer: float = None,
    mongo_port: str = "1234",
    mongo_host: str = "localhost",
    db_name: str = "scvi_db",
    multiple_hosts: bool = False,
) -> (Type[Trainer], Trials):
    """Perform automatic hyperparameter optimization of an scVI model
    and return best model and hyperopt Trials object.

    ``Trials`` object contains hyperparameter space and loss history for each trial.
    We provide a default hyperparameter search space (see source code),
    but we recommend the user to build a custom one for each application.
    Convention: fixed parameters (no default) have precedence over tunable parameters (default).
    Note that the verbosity of this function has to be set using the logging module.
    In particular, for the parallel case, only a progress bar is shown if the
    logging level is equal or higher to ``logging.WARNING``.

    :param exp_key: Name of the experiment in MongoDb.
        If already exists in db, ``hyperopt`` will run a number of trainings equal to
        the difference between current and previous ``max_evals``.
    :param gene_dataset: scVI gene expression dataset.
    :param delayed_populating: Switch for the delayed populating mechanism
        of scvi.dataset.dataset.DownloadableDataset. Useful for large datasets
        which have to be instantiated inside the workers.
    :param custom_objective_hyperopt: A custom objective function respecting the ``hyperopt`` format.
        Roughly, it needs to return the quantity to optimize for, either directly
        or in a ``dict`` under the "loss" key.
        See https://github.com/hyperopt/hyperopt/wiki for a more detailed explanation.
        By default, we provide an objective function which can be parametrized
        through the various arguments of this function (``gene_dataset``, ``model_class``, etc.)
    :param objective_kwargs: Dictionnary containaing the fixed keyword arguments `
        to the custom `objective_hyperopt.
    :param model_class: scVI model class (e.g ``VAE``, ``VAEC``, ``SCANVI``)
    :param trainer_class: ``Trainer`` sub-class (e.g ``UnsupervisedTrainer``)
    :param metric_name: Name of the metric to optimize for. If `None` defaults to "marginal_ll"
    :param metric_kwargs: keyword arguments for the metric method.
        If `metric_name` is None, defaults to {"n_mc_samples": 100}.
    :param posterior_name: Name of the posterior distribution to compute the metric with.
    :param model_specific_kwargs: ``dict`` of fixed parameters which will be passed to the model.
    :param trainer_specific_kwargs: ``dict`` of fixed parameters which will be passed to the trainer.
    :param train_func_specific_kwargs: dict of fixed parameters which will be passed to the train method.
    :param space: dict containing up to three sub-dicts with keys "model_tunable_kwargs",
        "trainer_tunable_kwargs" or "train_func_tunable_kwargs".
        Each of those dict contains ``hyperopt`` defined parameter spaces (e.g. ``hp.choice(..)``)
        which will be passed to the corresponding object : model, trainer or train method
        when performing hyper-optimization. Default: mutable, see source code.
    :param max_evals: Maximum number of evaluations of the objective.
    :param train_best: If ``True``, train best model and return it.
    :param pickle_result: If ``True``, pickle ``Trials`` and  ``Trainer`` objects using ``save_path``.
    :param save_path: Path where to save best model, trainer, trials and mongo files.
    :param use_batches: If ``False``, pass ``n_batch=0`` to model else pass ``gene_dataset.n_batches``.
    :param parallel: If ``True``, use ``MongoTrials`` object to run trainings in parallel.
    :param n_cpu_workers: Number of cpu workers to launch. If None, and no GPUs are found,
        defaults to ``os.cpucount() - 1``. Else, defaults to 0.
    :param gpu_ids: Ids of the GPUs to use. If None defaults to all GPUs found by ``torch``.
        Note that considered gpu ids are int from 0 to ``torch.cuda.device_count()``.
    :param n_workers_per_gpu: Number of workers to launch per gpu found by ``torch``.
    :param reserve_timeout: Amount of time, in seconds, a worker tries to reserve a job for
        before throwing a ``ReserveTimeout`` Exception.
    :param fmin_timeout: Amount of time, in seconds, fmin_process has to terminate
        after all workers have died - before throwing a ``FminTimeoutError``.
        If ``multiple_hosts`` is set to ``True``, this is set to ``None`` to prevent timing out.
    :param fmin_timer: Global amount of time allowed for fmin_process.
        If not None, the minimization procedure will be stopped after ``fmin_timer`` seconds.
        Used only if ``parallel`` is set to ``True``.
    :param mongo_port: Port to the Mongo db.
    :param mongo_host: Hostname used with ``mongo_port`` to indicate the prefix of the mongodb address.
        The prefix of the address passed onto the workers and ``MongoTrials`` object
        is ``'{mongo_host}:{mongo_port}'``.
    :param db_name: Name to use when creating the Mongo database. Suffix of the Mongo address.
    :param multiple_hosts: If ``True``, user is considered to have workers launched on several machines.
        Therefore, setting this to ``True`` disables the ``fmin_timeout`` behaviour.
    :return: ``Trainer`` object for the best model and ``(Mongo)Trials`` object containing logs for the different runs.

    Examples:
        >>> from scvi.dataset import CortexDataset
        >>> gene_dataset = CortexDataset()
        >>> best_trainer, trials = auto_tune_scvi_model("cortex", gene_dataset)
    """
    global fh_autotune

    # add file handler
    fh_autotune = logging.handlers.RotatingFileHandler(
        os.path.join(save_path, "scvi_autotune_logfile.txt")
    )
    fh_autotune.setFormatter(autotune_formatter)
    fh_autotune.setLevel(logging.DEBUG)
    logger_all.addHandler(fh_autotune)

    if delayed_populating and not isinstance(gene_dataset, DownloadableDataset):
        raise ValueError(
            "The delayed_population mechanism requires an "
            "instance of scvi.dataset.dataset.DownloadableDataset."
        )

    if fmin_timer and train_best:
        logger_all.warning(
            "fmin_timer and train_best are both set to True. "
            "This means that runtime will exceed fmin_timer "
            "by at least the time it takes to complete a full training."
        )

    logger_all.info("Starting experiment: {exp_key}".format(exp_key=exp_key))

    # default search space
    if space is None:
        logger_all.debug("Using default parameter search space.")
        space = {
            "model_tunable_kwargs": {
                "n_latent": 5 + hp.randint("n_latent", 11),  # [5, 15]
                "n_hidden": hp.choice("n_hidden", [64, 128, 256]),
                "n_layers": 1 + hp.randint("n_layers", 5),
                "dropout_rate": hp.choice("dropout_rate", [0.1, 0.3, 0.5, 0.7]),
                "reconstruction_loss": hp.choice("reconstruction_loss", ["zinb", "nb"]),
            },
            "train_func_tunable_kwargs": {
                "lr": hp.choice("lr", [0.01, 0.005, 0.001, 0.0005, 0.0001])
            },
        }

    # default metric
    if metric_name is None:
        metric_name = "marginal_ll"
        metric_kwargs = {"n_mc_samples": 100}

    # build a partial objective function restricted to the search space
    if custom_objective_hyperopt is None:
        # default specific kwargs
        model_specific_kwargs = model_specific_kwargs if model_specific_kwargs else {}
        trainer_specific_kwargs = (
            trainer_specific_kwargs if trainer_specific_kwargs else {}
        )
        train_func_specific_kwargs = (
            train_func_specific_kwargs if train_func_specific_kwargs else {}
        )

        # default early stopping
        if "early_stopping_kwargs" not in trainer_specific_kwargs:
            logger_all.debug("Adding default early stopping behaviour.")
            early_stopping_kwargs = {
                "early_stopping_metric": "elbo",
                "save_best_state_metric": "elbo",
                "patience": 50,
                "threshold": 0,
                "reduce_lr_on_plateau": True,
                "lr_patience": 25,
                "lr_factor": 0.2,
            }
            trainer_specific_kwargs["early_stopping_kwargs"] = early_stopping_kwargs
            # add elbo to metrics to monitor
            metrics_to_monitor = trainer_specific_kwargs.get("metrics_to_monitor", [])
            metrics_to_monitor.append("elbo")
            trainer_specific_kwargs["metrics_to_monitor"] = metrics_to_monitor

        logger_all.info(
            "Fixed parameters: \n"
            "model: \n"
            + str(model_specific_kwargs)
            + "\n"
            + "trainer: \n"
            + str(trainer_specific_kwargs)
            + "\n"
            + "train method: \n"
            + str(train_func_specific_kwargs)
        )
        objective_hyperopt = partial(
            _objective_function,
            **{
                "gene_dataset": gene_dataset,
                "delayed_populating": delayed_populating,
                "model_class": model_class,
                "trainer_class": trainer_class,
                "metric_name": metric_name,
                "metric_kwargs": metric_kwargs,
                "posterior_name": posterior_name,
                "model_specific_kwargs": model_specific_kwargs,
                "trainer_specific_kwargs": trainer_specific_kwargs,
                "train_func_specific_kwargs": train_func_specific_kwargs,
                "use_batches": use_batches,
            },
        )
    else:
        logger_all.info("Using custom objective function.")
        objective_hyperopt = partial(custom_objective_hyperopt, **objective_kwargs)

    if parallel:
        logger_all.info("Starting parallel hyperoptimization")
        trials = _auto_tune_parallel(
            objective_hyperopt=objective_hyperopt,
            exp_key=exp_key,
            space=space,
            max_evals=max_evals,
            save_path=save_path,
            n_cpu_workers=n_cpu_workers,
            gpu_ids=gpu_ids,
            n_workers_per_gpu=n_workers_per_gpu,
            reserve_timeout=reserve_timeout,
            fmin_timeout=fmin_timeout,
            fmin_timer=fmin_timer,
            mongo_port=mongo_port,
            mongo_host=mongo_host,
            db_name=db_name,
            multiple_hosts=multiple_hosts,
        )

    else:
        logger_all.info("Starting sequential hyperoptimization")
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
    if train_best:
        logger_all.debug("Training best model with full training set")
        best_space = trials.best_trial["result"]["space"]
        best_trainer = objective_hyperopt(best_space, is_best_training=True)

    if pickle_result:
        if train_best:
            logger_all.debug("Pickling best model and trainer")
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
        logger_all.debug("Pickling Trials object")
        if hasattr(trials, "handle"):
            del trials.handle
        with open(
            os.path.join(save_path, "trials_{key}".format(key=exp_key)), "wb"
        ) as f:
            pickle.dump(trials, f)

    # remove added logging handlers
    _cleanup_logger()

    if train_best:
        return best_trainer, trials
    else:
        return trials


def _auto_tune_parallel(
    objective_hyperopt: Callable,
    exp_key: str,
    space: dict = None,
    max_evals: int = 100,
    save_path: str = ".",
    n_cpu_workers: int = None,
    gpu_ids: List[int] = None,
    n_workers_per_gpu: int = 1,
    reserve_timeout: float = 180.0,
    fmin_timeout: float = 300.0,
    fmin_timer: float = None,
    mongo_port: str = "1234",
    mongo_host: str = "localhost",
    db_name: str = "scvi_db",
    multiple_hosts: bool = False,
) -> MongoTrials:
    """Parallel version of the hyperoptimization procedure.
    Called by ``auto_tune_scvi_model`` when ``parallel=True``.
    Specifically, first the MongoDb service is launched in its own forked process.
    Then, the call to the minimization process is made in its own forked process.
    Then, the call ``worker_launcher`` is made in its own Thread.
    After that, the program waits for either the minimization
    process to finish or for the workers to all timeout.
    When one of these conditions is verified the program kills the waiter for the other
    and tries to dequeue the results from the minimization process.
    At that point, if ``multiple_hosts`` is set to True, the program waits indefinitely
    for the minimization process to put the results in the queue.
    If not, the minimisation process has ``fmin_timeout`` seconds to finish.
    This mechanism ensures that the program does not hang if, for any reason,
    the workers die before completing all the jobs.
    Note that logs to the ``hyperopt`` package are automatically stored in ``./hyperopt_logfile.txt``.
    Note that the progress bar is automatically disabled if the logging level
    for ``scvi.inference.autotune`` is lower than logging.WARNING.

    :param objective_hyperopt: Callable, the objective function to minimize.
    :param exp_key: Name of the experiment in MongoDb.
    :param space: ``dict`` containing up to three sub-dicts with keys "model_tunable_kwargs",
        "trainer_tunable_kwargs" or "train_func_tunable_kwargs".
        Each of those dict contains ``hyperopt`` defined parameter spaces (e.g. ``hp.choice(..)``)
        which will be passed to the corresponding object : model, trainer or train method
        when performing hyperoptimization. Default: mutable, see source code.
    :param max_evals: Maximum number of evaluations of the objective.
    :param save_path: Path where to save best model, trainer, trials and mongo files.
    :param n_cpu_workers: Number of cpu workers to launch. If None, and no GPUs are found,
        defaults to ``os.cpucount() - 1``. Else, defaults to 0.
    :param gpu_ids: Ids of the GPUs to use. If None defaults to all GPUs found by ``torch``.
        Note that considered gpu ids are int from ``0`` to ``torch.cuda.device_count()``.
    :param n_workers_per_gpu: Number of workers ton launch per gpu found by ``torch``.
    :param reserve_timeout: Amount of time, in seconds, a worker tries to reserve a job for
        before throwing a ``ReserveTimeout`` Exception.
    :param fmin_timeout: Amount of time, in seconds, ``fmin_process`` has to terminate
        after all workers have died - before throwing a ``FminTimeoutError``.
        If ``multiple_hosts`` is set to ``True``, this is set to None to disable the timeout behaviour.
    :param fmin_timer: Global amount of time allowed for fmin_process.
        If not None, the minimization procedure will be stopped after ``fmin_timer`` seconds.
        Used only if ``parallel`` is set to ``True``.
    :param mongo_port: Port to the mongo db.
    :param mongo_host: Hostname used with mongo_port to indicate the prefix of the mongodb address.
        The prefix of the address passed onto the workers and MongoTrials object is ``'{mongo_host}:{mongo_port}'``.
    :param db_name: Name to use when creating the Mongo database. Suffix of the mongo address.
    :param multiple_hosts: If ``True``, user is considered to have workers launched on several machines.
        Therefore, setting this to ``True`` disables the ``fmin_timeout`` behaviour.
    :return: ``MongoTrials`` object containing the results of the program.
    """
    global started_processes
    global started_threads
    global started_queues
    global fh_hyperopt

    # prepare parallel logging
    logging_queue = multiprocessing.Queue()
    started_queues.append(logging_queue)
    listener = QueueListener(logging_queue, DispatchHandler())
    listener.start()
    started_processes.append(listener)

    # run mongod bash script
    mongo_path = os.path.join(save_path, "mongo")
    if not os.path.exists(mongo_path):
        os.makedirs(mongo_path)
    mongo_logfile = open(os.path.join(mongo_path, "mongo_logfile.txt"), "w")
    open_files.append(mongo_logfile)
    logger_all.debug(
        "Starting MongoDb process, logs redirected to "
        "{name}.".format(name=mongo_logfile.name)
    )
    mongod_process = Popen(
        [
            "mongod",
            "--quiet",
            "--dbpath={path}".format(path=mongo_path),
            "--port={port}".format(port=mongo_port),
        ],
        stdout=mongo_logfile,
    )
    # let mongo server start and check it did
    time.sleep(5)
    client = pymongo.MongoClient(
        mongo_host + ":" + mongo_port, serverSelectionTimeoutMS=100
    )
    try:
        client.server_info()
        client.close()
    except pymongo.mongo_client.ServerSelectionTimeoutError:
        logger_all.error("Failed to connect to mongo agent.")
        mongo_logfile.close()
        mongo_logfile = open(os.path.join(mongo_path, "mongo_logfile.txt"), "r")
        logger_all.error(
            "Logs for the mongod subprocess: \n" + "".join(mongo_logfile.readlines())
        )
        raise

    mongo_url = os.path.join(mongo_host + ":" + mongo_port, db_name)
    started_processes.append(mongod_process)

    # log hyperopt only to file
    hp_logger = logging.getLogger("hyperopt")
    hp_logger.propagate = False
    fh_hyperopt = logging.handlers.RotatingFileHandler(
        os.path.join(save_path, "hyperopt_logfile.txt")
    )
    fh_hyperopt.setFormatter(autotune_formatter)
    hp_logger.addHandler(fh_hyperopt)

    # start fmin launcher thread
    logger_all.debug("Starting minimization procedure")
    queue = multiprocessing.Queue()
    started_queues.append(queue)
    fmin_launcher_thread = FminLauncherThread(
        logging_queue=logging_queue,
        queue=queue,
        objective_hyperopt=objective_hyperopt,
        exp_key=exp_key,
        space=space,
        algo=tpe.suggest,
        max_evals=max_evals,
        fmin_timer=fmin_timer,
        mongo_url=mongo_url,
    )
    fmin_launcher_thread.start()
    started_threads.append(fmin_launcher_thread)

    # start worker launcher
    logger_all.debug("Starting worker launcher")
    worker_launcher_thread = WorkerLauncherThread(
        logging_queue=logging_queue,
        exp_key=exp_key,
        n_cpu_workers=n_cpu_workers,
        gpu_ids=gpu_ids,
        n_workers_per_gpu=n_workers_per_gpu,
        reserve_timeout=reserve_timeout,
        workdir=mongo_path,
        mongo_url=mongo_url,
        multiple_hosts=multiple_hosts,
        max_evals=max_evals,
    )
    worker_launcher_thread.start()
    started_threads.append(worker_launcher_thread)

    # wait for one to finish
    while worker_launcher_thread.is_alive() and fmin_launcher_thread.is_alive():
        time.sleep(5)

    if not fmin_launcher_thread.is_alive():
        logger_all.debug("Setting worker launcher stop event.")
        worker_launcher_thread.stop_event.set()
    try:
        if multiple_hosts:
            # if using multiple_hosts, there could still be workers -> disable fmin timeout
            fmin_timeout = None
            logger_all.debug(
                "multiple_hosts set to True, fmin will block until all trials have been completed."
            )
        else:
            logger_all.debug(
                "multiple_hosts set to false, Fmin has {time} seconds to finish".format(
                    time=fmin_timeout
                )
            )
        trials = queue.get(timeout=fmin_timeout)
        queue.close()
    except Empty:
        logger_all.error(
            "Queue still empty {fmin_timeout} seconds after all workers have died."
            "\n".format(fmin_timeout=fmin_timeout) + "Terminating minimization process."
        )
        raise FminTimeoutError(
            "Queue still empty {fmin_timeout} seconds after all workers "
            "have died. Check that you have used a new exp_key or allowed "
            "a higher max_evals".format(fmin_timeout=fmin_timeout)
        )

    # sanity: wait for fmin, terminate workers and wait for launcher
    fmin_launcher_thread.join()
    worker_launcher_thread.join()
    logger_all.info(
        "Finished minimization procedure for experiment {exp_key}.".format(
            exp_key=exp_key
        )
    )
    logger_all.debug("Terminating mongod process.")
    mongod_process.terminate()
    # wait for process to actually terminate, avoid issues with unreleased mongod.lock
    mongod_process.wait()
    mongo_logfile.close()
    mongo_logfile = open(os.path.join(mongo_path, "mongo_logfile.txt"), "r")
    logger_all.info(
        "Logs for the mongod subprocess: \n" + "".join(mongo_logfile.readlines())
    )
    logger_all.debug("Stopping asynchronous logging listener.")
    listener.stop()
    logging_queue.close()

    # cleanup queues, processes, threads, files and logger
    _cleanup_processes_files()
    _cleanup_logger()

    return trials


class FminLauncherThread(StoppableThread):
    """Starts the process which ultimately call the minimzation procedure.

    Is encapsulated in a ``threading.Thread`` to allow for the ``fmin_timer`` mechanism.

    :param logging_queue: Queue to send logs to main process using a ``QueueHandler``.
        Here to be passed on to `FminProcess`.
    :param queue: Queue to put trials in. Here to be passed on to `FminProcess`.
    :param objective_hyperopt: Callable, the objective function to minimize
    :param exp_key: Name of the experiment in MongoDb.
    :param space: ``dict`` containing up to three sub-dicts with keys "model_tunable_kwargs",
        "trainer_tunable_kwargs" or "train_func_tunable_kwargs".
        Each of those dict contains ``hyperopt`` defined parameter spaces (e.g. ``hp.choice(..)``)
        which will be passed to the corresponding object : model, trainer or train method
        when performing hyperoptimization. Default: mutable, see source code.
    :param algo: Bayesian optimization algorithm from ``hyperopt`` to use.
    :param max_evals: Maximum number of evaluations of the objective.
    :param fmin_timer: Global amount of time allowed for fmin_process.
        If not None, the minimization procedure will be stopped after ``fmin_timer`` seconds.
        Used only if ``parallel`` is set to ``True``.
    :param mongo_url: String of the form mongo_host:mongo_port/db_name.
    """

    def __init__(
        self,
        logging_queue: multiprocessing.Queue,
        queue: multiprocessing.Queue,
        objective_hyperopt: Callable,
        exp_key: str,
        space: dict,
        algo: Callable = tpe.suggest,
        max_evals: int = 100,
        fmin_timer: float = None,
        mongo_url: str = "localhost:1234/scvi_db",
    ):
        super().__init__(name="Fmin Launcher")
        self.logging_queue = logging_queue
        self.queue = queue
        self.objective_hyperopt = objective_hyperopt
        self.exp_key = exp_key
        self.space = space
        self.algo = algo
        self.max_evals = max_evals
        self.fmin_timer = fmin_timer
        self.mongo_url = mongo_url

    @_error_logger_decorator
    def run(self):
        """Launches a ``hyperopt`` minimization procedure."""
        # call fmin in a process to enable termination
        fmin_process = FminProcess(
            logging_queue=self.logging_queue,
            queue=self.queue,
            objective_hyperopt=self.objective_hyperopt,
            space=self.space,
            mongo_url=self.mongo_url,
            exp_key=self.exp_key,
            algo=self.algo,
            max_evals=self.max_evals,
        )
        logger_all.debug("Starting FminProcess.")
        fmin_process.start()
        started_processes.append(fmin_process)
        if self.fmin_timer is not None:
            logger_all.info(
                "Timer set, fmin will run for at most {timer}.".format(
                    timer=self.fmin_timer
                )
            )
            start_time = time.monotonic()
            run_time = 0
            while (
                run_time < self.fmin_timer
                and fmin_process.is_alive()
                and not self.stop_event.is_set()
            ):
                time.sleep(10)
                run_time = time.monotonic() - start_time
            if self.stop_event.is_set():
                logger_all.debug("Stop event set.")
            elif run_time > self.fmin_timer and fmin_process.is_alive():
                logger_all.debug(
                    "Timer ran out. Terminating FminProcess and putting current Trials in queue."
                )
                fmin_process.terminate()
                # queue.put uses pickle so remove attribute containing thread.lock
                trials = MongoTrials(
                    as_mongo_str(os.path.join(self.mongo_url, "jobs")),
                    exp_key=self.exp_key,
                )
                if hasattr(trials, "handle"):
                    logger_all.debug("Deleting Trial handle for pickling.")
                    del trials.handle
                logger_all.debug("Putting Trials in Queue.")
                self.queue.put(trials)
            else:
                logger_all.debug("fmin finished.")
        else:
            logger_all.debug("No timer, waiting for fmin...")
            while fmin_process.is_alive() and not self.stop_event.is_set():
                time.sleep(10)
            logger_all.debug("fmin finished.")


class FminProcess(multiprocessing.Process):
    """Call ``hyperopt``'s fmin.

    Is encapsulated in a ``multiprocessing.Process`` in order to
    allow for termination in case cleanup is required.

    :param logging_queue: Queue to send logs to main process using a ``QueueHandler``.
    :param queue: Queue to put trials in.
    :param objective_hyperopt: Callable, the objective function to minimize
    :param space: ``dict`` containing up to three sub-dicts with keys "model_tunable_kwargs",
        "trainer_tunable_kwargs" or "train_func_tunable_kwargs".
        Each of those dict contains ``hyperopt`` defined parameter spaces (e.g. ``hp.choice(..)``)
        which will be passed to the corresponding object : model, trainer or train method
        when performing hyperoptimization. Default: mutable, see source code.
    :param exp_key: Name of the experiment in MongoDb.
    :param mongo_url: String of the form mongo_host:mongo_port/db_name
    :param algo: Bayesian optimization algorithm from ``hyperopt`` to use.
    :param max_evals: Maximum number of evaluations of the objective.
    :param show_progressbar: Whether or not to show the ``hyperopt`` progress bar.
    """

    def __init__(
        self,
        logging_queue: multiprocessing.Queue,
        queue: multiprocessing.Queue,
        objective_hyperopt: Callable,
        space: dict,
        exp_key: str,
        mongo_url: str = "localhost:1234/scvi_db",
        algo: Callable = tpe.suggest,
        max_evals: int = 100,
        show_progressbar: bool = False,
    ):
        super().__init__(name="Fmin")
        self.logging_queue = logging_queue
        self.queue = queue
        self.objective_hyperopt = objective_hyperopt
        self.space = space
        self.mongo_url = mongo_url
        self.exp_key = exp_key
        self.algo = algo
        self.max_evals = max_evals
        self.show_progressbar = show_progressbar

    @_asynchronous_logging_method_decorator
    @_error_logger_decorator
    def run(self):
        logger_all.debug("Instantiating MongoTrials object.")
        trials = MongoTrials(
            as_mongo_str(os.path.join(self.mongo_url, "jobs")), exp_key=self.exp_key
        )
        logger_all.debug("Calling fmin.")
        fmin(
            fn=self.objective_hyperopt,
            space=self.space,
            algo=self.algo,
            max_evals=self.max_evals,
            trials=trials,
            show_progressbar=self.show_progressbar,
        )
        # queue.put uses pickle so remove attribute containing thread.lock
        if hasattr(trials, "handle"):
            logger_all.debug("fmin returned. Deleting Trial handle for pickling.")
            del trials.handle
        logger_all.debug("Putting Trials in Queue.")
        self.queue.put(trials)


class WorkerLauncherThread(StoppableThread):
    """Launches the local workers which are going to run the jobs required by the minimization process.
    Terminates when the worker_watchdog call finishes.
    Specifically, first ``n_gpu_workers`` are launched per GPU in ``gpu_ids`` in their own spawned process.
    Then, ``n_cpu_workers`` CPU workers are launched, also in their own spawned process.
    The use of spawned processes (each have their own python interpreter) is mandatory for compatiblity with CUDA.
    See https://pytorch.org/docs/stable/notes/multiprocessing.html for more information.

    :param logging_queue: Queue to send logs to main process using a ``QueueHandler``.
        Here to be passed on to the `HyperoptWorker` processes.
    :param exp_key: This key is used by hyperopt as a suffix to the part of the MongoDb
        which corresponds to the current experiment. In particular, it has to be passed to ``MongoWorker``.
    :param n_cpu_workers: Number of cpu workers to launch. If None, and no GPUs are found,
        defaults to ``os.cpu_count() - 1``. Else, defaults to 0.
    :param gpu_ids: Ids of the GPUs to use. If None defaults to all GPUs found by ``torch``.
        Note that considered gpu ids are int from ``0`` to ``torch.cuda.device_count()``.
    :param n_workers_per_gpu: Number of workers ton launch per gpu found by ``torch``.
    :param reserve_timeout: Amount of time, in seconds, a worker tries to reserve a job for
        before throwing a ``ReserveTimeout`` Exception.
    :param workdir: Directory where the workers
    :param mongo_url: Address to the running MongoDb service.
    :param multiple_hosts: ``True`` if launching workers form multiple hosts.
    :param max_evals: Maximum number of evaluations of the objective.
        Useful for instantiating a progress bar.
    """

    def __init__(
        self,
        logging_queue: multiprocessing.Queue,
        exp_key: str,
        n_cpu_workers: int = None,
        gpu_ids: List[int] = None,
        n_workers_per_gpu: int = 1,
        reserve_timeout: float = 30.0,
        workdir: str = ".",
        mongo_url: str = "localhost:1234/scvi_db",
        multiple_hosts: bool = False,
        max_evals: int = 100,
    ):
        super().__init__(name="Worker Launcher")
        self.logging_queue = logging_queue
        self.exp_key = exp_key
        self.n_cpu_workers = n_cpu_workers
        self.gpu_ids = gpu_ids
        self.n_workers_per_gpu = n_workers_per_gpu
        self.reserve_timeout = reserve_timeout
        self.workdir = workdir
        self.mongo_url = mongo_url
        self.multiple_hosts = multiple_hosts
        self.max_evals = max_evals

    @_error_logger_decorator
    def run(self):
        global started_processes

        if self.gpu_ids is None:
            n_gpus = torch.cuda.device_count()
            logger_all.debug(
                "gpu_ids is None, defaulting to all {n_gpus} GPUs found by torch.".format(
                    n_gpus=n_gpus
                )
            )
            self.gpu_ids = list(range(n_gpus))
            if n_gpus and self.n_cpu_workers is None:
                self.n_cpu_workers = 0
                logger_all.debug(
                    "Some GPU.s found and n_cpu_wokers is None, defaulting to n_cpu_workers = 0"
                )
            if not n_gpus and self.n_cpu_workers is None:
                self.n_cpu_workers = os.cpu_count() - 1
                logger_all.debug(
                    "No GPUs found and n_cpu_wokers is None, defaulting to n_cpu_workers = "
                    "{n_cpu_workers} (os.cpu_count() - 1)".format(
                        n_cpu_workers=self.n_cpu_workers
                    )
                )
        if (
            self.gpu_ids is None
            and (self.n_cpu_workers == 0 or self.n_cpu_workers is None)
            and not self.multiple_hosts
        ):
            raise ValueError("No hardware (cpu/gpu) selected/found.")

        # log progress with queue and progress_listener
        pbar = None
        if not self.multiple_hosts and logger.level >= logging.WARNING:
            pbar = tqdm.tqdm(total=self.max_evals)
        progress_queue = multiprocessing.Queue()
        started_queues.append(progress_queue)
        prog_listener = ProgressListener(progress_queue=progress_queue, pbar=pbar)
        prog_listener.start()
        started_threads.append(prog_listener)

        running_workers = []
        # launch gpu workers
        logger_all.info(
            "Starting {n_workers_per_gpu} worker.s for each of the {n_gpus} gpu.s set for use/"
            "found.".format(
                n_workers_per_gpu=self.n_workers_per_gpu, n_gpus=len(self.gpu_ids)
            )
        )
        for gpu_id in self.gpu_ids:
            for sub_id in range(self.n_workers_per_gpu):
                worker = HyperoptWorker(
                    progress_queue=progress_queue,
                    logging_queue=self.logging_queue,
                    exp_key=self.exp_key,
                    workdir=self.workdir,
                    gpu=True,
                    hw_id=str(gpu_id),
                    reserve_timeout=self.reserve_timeout,
                    mongo_url=self.mongo_url,
                    name="Worker GPU " + str(gpu_id) + ":" + str(sub_id),
                )
                worker.start()
                running_workers.append(worker)

        # launch cpu workers
        logger_all.info(
            "Starting {n_cpu_workers} cpu worker.s".format(
                n_cpu_workers=self.n_cpu_workers
            )
        )
        for cpu_id in range(self.n_cpu_workers):
            worker = HyperoptWorker(
                progress_queue=progress_queue,
                logging_queue=self.logging_queue,
                exp_key=self.exp_key,
                workdir=self.workdir,
                gpu=False,
                hw_id=str(cpu_id),
                reserve_timeout=self.reserve_timeout,
                mongo_url=self.mongo_url,
                name="Worker CPU " + str(cpu_id),
            )
            worker.start()
            running_workers.append(worker)
        started_processes.extend(running_workers)

        # wait or return if all workers have died
        while not self.stop_event.is_set():
            n_alive = 0
            for worker in running_workers:
                n_alive += 1 if worker.is_alive() else n_alive
            if n_alive == 0:
                logger_all.debug(
                    "All workers have died, check stdout/stderr for error tracebacks."
                )
                break
        logger_all.debug(
            "Worker watchdog finished, terminating workers and stopping listener."
        )
        for worker in running_workers:
            if worker.is_alive():
                worker.terminate()
        prog_listener.stop_event.set()
        prog_listener.join()


class ProgressListener(StoppableThread):
    """Listens to workers when they finish a job and logs progress.

    Workers put in the progress_queue when they finish a job
    and when they do this function sends a log to the progress logger.
    """

    def __init__(self, progress_queue: multiprocessing.Queue, pbar: tqdm.tqdm = None):
        super().__init__(name="Progress Listener")
        self.progress_queue = progress_queue
        self.pbar = pbar

    @_error_logger_decorator
    def run(self):
        logger_all.debug("Listener listening...")

        i = 0
        while not self.stop_event.is_set():
            # get job done signal
            try:
                self.progress_queue.get(block=False)
                i += 1
                logger_all.info("{i} job.s done".format(i=i))
                # update progress bar through ProgressHandler
                if self.pbar is not None:
                    self.pbar.update()
            except Empty:
                pass
            time.sleep(5)
        if self.pbar is not None:
            self.pbar.close()
        self.progress_queue.close()


class HyperoptWorker(multiprocessing.Process):
    """Launches a ``hyperopt`` ``MongoWorker`` which runs jobs until ``ReserveTimeout`` is raised.

    :param progress_queue: Queue in which to put None when a job is done.
    :param logging_queue: Queue to send logs to main process using a ``QueueHandler``.
    :param exp_key: This key is used by hyperopt as a suffix to the part of the MongoDb
        which corresponds to the current experiment. In particular, it has to be passed to ``MongoWorker``.
    :param workdir:
    :param gpu: If ``True`` means a GPU is to be used.
    :param hw_id: Id of the GPU to use. set via env variable ``CUDA_VISIBLE_DEVICES``.
    :param poll_interval: Time to wait between attempts to reserve a job.
    :param reserve_timeout: Amount of time, in seconds, a worker tries to reserve a job for
        before throwing a ``ReserveTimeout`` Exception.
    :param mongo_url: Address to the running MongoDb service.
    """

    def __init__(
        self,
        name: str,
        progress_queue: multiprocessing.Queue,
        logging_queue: multiprocessing.Queue,
        exp_key: str,
        workdir: str = ".",
        gpu: bool = True,
        hw_id: str = None,
        poll_interval: float = 1.0,
        reserve_timeout: float = 30.0,
        mongo_url: str = "localhost:1234/scvi_db",
    ):
        super().__init__(name=name)
        self.progress_queue = progress_queue
        self.logging_queue = logging_queue
        self.exp_key = exp_key
        self.workdir = workdir
        self.gpu = gpu
        self.hw_id = hw_id
        self.poll_interval = poll_interval
        self.reserve_timeout = reserve_timeout
        self.mongo_url = mongo_url

    @_asynchronous_logging_method_decorator
    @_error_logger_decorator
    def run(self):
        logger_all.debug("Worker working...")

        os.environ["CUDA_VISIBLE_DEVICES"] = self.hw_id if self.gpu else str()

        mjobs = MongoJobs.new_from_connection_str(
            os.path.join(as_mongo_str(self.mongo_url), "jobs")
        )
        mworker = MongoWorker(
            mjobs, float(self.poll_interval), workdir=self.workdir, exp_key=self.exp_key
        )

        while True:
            try:
                mworker.run_one(reserve_timeout=float(self.reserve_timeout))
                self.progress_queue.put(None)
            except ReserveTimeout:
                logger_all.debug(
                    "Caught ReserveTimeout. "
                    "Exiting after failing to reserve job for {time} seconds.".format(
                        time=self.reserve_timeout
                    )
                )
                break


@_error_logger_decorator
def _objective_function(
    space: dict,
    gene_dataset: GeneExpressionDataset,
    delayed_populating: bool = False,
    model_class: Type[VAE] = VAE,
    trainer_class: Type[Trainer] = UnsupervisedTrainer,
    metric_name: str = None,
    metric_kwargs: Dict[str, Any] = None,
    posterior_name: str = "test_set",
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
    :param metric_name: Name of the metric to optimize for. If `None` defaults to "marginal_ll"
    :param metric_kwargs: keyword arguments for the metric method.
        If `metric_name` is None, defaults to {"n_mc_samples": 100}.
    :param posterior_name: Name of the posterior distribution to compute the metric with.
    :param model_specific_kwargs: ``dict`` of fixed parameters which will be passed to the model.
    :param model_specific_kwargs: dict of fixed parameters which will be passed to the model.
    :param trainer_specific_kwargs: dict of fixed parameters which will be passed to the trainer.
    :param train_func_specific_kwargs: dict of fixed parameters which will be passed to the train method.
    :param use_batches: If False, pass n_batch=0 to model else pass gene_dataset.n_batches
    :param is_best_training: True if training the model with the best hyperparameters
    :return: best value of the early stopping metric, and best model if is_best_training
    """
    # handle mutable defaults
    metric_kwargs = metric_kwargs if metric_kwargs is not None else {}

    if delayed_populating and isinstance(gene_dataset, DownloadableDataset):
        gene_dataset.populate()

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
        logger_all.info(
            "Parameters being tested: \n"
            "model: \n"
            + str(model_tunable_kwargs)
            + "\n"
            + "trainer: \n"
            + str(trainer_tunable_kwargs)
            + "\n"
            + "train method: \n"
            + str(train_func_tunable_kwargs)
        )

    # define model
    logger_all.debug("Instantiating model")
    model = model_class(
        n_input=gene_dataset.nb_genes,
        n_batch=gene_dataset.n_batches * use_batches,
        **model_tunable_kwargs,
    )

    # define trainer
    logger_all.debug("Instantiating trainer")
    trainer = trainer_class(model, gene_dataset, **trainer_tunable_kwargs)

    # train model
    logger_all.debug("Starting training")
    trainer.train(**train_func_tunable_kwargs)
    logger_all.debug("Finished training")
    elapsed_time = time.monotonic() - start_time
    # if training the best model, return model else return criterion
    if is_best_training:
        return trainer
    else:
        # select metric from early stopping kwargs if possible
        metric = None
        save_best_state_metric = None
        early_stopping_kwargs = trainer_specific_kwargs.get(
            "early_stopping_kwargs", None
        )
        if early_stopping_kwargs is not None:
            metric = early_stopping_kwargs.get("early_stopping_metric", None)
            save_best_state_metric = early_stopping_kwargs.get(
                "save_best_state_metric", None
            )

        # store run results
        if metric is not None:
            early_stopping_loss_is_best = True
            best_epoch = trainer.best_epoch
            # add actual number of epochs to be used when training best model
            space["train_func_tunable_kwargs"]["n_epochs"] = best_epoch
            early_stopping_loss = trainer.early_stopping.best_performance
            metric += "_" + trainer.early_stopping.on
        # default to elbo
        else:
            early_stopping_loss_is_best = False
            metric = "elbo_test_set"
            early_stopping_loss = trainer.history[metric][-1]
            best_epoch = len(trainer.history[metric])

        # load best state
        if save_best_state_metric is not None:
            model.load_state_dict(trainer.best_state_dict)

        # compute optimized metric
        loss = getattr(getattr(trainer, posterior_name), metric_name)(**metric_kwargs)

        logger_all.debug(
            "Training of {n_epochs} epochs finished in {time} with loss = {loss}".format(
                n_epochs=len(trainer.history[metric]),
                time=str(datetime.timedelta(seconds=elapsed_time)),
                loss=loss,
            )
        )

        # check status
        status = STATUS_OK
        if np.isnan(loss):
            status = STATUS_FAIL

        return {
            "loss": loss,
            "early_stopping_loss": early_stopping_loss,
            "early_stopping_loss_is_best": early_stopping_loss_is_best,
            "best_epoch": best_epoch,
            "elapsed_time": elapsed_time,
            "status": status,
            "history": trainer.history,
            "space": space,
            "worker_name": multiprocessing.current_process().name,
        }
