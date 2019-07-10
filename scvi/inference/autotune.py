import datetime
import logging
import multiprocessing
import os
import pickle
import sys
import time
import threading

from collections import defaultdict
from functools import partial, wraps
from io import IOBase
from logging.handlers import QueueListener, QueueHandler
from subprocess import Popen
from typing import Any, Callable, Dict, List, Type, Union
from queue import Empty

from hyperopt import fmin, tpe, Trials, hp, STATUS_OK, STATUS_FAIL
from hyperopt.mongoexp import (
    as_mongo_str,
    MongoJobs,
    MongoTrials,
    MongoWorker,
    ReserveTimeout,
)

import numpy as np
import torch
import tqdm

from scvi.dataset import GeneExpressionDataset
from scvi.models import VAE

from . import Trainer, UnsupervisedTrainer

# TODO: add database watcher and visualizations
# TODO: make worker_launcher a subclass of threading.Thread
# TODO: and hyperopt_worker a subclass of multiprocessing.Process

# spawning is required for processes relying on cuda
spawn_ctx = multiprocessing.get_context("spawn")
fork_ctx = multiprocessing.get_context("fork")

# register running process and open files to terminate/close at exit
started_processes: List[Union[multiprocessing.Process, Popen, QueueListener]] = []
started_threads: List[threading.Thread] = []
open_files: List[IOBase] = []

# instantiate logger, handler and formatter
logger = logging.getLogger(__name__)
formatter = logging.Formatter(
    "[%(asctime)s - %(processName)s - %(threadName)s] %(levelname)s - %(name)s\n%(message)s"
)
ch = logging.StreamHandler()
ch.setFormatter(formatter)
# instantiate hyperopt and autotune file handlers as global variables for clean up
fh_hyperopt = None
fh_autotune = None

# global Event to stop threads when cleaning up
cleanup_event = threading.Event()


class FminTimeoutError(Exception):
    """Thrown if fmin process hasn't finished in the allotted
    time after all workers have died.
    """


class DispatchHandler:
    """A simple handler for logging events. It dispatches events to loggers
    based on the name in the received record, which then get dispatched,
    by the logging system, to the handlers, configured for those loggers.
    """

    def handle(self, record: logging.LogRecord):
        logger = logging.getLogger(record.name)
        if record.levelno >= logger.level:
            logger.handle(record)


class ProgressHandler:
    """A simple handler for keeping track of the worker's progress.
    When assigned to a logger, logs sent using that logger trigger
    an update of the progress bar associated with this handler.
    """

    def __init__(self, pbar: tqdm.tqdm, disable: bool):
        self.level = 0
        self.pbar = pbar
        self.disabled = disable

    def handle(self, record: logging.LogRecord):
        if not self.disabled:
            self.pbar.update()


# cleanup helpers
def _cleanup_processes_files():
    """Cleanup function, starts with latest processes/files.
    Terminates processes, sets cleanup_event to stop threads, closes open files."""
    logger.info("Cleaning up")
    logger.debug("Cleaning up: closing files")
    for f in open_files[::-1]:
        if not f.closed:
            f.close()
    logger.debug("Cleaning up: setting cleanup_event and joining threads")
    cleanup_event.is_set()
    for t in started_threads[::-1]:
        if t.is_alive():
            t.join()
    logger.debug("Cleaning up: terminating processes")
    for p in started_processes[::-1]:
        if isinstance(p, Popen):
            if p.poll() is not None:
                p.terminate()
        if isinstance(p, multiprocessing.Process):
            if p.is_alive():
                p.terminate()
        if isinstance(p, QueueListener):
            if p._thread is not None:
                p.stop()


def _cleanup_logger():
    """Removes added handlers."""
    logger.debug("Cleaning up: removing added logging handler")
    for handler in logger.handlers:
        if handler == ch:
            logger.removeHandler(ch)
    for handler in logging.getLogger("hyperopt").handlers:
        if handler == fh_hyperopt:
            logger.removeHandler(fh_hyperopt)
        if handler == fh_autotune:
            logger.removeHandler(fh_autotune)


def _cleanup_decorator(func: Callable):
    """Decorates top-level calls in order to launch cleanup when an Exception is caught."""

    @wraps(func)
    def decorated(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logger.exception(
                "Caught {exception} in {func}, starting cleanup".format(
                    exception=e.args, func=func.__name__
                )
            )
            _cleanup_processes_files()
            _cleanup_logger()
            raise

    return decorated


def auto_tune_scvi_model(
    exp_key: str,
    gene_dataset: GeneExpressionDataset,
    objective_hyperopt: Callable = None,
    model_class: VAE = VAE,
    trainer_class: Trainer = UnsupervisedTrainer,
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
    reserve_timeout: float = 30.0,
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
    :param gene_dataset: scVI gene dataset.
    :param objective_hyperopt: A custom objective function respecting the ``hyperopt`` format.
        Roughly, it needs to return the quantity to optimize for, either directly
        or in a ``dict`` under the "loss" key.
        See https://github.com/hyperopt/hyperopt/wiki for a more detailed explanation.
        By default, we provide an objective function which can be parametrized
        through the various arguments of this function (``gene_dataset``, ``model_class``, etc.)
    :param model_class: scVI model class (e.g ``VAE``, ``VAEC``, ``SCANVI``)
    :param trainer_class: ``Trainer`` sub-class (e.g ``UnsupervisedTrainer``)
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
        >>> best_trainer, trials = auto_tune_scvi_model(gene_dataset)
    """
    if fmin_timer and train_best:
        logger.warning(
            "fmin_timer and train_best are both set to True. "
            "This means that runtime will exceed fmin_timer "
            "by at least the time it takes to complete a full training."
        )

    # if no handlers add console handler, add formatter to handlers
    if len(logger.handlers) < 1:
        logger.addHandler(ch)
    else:
        # if no formatter add default module formatter
        for handler in logger.handlers:
            if not handler.formatter:
                handler.setFormatter(formatter)

    # also add file handler
    fh_autotune = logging.FileHandler(
        os.path.join(save_path, "scvi_autotune_logfile.txt")
    )
    fh_autotune.setFormatter(formatter)
    fh_autotune.setLevel(logging.DEBUG)
    logger.addHandler(fh_autotune)

    logger.info("Starting experiment: {exp_key}".format(exp_key=exp_key))
    # default specific kwargs
    model_specific_kwargs = model_specific_kwargs if model_specific_kwargs else {}
    trainer_specific_kwargs = trainer_specific_kwargs if trainer_specific_kwargs else {}
    train_func_specific_kwargs = (
        train_func_specific_kwargs if train_func_specific_kwargs else {}
    )

    # default early stopping
    if "early_stopping_kwargs" not in trainer_specific_kwargs:
        logger.debug("Adding default early stopping behaviour.")
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

    # default search space
    if space is None:
        logger.debug("Using default parameter search space.")
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
    )

    # build a partial objective function restricted to the search space
    if objective_hyperopt is None:
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
    if train_best:
        logger.debug("Training best model with full training set")
        best_space = trials.best_trial["result"]["space"]
        best_trainer = objective_hyperopt(best_space, is_best_training=True)

    if pickle_result:
        if train_best:
            logger.debug("Pickling best model and trainer")
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
        logger.debug("Pickling Trials object")
        if hasattr(trials, "handle"):
            del trials.handle
        with open(
            os.path.join(save_path, "trials_{key}".format(key=exp_key)), "wb"
        ) as f:
            pickle.dump(trials, f)

    # remove added logging handlers/formatters
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
    reserve_timeout: float = 30.0,
    fmin_timeout: float = 60.0,
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

    :param objective_hyperopt: Callable, the objective function to minimize
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
        If ``multiple_hosts`` is set to ``True``, this is set to None to disable the timineout behaviour.
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
    # run mongod bash script
    mongo_path = os.path.join(save_path, "mongo")
    if not os.path.exists(mongo_path):
        os.makedirs(mongo_path)
    mongo_logfile = open(os.path.join(mongo_path, "mongo_logfile.txt"), "w")
    open_files.append(mongo_logfile)
    logger.debug(
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
    mongo_port_address = os.path.join(mongo_host + ":" + mongo_port, db_name)
    started_processes.append(mongod_process)

    # log hyperopt to file
    hp_logger = logging.getLogger("hyperopt")
    fh_hyperopt = logging.FileHandler(os.path.join(save_path, "hyperopt_logfile.txt"))
    fh_hyperopt.setFormatter(formatter)
    hp_logger.addHandler(fh_hyperopt)

    # add progress handler to progress logger
    progress_logger = logging.getLogger("progress_logger")
    disable = multiple_hosts or (logger.level < logging.WARNING)
    pbar = tqdm.tqdm(total=max_evals, disable=disable)
    progress_logger.addHandler(ProgressHandler(pbar=pbar, disable=disable))

    # start by running fmin process so that workers don't timeout
    # run hyperoptimization, in a forked process
    # this allows to warn if the workers crash
    # since mongo is not thread-safe, trials must be instantiated in each child
    logger.debug("Starting minimization procedure")
    queue = fork_ctx.Queue()
    fmin_kwargs = {
        "queue": queue,
        "fn": objective_hyperopt,
        "exp_key": exp_key,
        "space": space,
        "algo": tpe.suggest,
        "max_evals": max_evals,
        "fmin_timer": fmin_timer,
        "show_progressbar": False,  # progbar useless in parallel mode
        "mongo_port_address": mongo_port_address,
    }
    fmin_process = fork_ctx.Process(
        target=_fmin_parallel, kwargs=fmin_kwargs, name="fmin Process"
    )
    fmin_process.start()
    started_processes.append(fmin_process)

    # start worker launcher
    logger.debug("Starting worker launcher")
    stop_watchdog_event = threading.Event()
    launcher_kwargs = {
        "stop_watchdog_event": stop_watchdog_event,
        "exp_key": exp_key,
        "n_cpu_workers": n_cpu_workers,
        "gpu_ids": gpu_ids,
        "n_workers_per_gpu": n_workers_per_gpu,
        "reserve_timeout": reserve_timeout,
        "workdir": mongo_path,
        "mongo_port_address": mongo_port_address,
        "multiple_hosts": multiple_hosts,
    }
    workers_thread = threading.Thread(
        target=launch_workers, kwargs=launcher_kwargs, name="Worker Launcher"
    )
    workers_thread.start()
    started_threads.append(workers_thread)

    # wait for workers and fmin process simultaneously
    workers_done_event = threading.Event()
    fmin_done_event = threading.Event()
    fmin_waiter = threading.Thread(
        target=_wait_for_process_or_thread,
        kwargs={"process": fmin_process, "event": fmin_done_event},
        name="Waiter fmin",
    )
    fmin_waiter.start()
    started_threads.append(fmin_waiter)
    workers_waiter = threading.Thread(
        target=_wait_for_process_or_thread,
        kwargs={"process": workers_thread, "event": workers_done_event},
        name="Waiter workers",
    )
    workers_waiter.start()
    started_threads.append(workers_waiter)
    while not workers_done_event.is_set() and not fmin_done_event.is_set():
        time.sleep(5)

    # when one of them finishes, if it is fmin -> trials should be in the queue
    # if not and not using multiple hosts we wait fmin_timeout seconds for fmin to finish
    # in any case, close waiter threads
    if fmin_done_event.is_set():
        logger.debug("Setting worker watchdog and waiter stop events.")
        stop_watchdog_event.set()
        workers_done_event.set()
    if workers_done_event.is_set() and not multiple_hosts:
        logger.debug("Setting fmin waiter stop event.")
        fmin_done_event.set()
    try:
        if multiple_hosts is not None:
            # if using multiple_hosts, there could still be workers -> disable fmin timeout
            fmin_timeout = None
            logger.debug(
                "multiple_hosts set to True, fmin will block until all trials have been completed."
            )
        else:
            logger.debug(
                "multiple_hosts set to false, Fmin has {time} seconds to finish".format(
                    time=fmin_timeout
                )
            )
        trials = queue.get(timeout=fmin_timeout)
    except Empty:
        logger.error(
            "Queue still empty {fmin_timeout} seconds after all workers have died."
            "\n".format(fmin_timeout=fmin_timeout) + "Terminating minimization process."
        )
        raise FminTimeoutError(
            "Queue still empty {fmin_timeout} seconds after all workers "
            "have died. Check that you have used a new exp_key or allowed "
            "a higher max_evals".format(fmin_timeout=fmin_timeout)
        )

    # sanity: wait for fmin, terminate workers and wait for launcher
    fmin_process.join()
    stop_watchdog_event.set()
    workers_thread.join()
    logger.info(
        "Finished minimization procedure for experiment {exp_key}.".format(
            exp_key=exp_key
        )
    )
    logger.debug("Terminating mongod process.")
    mongod_process.terminate()

    # cleanup processes, threads and files
    _cleanup_processes_files()

    return trials


@_cleanup_decorator
def _fmin_parallel(
    queue: multiprocessing.Queue,
    fn: Callable,
    exp_key: str,
    space: dict,
    algo: Callable = tpe.suggest,
    max_evals: int = 100,
    fmin_timer: float = None,
    show_progressbar: bool = False,
    mongo_port_address: str = "localhost:1234/scvi_db",
):
    """Launches a ``hyperopt`` minimization procedure.
    """
    logger.debug("Instantiating trials object.")
    # instantiate Trials object
    trials = MongoTrials(
        as_mongo_str(os.path.join(mongo_port_address, "jobs")), exp_key=exp_key
    )

    # run hyperoptimization in another fork to enable the use of fmin_timer
    fmin_kwargs = {
        "fn": fn,
        "space": space,
        "algo": algo,
        "max_evals": max_evals,
        "trials": trials,
        "show_progressbar": show_progressbar,
    }
    fmin_thread = threading.Thread(target=fmin, kwargs=fmin_kwargs)
    logger.debug("Calling fmin.")
    # set fmin thread as daemon so it stops when the main process terminates
    fmin_thread.daemon = True
    fmin_thread.start()
    started_threads.append(fmin_thread)
    if fmin_timer is not None:
        logging.debug(
            "Timer set, fmin will run for at most {timer}".format(timer=fmin_timer)
        )
        start_time = time.monotonic()
        run_time = 0
        while run_time < fmin_timer and fmin_thread.is_alive():
            time.sleep(10)
            run_time = time.monotonic() - start_time
    else:
        logging.debug("No timer, waiting for fmin")
        while True:
            if not fmin_thread.is_alive():
                break
            else:
                time.sleep(10)
    logger.debug("fmin returned or timer ran out.")
    # queue.put uses pickle so remove attribute containing thread.lock
    if hasattr(trials, "handle"):
        logger.debug("Deleting Trial handle for pickling.")
        del trials.handle
    logger.debug("Putting Trials in Queue.")
    queue.put(trials)


def _wait_for_process_or_thread(
    process: Union[multiprocessing.Process, threading.Thread], event: threading.Event
):
    """Waits for a process to finish - breaks and sets ``event`` when it does.
    Can be terminated by setting event from outside or by setting the global ``cleanup_event`` of this module.
    """
    logger.debug("Started waiting for {name}.".format(name=process.name))
    while True:
        # set event and break is process is dead
        if not process.is_alive():
            logger.debug("{name} died. Terminating waiter.".format(name=process.name))
            event.set()
            break
        # break if event was set
        if event.is_set():
            logger.debug(
                "Waiting event for {name} set from outside. "
                "Terminating waiter.".format(name=process.name)
            )
            break
        if cleanup_event.is_set():
            logger.debug(
                "Waiting thread for {name} cleaned up.".format(name=process.name)
            )
            event.set()
            break
        time.sleep(5)


@_cleanup_decorator
def launch_workers(
    stop_watchdog_event: threading.Event(),
    exp_key: str,
    n_cpu_workers: int = None,
    gpu_ids: List[int] = None,
    n_workers_per_gpu: int = 1,
    reserve_timeout: float = 30.0,
    workdir: str = ".",
    mongo_port_address: str = "localhost:1234/scvi_db",
    multiple_hosts: bool = False,
):
    """Launches the local workers which are going to run the jobs required by the minimization process.
    Terminates when the worker_watchdog call finishes.
    Specifically, first ``n_gpu_workers`` are launched per GPU in ``gpu_ids`` in their own spawned process.
    Then, ``n_cpu_workers`` CPU workers are launched, also in their own spawned process.
    The use of spawned processes (each have their own python interpreter) is mandatory for compatiblity with CUDA.
    See https://pytorch.org/docs/stable/notes/multiprocessing.html for more information.

    :param stop_watchdog_event: When set, this event stops the watchdog Thread
        which checks that local workers are still running.
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
    :param mongo_port_address: Address to the running MongoDb service.
    :param multiple_hosts: ``True`` if launching workers form multiple hosts.
    """
    # prepare parallel logging
    _logging_queue = spawn_ctx.Queue()
    listener = QueueListener(_logging_queue, DispatchHandler())
    listener.start()
    started_processes.append(listener)

    if gpu_ids is None:
        n_gpus = torch.cuda.device_count()
        logger.debug(
            "gpu_ids is None, defaulting to all {n_gpus} GPUs found by torch.".format(
                n_gpus=n_gpus
            )
        )
        gpu_ids = list(range(n_gpus))
        if n_gpus and n_cpu_workers is None:
            n_cpu_workers = 0
            logging.debug(
                "Some GPU.s found and n_cpu_wokers is None, defaulting to n_cpu_workers = 0"
            )
        if not n_gpus and n_cpu_workers is None:
            n_cpu_workers = os.cpu_count() - 1
            logging.debug(
                "No GPUs found and n_cpu_wokers is None, defaulting to n_cpu_workers = "
                "{n_cpu_workers} (os.cpu_count() - 1)".format(
                    n_cpu_workers=n_cpu_workers
                )
            )
    if (
        gpu_ids is None
        and (n_cpu_workers == 0 or n_cpu_workers is None)
        and not multiple_hosts
    ):
        raise ValueError("No hardware (cpu/gpu) selected/found.")

    # log progress with queue and progress_listener
    progress_queue = spawn_ctx.Queue()
    prog_listener_kwargs = {
        "progress_queue": progress_queue,
        "logging_queue": _logging_queue,
    }
    prog_listener = spawn_ctx.Process(
        target=progress_listener, kwargs=prog_listener_kwargs, name="Progress listener"
    )
    prog_listener.start()
    started_processes.append(prog_listener)

    running_workers = []
    # launch gpu workers
    logger.info(
        "Starting {n_workers_per_gpu} worker.s for each of the {n_gpus} gpu.s set for use/"
        "found.".format(n_workers_per_gpu=n_workers_per_gpu, n_gpus=len(gpu_ids))
    )
    for gpu_id in gpu_ids:
        for sub_id in range(n_workers_per_gpu):
            worker_kwargs = {
                "progress_queue": progress_queue,
                "logging_queue": _logging_queue,
                "exp_key": exp_key,
                "workdir": workdir,
                "gpu": True,
                "hw_id": str(gpu_id),
                "reserve_timeout": reserve_timeout,
                "mongo_port_address": mongo_port_address,
            }
            p = spawn_ctx.Process(
                target=hyperopt_worker,
                kwargs=worker_kwargs,
                name="Worker GPU " + str(gpu_id) + ":" + str(sub_id),
            )
            p.start()
            running_workers.append(p)

    # launch cpu workers
    # TODO: add cpu affinity?
    logger.info(
        "Starting {n_cpu_workers} cpu worker.s".format(n_cpu_workers=n_cpu_workers)
    )
    for cpu_id in range(n_cpu_workers):
        worker_kwargs = {
            "progress_queue": progress_queue,
            "logging_queue": _logging_queue,
            "exp_key": exp_key,
            "workdir": workdir,
            "gpu": False,
            "hw_id": str(cpu_id),
            "reserve_timeout": reserve_timeout,
            "mongo_port_address": mongo_port_address,
        }
        p = spawn_ctx.Process(
            target=hyperopt_worker,
            kwargs=worker_kwargs,
            name="Worker CPU " + str(cpu_id),
        )
        # FIXME won't terminate if parent is killed (SIGKILL)
        p.start()
        running_workers.append(p)
    started_processes.extend(running_workers)

    # wait or return if all workers have died
    workers_watchdog(running_workers=running_workers, stop_event=stop_watchdog_event)
    logger.debug("Worker watchdog finished, terminating workers and closing listener.")
    for worker in running_workers:
        if worker.is_alive():
            worker.terminate()
    listener.stop()
    prog_listener.terminate()


@_cleanup_decorator
def progress_listener(progress_queue, logging_queue):
    """Listens to workers when they finish a job and logs progress.
    Workers put in the progress_queue when they finish a job
    and when they do this function sends a log to the progress logger.
    """
    # write all logs to queue
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    queue_handler = QueueHandler(logging_queue)
    queue_handler.setLevel(logging.DEBUG)
    root_logger.addHandler(queue_handler)
    logger.debug("Listener listening...")

    progress_logger = logging.getLogger("progress_logger")

    i = 0
    while True:
        # get job done signal
        progress_queue.get()
        i += 1
        logger.info("{i} job.s done".format(i=i))
        # update progress bar through ProgressHandler
        progress_logger.info(None)
        if cleanup_event.is_set():
            break


def hyperopt_worker(
    progress_queue: multiprocessing.Queue,
    logging_queue: multiprocessing.Queue,
    exp_key: str,
    workdir: str = ".",
    gpu: bool = True,
    hw_id: str = None,
    poll_interval: float = 1.0,
    reserve_timeout: float = 30.0,
    mongo_port_address: str = "localhost:1234/scvi_db",
):
    """Launches a ``hyperopt`` ``MongoWorker`` which runs jobs until ``ReserveTimeout`` is raised.

    :param progress_queue: Queue in which to put None when a job is done.
    :param logging_queue: Queue to send logs to using a ``QueueHandler``.
    :param exp_key: This key is used by hyperopt as a suffix to the part of the MongoDb
        which corresponds to the current experiment. In particular, it has to be passed to ``MongoWorker``.
    :param workdir:
    :param gpu: If ``True`` means a GPU is to be used.
    :param hw_id: Id of the GPU to use. set via env variable ``CUDA_VISIBLE_DEVICES``.
    :param poll_interval: Time to wait between attempts to reserve a job.
    :param reserve_timeout: Amount of time, in seconds, a worker tries to reserve a job for
        before throwing a ``ReserveTimeout`` Exception.
    :param mongo_port_address: Addres to the running MongoDb service.
    """
    # write all logs to queue
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    queue_handler = QueueHandler(logging_queue)
    queue_handler.setLevel(logging.DEBUG)
    root_logger.addHandler(queue_handler)
    logger.debug("Worker working...")

    os.environ["CUDA_VISIBLE_DEVICES"] = hw_id if gpu else str()

    # FIXME is this stil necessary?
    sys.path.append(".")

    mjobs = MongoJobs.new_from_connection_str(
        os.path.join(as_mongo_str(mongo_port_address), "jobs")
    )
    mworker = MongoWorker(mjobs, float(poll_interval), workdir=workdir, exp_key=exp_key)

    while True:
        # FIXME we don't protect ourselves from memory leaks, bad cleanup, etc.
        try:
            mworker.run_one(reserve_timeout=float(reserve_timeout))
            progress_queue.put(None)
        except ReserveTimeout:
            logger.debug(
                "Caught ReserveTimeout. "
                "Exiting after failing to reserve job for {time} seconds.".format(
                    time=reserve_timeout
                )
            )
            break


def workers_watchdog(
    running_workers: List[multiprocessing.Process], stop_event: threading.Event()
):
    """Checks that workers in running_workers are stil running.
    If none are running anymore, inform user and finish.
    """
    while True:
        one_alive = False
        for worker in running_workers:
            one_alive = one_alive or worker.is_alive()
        # if all workers are dead, inform user
        if not one_alive:
            logger.debug(
                "All workers have died, check stdout/stderr for error tracebacks."
            )
            break
        if stop_event.is_set():
            logger.debug("Stopping Event set, stopping worker watchdog.")
            break
        if cleanup_event.is_set():
            logger.debug("Cleaning up Event set, stopping worker watchdog.")
            stop_event.set()
            break
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
    logger.debug("Instantiating model")
    model = model_class(
        n_input=gene_dataset.nb_genes,
        n_batch=gene_dataset.n_batches * use_batches,
        **model_tunable_kwargs,
    )

    # define trainer
    logger.debug("Instantiating trainer")
    trainer = trainer_class(model, gene_dataset, **trainer_tunable_kwargs)

    # train model
    logger.debug("Starting training")
    trainer.train(**train_func_tunable_kwargs)
    logger.debug("Finished training")
    elapsed_time = time.monotonic() - start_time
    # if training the best model, return model else return criterion
    if is_best_training:
        return trainer
    else:
        # select metric from early stopping kwargs if possible
        metric = None
        early_stopping_kwargs = trainer_specific_kwargs.get(
            "early_stopping_kwargs", None
        )
        if early_stopping_kwargs is not None:
            metric = early_stopping_kwargs.get("early_stopping_metric", None)

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

        # compute true ll
        loss = trainer.test_set.marginal_ll(n_mc_samples=100)

        logger.debug(
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
