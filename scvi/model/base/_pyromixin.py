from __future__ import annotations

import logging
from typing import Callable

import numpy as np
import torch
from lightning.pytorch.callbacks import Callback
from pyro import poutine

from scvi import settings
from scvi.dataloaders import AnnDataLoader, DataSplitter, DeviceBackedDataSplitter
from scvi.model._utils import get_max_epochs_heuristic, parse_device_args
from scvi.train import PyroTrainingPlan, TrainRunner
from scvi.utils import track
from scvi.utils._docstrings import devices_dsp

logger = logging.getLogger(__name__)


class PyroJitGuideWarmup(Callback):
    """A callback to warmup a Pyro guide.

    This helps initialize all the relevant parameters by running
    one minibatch through the Pyro model.
    """

    def __init__(self, dataloader: AnnDataLoader = None) -> None:
        super().__init__()
        self.dataloader = dataloader

    def on_train_start(self, trainer, pl_module):
        """Way to warmup Pyro Guide in an automated way.

        Also device agnostic.
        """
        # warmup guide for JIT
        pyro_guide = pl_module.module.guide
        if self.dataloader is None:
            dl = trainer.datamodule.train_dataloader()
        else:
            dl = self.dataloader
        for tensors in dl:
            tens = {k: t.to(pl_module.device) for k, t in tensors.items()}
            args, kwargs = pl_module.module._get_fn_args_from_batch(tens)
            pyro_guide(*args, **kwargs)
            break


class PyroModelGuideWarmup(Callback):
    """A callback to warmup a Pyro guide and model.

    This helps initialize all the relevant parameters by running
    one minibatch through the Pyro model. This warmup occurs on the CPU.
    """

    def __init__(self, dataloader: AnnDataLoader) -> None:
        super().__init__()
        self.dataloader = dataloader

    def setup(self, trainer, pl_module, stage=None):
        """Way to warmup Pyro Model and Guide in an automated way.

        Setup occurs before any device movement, so params are iniitalized on CPU.
        """
        if stage == "fit":
            pyro_guide = pl_module.module.guide
            dl = self.dataloader
            for tensors in dl:
                tens = {k: t.to(pl_module.device) for k, t in tensors.items()}
                args, kwargs = pl_module.module._get_fn_args_from_batch(tens)
                pyro_guide(*args, **kwargs)
                break


class PyroSviTrainMixin:
    """Mixin class for training Pyro models.

    Training using minibatches and using full data (copies data to GPU only once).
    """

    _data_splitter_cls = DataSplitter
    _training_plan_cls = PyroTrainingPlan
    _train_runner_cls = TrainRunner

    @devices_dsp.dedent
    def train(
        self,
        max_epochs: int | None = None,
        accelerator: str = "auto",
        device: int | str = "auto",
        train_size: float = 0.9,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        batch_size: int = 128,
        early_stopping: bool = False,
        lr: float | None = None,
        training_plan: PyroTrainingPlan | None = None,
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        **trainer_kwargs,
    ):
        """Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset. If `None`, defaults to
            `np.min([round((20000 / n_cells) * 400), 400])`
        %(param_accelerator)s
        %(param_device)s
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        shuffle_set_split
            Whether to shuffle indices before splitting. If `False`, the val, train, and test set are split in the
            sequential order of the data according to `validation_size` and `train_size` percentages.
        batch_size
            Minibatch size to use during training. If `None`, no minibatching occurs and all
            data is copied to device (e.g., GPU).
        early_stopping
            Perform early stopping. Additional arguments can be passed in `**kwargs`.
            See :class:`~scvi.train.Trainer` for further options.
        lr
            Optimiser learning rate (default optimiser is :class:`~pyro.optim.ClippedAdam`).
            Specifying optimiser via plan_kwargs overrides this choice of lr.
        training_plan
            Training plan :class:`~scvi.train.PyroTrainingPlan`.
        datasplitter_kwargs
            Additional keyword arguments passed into :class:`~scvi.dataloaders.DataSplitter`.
        plan_kwargs
            Keyword args for :class:`~scvi.train.PyroTrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        if max_epochs is None:
            max_epochs = get_max_epochs_heuristic(self.adata.n_obs, epochs_cap=1000)

        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else {}
        if lr is not None and "optim" not in plan_kwargs.keys():
            plan_kwargs.update({"optim_kwargs": {"lr": lr}})

        datasplitter_kwargs = datasplitter_kwargs or {}

        if batch_size is None:
            # use data splitter which moves data to GPU once
            data_splitter = DeviceBackedDataSplitter(
                self.adata_manager,
                train_size=train_size,
                validation_size=validation_size,
                batch_size=batch_size,
                accelerator=accelerator,
                device=device,
                **datasplitter_kwargs,
            )
        else:
            data_splitter = self._data_splitter_cls(
                self.adata_manager,
                train_size=train_size,
                validation_size=validation_size,
                shuffle_set_split=shuffle_set_split,
                batch_size=batch_size,
                **datasplitter_kwargs,
            )

        if training_plan is None:
            training_plan = self._training_plan_cls(self.module, **plan_kwargs)

        es = "early_stopping"
        trainer_kwargs[es] = (
            early_stopping if es not in trainer_kwargs.keys() else trainer_kwargs[es]
        )

        if "callbacks" not in trainer_kwargs.keys():
            trainer_kwargs["callbacks"] = []
        trainer_kwargs["callbacks"].append(PyroJitGuideWarmup())

        runner = self._train_runner_cls(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=device,
            **trainer_kwargs,
        )
        return runner()


class PyroSampleMixin:
    """Mixin class for generating samples from posterior distribution.

    Works using both minibatches and full data.
    """

    @torch.inference_mode()
    def _get_one_posterior_sample(
        self,
        args,
        kwargs,
        return_sites: list | None = None,
        return_observed: bool = False,
    ):
        """Get one sample from posterior distribution.

        Parameters
        ----------
        args
            arguments to model and guide
        kwargs
            arguments to model and guide
        return_sites
            List of variables for which to generate posterior samples, defaults to all variables.
        return_observed
            Record samples of observed variables.

        Returns
        -------
        Dictionary with a sample for each variable
        """
        if isinstance(self.module.guide, poutine.messenger.Messenger):
            # This already includes trace-replay behavior.
            sample = self.module.guide(*args, **kwargs)
        else:
            guide_trace = poutine.trace(self.module.guide).get_trace(*args, **kwargs)
            model_trace = poutine.trace(
                poutine.replay(self.module.model, guide_trace)
            ).get_trace(*args, **kwargs)
            sample = {
                name: site["value"]
                for name, site in model_trace.nodes.items()
                if (
                    (site["type"] == "sample")  # sample statement
                    and (
                        (return_sites is None) or (name in return_sites)
                    )  # selected in return_sites list
                    and (
                        (
                            (not site.get("is_observed", True)) or return_observed
                        )  # don't save observed unless requested
                        or (site.get("infer", False).get("_deterministic", False))
                    )  # unless it is deterministic
                    and not isinstance(
                        site.get("fn", None), poutine.subsample_messenger._Subsample
                    )  # don't save plates
                )
            }

        sample = {name: site.cpu().numpy() for name, site in sample.items()}

        return sample

    def _get_posterior_samples(
        self,
        args,
        kwargs,
        num_samples: int = 1000,
        return_sites: list | None = None,
        return_observed: bool = False,
        show_progress: bool = True,
    ):
        """Get many (num_samples=N) samples from posterior distribution.

        Parameters
        ----------
        args
            arguments to model and guide
        kwargs
            keyword arguments to model and guide
        return_sites
            List of variables for which to generate posterior samples, defaults to all variables.
        return_observed
            Record samples of observed variables.
        show_progress
            show progress bar

        Returns
        -------
        Dictionary with array of samples for each variable
        dictionary {variable_name: [array with samples in 0 dimension]}
        """
        samples = self._get_one_posterior_sample(
            args, kwargs, return_sites=return_sites, return_observed=return_observed
        )
        samples = {k: [v] for k, v in samples.items()}

        for _ in track(
            range(1, num_samples),
            style="tqdm",
            description="Sampling global variables, sample: ",
            disable=not show_progress,
        ):
            # generate new sample
            samples_ = self._get_one_posterior_sample(
                args, kwargs, return_sites=return_sites, return_observed=return_observed
            )

            # add new sample
            samples = {k: samples[k] + [samples_[k]] for k in samples.keys()}

        return {k: np.array(v) for k, v in samples.items()}

    def _get_obs_plate_return_sites(self, return_sites, obs_plate_sites):
        """Check return_sites for overlap with observation/minibatch plate sites."""
        # check whether any variable requested in return_sites are in obs_plate
        if return_sites is not None:
            return_sites = np.array(return_sites)
            return_sites = return_sites[np.isin(return_sites, obs_plate_sites)]
            if len(return_sites) == 0:
                return [return_sites]
            else:
                return list(return_sites)
        else:
            return obs_plate_sites

    def _get_obs_plate_sites(
        self,
        args: list,
        kwargs: dict,
        return_observed: bool = False,
    ):
        """Automatically guess which model sites belong to observation/minibatch plate.

        This function requires minibatch plate name specified in `self.module.list_obs_plate_vars["name"]`.

        Parameters
        ----------
        args
            Arguments to the model.
        kwargs
            Keyword arguments to the model.
        return_observed
            Record samples of observed variables.

        Returns
        -------
        Dictionary with keys corresponding to site names and values to plate dimension.
        """
        plate_name = self.module.list_obs_plate_vars["name"]

        # find plate dimension
        trace = poutine.trace(self.module.model).get_trace(*args, **kwargs)
        obs_plate = {
            name: site["cond_indep_stack"][0].dim
            for name, site in trace.nodes.items()
            if (
                (site["type"] == "sample")  # sample statement
                and (
                    (
                        (not site.get("is_observed", True)) or return_observed
                    )  # don't save observed unless requested
                    or (site.get("infer", False).get("_deterministic", False))
                )  # unless it is deterministic
                and not isinstance(
                    site.get("fn", None), poutine.subsample_messenger._Subsample
                )  # don't save plates
            )
            if any(f.name == plate_name for f in site["cond_indep_stack"])
        }

        return obs_plate

    @devices_dsp.dedent
    def _posterior_samples_minibatch(
        self,
        accelerator: str = "auto",
        device: int | str = "auto",
        batch_size: int | None = None,
        **sample_kwargs,
    ):
        """Generate samples of the posterior distribution in minibatches.

        Generate samples of the posterior distribution of each parameter, separating local (minibatch) variables
        and global variables, which is necessary when performing minibatch inference.

        Parameters
        ----------
        %(param_accelerator)s
        %(param_device)s
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        dictionary {variable_name: [array with samples in 0 dimension]}
        """
        samples = {}

        _, _, device = parse_device_args(
            accelerator=accelerator,
            devices=device,
            return_device="torch",
            validate_single_device=True,
        )

        batch_size = batch_size if batch_size is not None else settings.batch_size

        train_dl = AnnDataLoader(
            self.adata_manager, shuffle=False, batch_size=batch_size
        )
        # sample local parameters
        i = 0
        for tensor_dict in track(
            train_dl,
            style="tqdm",
            description="Sampling local variables, batch: ",
        ):
            args, kwargs = self.module._get_fn_args_from_batch(tensor_dict)
            args = [a.to(device) for a in args]
            kwargs = {k: v.to(device) for k, v in kwargs.items()}
            self.to_device(device)

            if i == 0:
                return_observed = getattr(sample_kwargs, "return_observed", False)
                obs_plate_sites = self._get_obs_plate_sites(
                    args, kwargs, return_observed=return_observed
                )
                if len(obs_plate_sites) == 0:
                    # if no local variables - don't sample
                    break
                obs_plate_dim = list(obs_plate_sites.values())[0]

                sample_kwargs_obs_plate = sample_kwargs.copy()
                sample_kwargs_obs_plate[
                    "return_sites"
                ] = self._get_obs_plate_return_sites(
                    sample_kwargs["return_sites"], list(obs_plate_sites.keys())
                )
                sample_kwargs_obs_plate["show_progress"] = False

                samples = self._get_posterior_samples(
                    args, kwargs, **sample_kwargs_obs_plate
                )
            else:
                samples_ = self._get_posterior_samples(
                    args, kwargs, **sample_kwargs_obs_plate
                )

                samples = {
                    k: np.array(
                        [
                            np.concatenate(
                                [samples[k][j], samples_[k][j]],
                                axis=obs_plate_dim,
                            )
                            for j in range(
                                len(samples[k])
                            )  # for each sample (in 0 dimension
                        ]
                    )
                    for k in samples.keys()  # for each variable
                }
            i += 1

        # sample global parameters
        global_samples = self._get_posterior_samples(args, kwargs, **sample_kwargs)
        global_samples = {
            k: v
            for k, v in global_samples.items()
            if k not in list(obs_plate_sites.keys())
        }

        for k in global_samples.keys():
            samples[k] = global_samples[k]

        self.module.to(device)

        return samples

    @devices_dsp.dedent
    def sample_posterior(
        self,
        num_samples: int = 1000,
        return_sites: list | None = None,
        accelerator: str = "auto",
        device: int | str = "auto",
        batch_size: int | None = None,
        return_observed: bool = False,
        return_samples: bool = False,
        summary_fun: dict[str, Callable] | None = None,
    ):
        """Summarise posterior distribution.

        Generate samples from posterior distribution for each parameter
        and compute mean, 5th/95th quantiles, standard deviation.

        Parameters
        ----------
        num_samples
            Number of posterior samples to generate.
        return_sites
            List of variables for which to generate posterior samples, defaults to all variables.
        %(param_accelerator)s
        %(param_device)s
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_observed
            Return observed sites/variables? Observed count matrix can be very large so not returned by default.
        return_samples
            Return all generated posterior samples in addition to sample mean, 5th/95th quantile and SD?
        summary_fun
             a dict in the form {"means": np.mean, "std": np.std} which specifies posterior distribution
             summaries to compute and which names to use. See below for default returns.

        Returns
        -------
        post_sample_means: Dict[str, :class:`np.ndarray`]
            Mean of the posterior distribution for each variable, a dictionary of numpy arrays for each variable;
        post_sample_q05: Dict[str, :class:`np.ndarray`]
            5th quantile of the posterior distribution for each variable;
        post_sample_q05: Dict[str, :class:`np.ndarray`]
            95th quantile of the posterior distribution for each variable;
        post_sample_q05: Dict[str, :class:`np.ndarray`]
            Standard deviation of the posterior distribution for each variable;
        posterior_samples: Optional[Dict[str, :class:`np.ndarray`]]
            Posterior distribution samples for each variable as numpy arrays of shape `(n_samples, ...)` (Optional).

        Notes
        -----
        Note for developers: requires overwritten :attr:`~scvi.module.base.PyroBaseModuleClass.list_obs_plate_vars` property.
        which lists observation/minibatch plate name and variables.
        See :attr:`~scvi.module.base.PyroBaseModuleClass.list_obs_plate_vars` for details of the variables it should contain.
        This dictionary can be returned by model class property `self.module.model.list_obs_plate_vars`
        to keep all model-specific variables in one place.
        """
        # sample using minibatches (if full data, data is moved to GPU only once anyway)
        samples = self._posterior_samples_minibatch(
            accelerator=accelerator,
            device=device,
            batch_size=batch_size,
            num_samples=num_samples,
            return_sites=return_sites,
            return_observed=return_observed,
        )

        param_names = list(samples.keys())
        results = {}
        if return_samples:
            results["posterior_samples"] = samples

        if summary_fun is None:
            summary_fun = {
                "means": np.mean,
                "stds": np.std,
                "q05": lambda x, axis: np.quantile(x, 0.05, axis=axis),
                "q95": lambda x, axis: np.quantile(x, 0.95, axis=axis),
            }
        for k, fun in summary_fun.items():
            results[f"post_sample_{k}"] = {
                v: fun(samples[v], axis=0) for v in param_names
            }

        return results
