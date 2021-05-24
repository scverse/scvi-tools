import logging
from typing import Optional, Union

import numpy as np
import torch
from pyro import poutine
from pytorch_lightning.callbacks import Callback

from scvi.dataloaders import AnnDataLoader, DataSplitter, DeviceBackedDataSplitter
from scvi.model._utils import parse_use_gpu_arg
from scvi.train import PyroTrainingPlan, TrainRunner
from scvi.utils import track

logger = logging.getLogger(__name__)

Number = Union[int, float]


class PyroJitGuideWarmup(Callback):
    def __init__(self, train_dl) -> None:
        super().__init__()
        self.dl = train_dl

    def on_train_start(self, trainer, pl_module):
        """
        Way to warmup Pyro Guide in an automated way.

        Also device agnostic.
        """

        # warmup guide for JIT
        pyro_guide = pl_module.module.guide
        for tensors in self.dl:
            tens = {k: t.to(pl_module.device) for k, t in tensors.items()}
            args, kwargs = pl_module.module._get_fn_args_from_batch(tens)
            pyro_guide(*args, **kwargs)
            break


class PyroSviTrainMixin:
    """
    This mixin class provides methods for:

    - training models using minibatches and using full data (copies data to GPU only once).
    """

    def train(
        self,
        max_epochs: Optional[int] = None,
        use_gpu: Optional[Union[str, int, bool]] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        early_stopping: bool = False,
        plan_kwargs: Optional[dict] = None,
        **trainer_kwargs,
    ):
        """
        Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset. If `None`, defaults to
            `np.min([round((20000 / n_cells) * 400), 400])`
        use_gpu
            Use default GPU if available (if None or True), or index of GPU to use (if int),
            or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training. If `None`, no minibatching occurs and all data is copied to device (e.g., GPU).
        early_stopping
            Perform early stopping. Additional arguments can be passed in `**kwargs`.
            See :class:`~scvi.train.Trainer` for further options.
        plan_kwargs
            Keyword args for :class:`~scvi.train.TrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        if max_epochs is None:
            n_obs = self.adata.n_obs
            max_epochs = np.min([round((20000 / n_obs) * 1000), 1000])

        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else dict()

        if batch_size is None:
            # use data splitter which moves data to GPU once
            data_splitter = DeviceBackedDataSplitter(
                self.adata,
                train_size=train_size,
                validation_size=validation_size,
                batch_size=batch_size,
                use_gpu=use_gpu,
            )
        else:
            data_splitter = DataSplitter(
                self.adata,
                train_size=train_size,
                validation_size=validation_size,
                batch_size=batch_size,
                use_gpu=use_gpu,
            )
        training_plan = PyroTrainingPlan(pyro_module=self.module, **plan_kwargs)

        es = "early_stopping"
        trainer_kwargs[es] = (
            early_stopping if es not in trainer_kwargs.keys() else trainer_kwargs[es]
        )

        data_splitter.setup()
        if "callbacks" not in trainer_kwargs.keys():
            trainer_kwargs["callbacks"] = []
        trainer_kwargs["callbacks"].append(
            PyroJitGuideWarmup(data_splitter.train_dataloader())
        )

        runner = TrainRunner(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            **trainer_kwargs,
        )
        return runner()


class PyroSampleMixin:
    """
    This mixin class provides methods for:

    - generating samples from posterior distribution using minibatches and full data
    """

    @torch.no_grad()
    def _get_one_posterior_sample(
        self,
        args,
        kwargs,
        return_sites: Optional[list] = None,
        sample_observed: bool = False,
    ):
        """Get one sample from posterior distribution.

        Parameters
        ----------
        args
            arguments to model and guide
        kwargs
            arguments to model and guide

        Returns
        -------
        Dictionary with a sample for each variable

        """

        guide_trace = poutine.trace(self.module.guide).get_trace(*args, **kwargs)
        model_trace = poutine.trace(
            poutine.replay(self.module.model, guide_trace)
        ).get_trace(*args, **kwargs)

        sample = {
            name: site["value"].cpu().numpy()
            for name, site in model_trace.nodes.items()
            if (
                (site["type"] == "sample")  # sample statement
                and (
                    (return_sites is None) or (name in return_sites)
                )  # selected in return_sites list
                and (
                    (
                        (not site.get("is_observed", True)) or sample_observed
                    )  # don't save observed unless requested
                    or (site.get("infer", False).get("_deterministic", False))
                )  # unless it is deterministic
                and not isinstance(
                    site.get("fn", None), poutine.subsample_messenger._Subsample
                )  # don't save plates
            )
        }

        return sample

    def _get_posterior_samples(
        self,
        args,
        kwargs,
        num_samples: int = 1000,
        return_sites: Optional[list] = None,
        return_observed: bool = False,
        show_progress: bool = True,
    ):
        """
        Get many samples from posterior distribution.

        Parameters
        ----------
        args
            arguments to model and guide
        kwargs
            arguments to model and guide
        show_progress
            show progress bar

        Returns
        -------
        Dictionary with array of samples for each variable
        dictionary {variable_name: [array with samples in 0 dimension]}

        """

        samples = self._get_one_posterior_sample(
            args, kwargs, return_sites=return_sites, sample_observed=return_observed
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
                args, kwargs, return_sites=return_sites, sample_observed=return_observed
            )

            # add new sample
            samples = {k: samples[k] + [samples_[k]] for k in samples.keys()}

        return {k: np.array(v) for k, v in samples.items()}

    def _get_obs_plate_return_sites(self, sample_kwargs, obs_plate_sites):

        # check whether any variable requested in return_sites are in obs_plate
        if ("return_sites" in sample_kwargs.keys()) and (
            sample_kwargs["return_sites"] is not None
        ):
            return_sites = np.array(sample_kwargs["return_sites"])
            return_sites = return_sites[np.isin(return_sites, obs_plate_sites)]
            if len(return_sites) == 0:
                return [return_sites]
            else:
                return list(return_sites)
        else:
            return obs_plate_sites

    def _get_obs_plate_sites(self, args, kwargs):

        plate_name = self.module.list_obs_plate_vars["name"]

        # find plate dimension
        trace = poutine.trace(self.module.model).get_trace(*args, **kwargs)
        obs_plate = {
            name: site["cond_indep_stack"][0].dim
            for name, site in trace.nodes.items()
            if site["type"] == "sample"
            if any(f.name == plate_name for f in site["cond_indep_stack"])
        }

        return obs_plate

    def _posterior_samples_minibatch(
        self, use_gpu: bool = None, batch_size: int = 128, **sample_kwargs
    ):
        """
        Generate samples of the posterior distribution in minibatches

        Generate samples of the posterior distribution of each parameter, separating local (minibatch) variables
        and global variables, which is necessary when performing minibatch inference.

        Parameters
        ----------
        use_gpu
            Load model on default GPU if available (if None or True),
            or index of GPU to use (if int), or name of GPU (if str), or use CPU (if False).

        Returns
        -------
        dictionary {variable_name: [array with samples in 0 dimension]}

        Notes
        -----
        Note for developers: requires scVI module property (a dictionary, self.module.list_obs_plate_vars)
        which lists observation/minibatch plate name and variables.
        See PyroBaseModuleClass.list_obs_plate_vars for details of the variables it should contain.
        This dictionary can be returned by model class method self.module.model.list_obs_plate_vars()
        to keep all model-specific variables in one place.

        """

        samples = dict()

        gpus, device = parse_use_gpu_arg(use_gpu)

        if batch_size is None:
            batch_size = self.adata.n_obs
        train_dl = AnnDataLoader(self.adata, shuffle=False, batch_size=batch_size)
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
                obs_plate_sites = self._get_obs_plate_sites(args, kwargs)
                if len(obs_plate_sites) == 0:
                    # if no local variables - don't sample
                    break
                obs_plate_dim = list(obs_plate_sites.values())[0]

                sample_kwargs_obs_plate = sample_kwargs.copy()
                sample_kwargs_obs_plate[
                    "return_sites"
                ] = self._get_obs_plate_return_sites(
                    sample_kwargs, list(obs_plate_sites.keys())
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
        tensor_dict = next(iter(train_dl))
        args, kwargs = self.module._get_fn_args_from_batch(tensor_dict)
        args = [a.to(device) for a in args]
        kwargs = {k: v.to(device) for k, v in kwargs.items()}
        self.to_device(device)

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

    def sample_posterior(
        self,
        num_samples: int = 1000,
        return_sites: Optional[list] = None,
        use_gpu: bool = None,
        batch_size: int = 128,
        sample_kwargs=None,
        return_samples: bool = False,
    ):
        """
        Generate samples from posterior distribution for each parameter
        and compute mean, 5%/95% quantiles, standard deviation.

        Parameters
        ----------
        num_samples
            number of posterior samples to generate.
        return_sites
            get samples for pyro model variable, default is all variables, otherwise list variable names).
        use_gpu
            Load model on default GPU if available (if None or True),
            or index of GPU to use (if int), or name of GPU (if str), or use CPU (if False).
        sample_kwargs
            dictionary with arguments to _get_posterior_samples (see below):
            return_observed
                return observed sites/variables?
        return_samples
            return samples in addition to sample mean, 5%/95% quantile and SD?

        Returns
        -------
        Posterior distribution samples, a dictionary with elements as follows,
         containing dictionaries of numpy arrays for each variable:
            post_sample_means - mean of the distribution for each variable;
            post_sample_q05 - 5% quantile;
            post_sample_q95 - 95% quantile;
            post_sample_sds - standard deviation
            posterior_samples - samples for each variable as numpy arrays of shape `(n_samples, ...)` (Optional)

        """

        sample_kwargs = sample_kwargs if isinstance(sample_kwargs, dict) else dict()
        sample_kwargs["num_samples"] = num_samples
        sample_kwargs["return_sites"] = return_sites

        # sample using minibatches (if full data, data is moved to GPU only once anyway)
        samples = self._posterior_samples_minibatch(
            use_gpu=use_gpu, batch_size=batch_size, **sample_kwargs
        )

        param_names = list(samples.keys())
        results = dict()
        if return_samples:
            results["posterior_samples"] = samples

        results["post_sample_means"] = {v: samples[v].mean(axis=0) for v in param_names}
        results["post_sample_q05"] = {
            v: np.quantile(samples[v], 0.05, axis=0) for v in param_names
        }
        results["post_sample_q95"] = {
            v: np.quantile(samples[v], 0.95, axis=0) for v in param_names
        }
        results["post_sample_sds"] = {v: samples[v].std(axis=0) for v in param_names}

        return results
