import logging
from typing import Union

from typing import Optional

import numpy as np
import pandas as pd
import pyro
import torch
from pyro import poutine
from pyro.infer import SVI
from pytorch_lightning.callbacks import Callback
from tqdm.auto import tqdm

from scvi.dataloaders import AnnDataLoader
from scvi.model._utils import parse_use_gpu_arg
from scvi.train import PyroTrainingPlan, Trainer

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

    @property
    def _plan_class(self):
        return PyroTrainingPlan

    def _train_full_data(
        self,
        max_epochs: Optional[int] = None,
        use_gpu: bool = False,
        plan_kwargs: Optional[dict] = None,
        lr: float = 0.01,
        autoencoding_lr: Optional[float] = None,
        clip_norm: float = 200,
        continue_training: bool = True,
    ):
        """
        Private method for training using full data.

        Parameters
        ----------
        max_epochs
            Number of training epochs / iterations
        use_gpu
            Bool, use gpu?
        plan_kwargs
            Training plan arguments such as optim and loss_fn
        continue_training
            When the model is already trained, should calling .train() continue training? (False = restart training)

        Returns
        -------
        ELBO history in self.module.history_

        """

        args, kwargs = self.module.model._get_fn_args_full_data(self.adata)
        gpus, device = parse_use_gpu_arg(use_gpu)

        args = [a.to(device) for a in args]
        kwargs = {k: v.to(device) for k, v in kwargs.items()}
        self.to_device(device)

        if not continue_training or not self.is_trained_:
            # models share param store, make sure it is cleared before training
            pyro.clear_param_store()
        # initialise guide params (warmup)
        self.module.guide(*args, **kwargs)

        svi = SVI(
            self.module.model,
            self.module.guide,
            # select optimiser, optionally choosing different lr for autoencoding guide
            pyro.optim.ClippedAdam(self._optim_param(lr, autoencoding_lr, clip_norm)),
            loss=plan_kwargs["loss_fn"],
        )

        iter_iterator = tqdm(range(max_epochs))
        hist = []
        for it in iter_iterator:

            loss = svi.step(*args, **kwargs)
            iter_iterator.set_description(
                "Epoch " + "{:d}".format(it) + ", -ELBO: " + "{:.4e}".format(loss)
            )
            hist.append(loss)

            if it % 500 == 0:
                torch.cuda.empty_cache()

        if continue_training and self.is_trained_:
            # add ELBO listory
            hist = self.module.history_ + hist
        self.module.history_ = hist
        self.module.is_trained_ = True
        self.history_ = hist
        self.is_trained_ = True

    def _train_minibatch(
        self,
        max_epochs: Optional[int] = None,
        max_steps: Optional[int] = None,
        use_gpu: bool = False,
        plan_kwargs: Optional[dict] = None,
        trainer_kwargs: Optional[dict] = None,
        lr: float = 0.01,
        optim_kwargs: Optional[dict] = None,
        early_stopping: bool = False,
        continue_training: bool = True,
    ):
        """
        Private method for training using minibatches (scVI interface and pytorch lightning).

        Parameters
        ----------
        max_epochs
            Number of training epochs / iterations
        max_steps
            Number of training steps
        use_gpu
            Bool, use gpu?
        plan_kwargs
            Training plan arguments such as optim and loss_fn
        trainer_kwargs
            Arguments for scvi.train.Trainer.
        optim_kwargs
            optimiser creation arguments to  such as autoencoding_lr, clip_norm, module_names
        early_stopping
            Bool, use early stopping? (not tested)
        continue_training
            When the model is already trained, should calling .train() continue training? (False = restart training)

        Returns
        -------
        ELBO history in self.module.history_

        """

        if not continue_training or not self.is_trained_:
            # models share param store, make sure it is cleared before training
            pyro.clear_param_store()

        gpus, device = parse_use_gpu_arg(use_gpu)
        if max_epochs is None:
            n_obs = self.adata.n_obs
            max_epochs = np.min([round((20000 / n_obs) * 400), 400])

        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else dict()
        trainer_kwargs = trainer_kwargs if isinstance(trainer_kwargs, dict) else dict()
        optim_kwargs = optim_kwargs if isinstance(optim_kwargs, dict) else dict()

        batch_size = self.module.model.batch_size
        # select optimiser, optionally choosing different lr for different parameters
        plan_kwargs["optim"] = pyro.optim.ClippedAdam(
            self._optim_param(lr, **optim_kwargs)
        )

        # create data loader for training
        train_dl = AnnDataLoader(self.adata, shuffle=True, batch_size=batch_size)
        plan = PyroTrainingPlan(self.module, **plan_kwargs)
        es = "early_stopping"
        trainer_kwargs[es] = (
            early_stopping if es not in trainer_kwargs.keys() else trainer_kwargs[es]
        )
        trainer = Trainer(
            gpus=gpus,
            max_epochs=max_epochs,
            max_steps=max_steps,
            callbacks=[PyroJitGuideWarmup(train_dl)],
            **trainer_kwargs
        )
        trainer.fit(plan, train_dl)
        self.module.to(device)

        try:
            if continue_training and self.is_trained_:
                # add ELBO listory
                index = range(
                    len(self.module.history_),
                    len(self.module.history_)
                    + len(trainer.logger.history["train_loss_epoch"]),
                )
                trainer.logger.history["train_loss_epoch"].index = index
                self.module.history_ = pd.concat(
                    [self.module.history_, trainer.logger.history["train_loss_epoch"]]
                )
            else:
                self.module.history_ = trainer.logger.history["train_loss_epoch"]
            self.history_ = self.module.history_
        except AttributeError:
            self.history_ = None

        self.module.is_trained_ = True
        self.is_trained_ = True

    def _optim_param(self, lr: float = 0.01, clip_norm: float = 200):
        # create function which fetches different lr for different parameters
        def optim_param(module_name, param_name):
            return {
                "lr": lr,
                # limit the gradient step from becoming too large
                "clip_norm": clip_norm,
            }

        return optim_param

    def train(self, **kwargs):
        """
        Train the model.

        Parameters
        ----------
        max_epochs
            Number of training epochs / iterations
        max_steps
            Number of training steps
        use_gpu
            Bool, use gpu?
        lr
            Learning rate.
        autoencoding_lr
            Optional, a separate learning rate for encoder network.
        clip_norm
            Gradient clipping norm (useful for preventing exploding gradients,
            which can lead to impossible values and NaN loss).
        trainer_kwargs
            Training plan arguments for scvi.train.PyroTrainingPlan (Excluding optim and loss_fn)
        early_stopping
            Bool, use early stopping? (not tested)

        Returns
        -------
        ELBO history in self.module.history_

        """

        plan_kwargs = {"loss_fn": pyro.infer.Trace_ELBO()}

        batch_size = self.module.model.batch_size

        if batch_size is None:
            # train using full data (faster for small datasets)
            self._train_full_data(plan_kwargs=plan_kwargs, **kwargs)
        else:
            # standard training using minibatches
            self._train_minibatch(plan_kwargs=plan_kwargs, **kwargs)


class PyroSampleMixin:
    """
    This mixin class provides methods for:

    - generating samples from posterior distribution using minibatches and full data
    """

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
            name: site["value"].detach().cpu().numpy()
            for name, site in model_trace.nodes.items()
            if (
                (site["type"] == "sample")  # sample statement
                and (
                    (return_sites is None) or (name in return_sites)
                )  # selected in return_sites list
                and (
                    (
                        (not site.get("is_observed", True)) or sample_observed
                    )  # don't save observed
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

        for _ in tqdm(
            range(1, num_samples),
            disable=not show_progress,
            desc="Sampling global variables, sample: ",
        ):
            # generate new sample
            samples_ = self._get_one_posterior_sample(
                args, kwargs, return_sites=return_sites, sample_observed=return_observed
            )

            # add new sample
            samples = {k: samples[k] + [samples_[k]] for k in samples.keys()}

        return {k: np.array(v) for k, v in samples.items()}

    def _posterior_samples_full_data(self, use_gpu: bool = True, **sample_kwargs):
        """
        Generate samples from posterior distribution using all data

        Parameters
        ----------
        sample_kwargs
            arguments to _get_posterior_samples

        Returns
        -------
        dictionary {variable_name: [array with samples in 0 dimension]}

        """

        self.module.eval()
        gpus, device = parse_use_gpu_arg(use_gpu)

        args, kwargs = self.module.model._get_fn_args_full_data(self.adata)
        args = [a.to(device) for a in args]
        kwargs = {k: v.to(device) for k, v in kwargs.items()}
        self.to_device(device)

        samples = self._get_posterior_samples(args, kwargs, **sample_kwargs)

        return samples

    def _posterior_samples_minibatch(self, use_gpu: bool = True, **sample_kwargs):
        """
        Generate samples of the posterior distribution of each parameter, separating local (minibatch) variables
        and global variables, which is necessary when performing minibatch inference.

        Note for developers: requires model class method which lists observation/minibatch plate
        variables (self.module.model.list_obs_plate_vars()).

        Parameters
        ----------
        use_gpu
            Bool, use gpu?

        Returns
        -------
        dictionary {variable_name: [array with samples in 0 dimension]}

        """

        gpus, device = parse_use_gpu_arg(use_gpu)

        self.module.eval()

        train_dl = AnnDataLoader(
            self.adata, shuffle=False, batch_size=self.module.model.batch_size
        )
        # sample local parameters
        i = 0
        with tqdm(train_dl, desc="Sampling local variables, batch: ") as tqdm_dl:
            for tensor_dict in tqdm_dl:
                if i == 0:
                    args, kwargs = self.module._get_fn_args_from_batch(tensor_dict)
                    args = [a.to(device) for a in args]
                    kwargs = {k: v.to(device) for k, v in kwargs.items()}
                    self.to_device(device)

                    # check whether any variable requested in return_sites are in obs_plate
                    sample_kwargs_obs_plate = sample_kwargs.copy()
                    if ("return_sites" in sample_kwargs.keys()) and (
                        sample_kwargs["return_sites"] is not None
                    ):
                        return_sites = np.array(sample_kwargs["return_sites"])
                        return_sites = return_sites[
                            np.isin(
                                return_sites,
                                list(
                                    self.module.model.list_obs_plate_vars()[
                                        "sites"
                                    ].keys()
                                ),
                            )
                        ]
                        if len(return_sites) == 0:
                            sample_kwargs_obs_plate["return_sites"] = [return_sites]
                        else:
                            sample_kwargs_obs_plate["return_sites"] = list(return_sites)
                    else:
                        sample_kwargs_obs_plate["return_sites"] = list(
                            self.module.model.list_obs_plate_vars()["sites"].keys()
                        )
                    sample_kwargs_obs_plate["show_progress"] = False
                    samples = self._get_posterior_samples(
                        args, kwargs, **sample_kwargs_obs_plate
                    )

                    # find plate dimension
                    trace = poutine.trace(self.module.model).get_trace(*args, **kwargs)
                    obs_plate = {
                        name: site["cond_indep_stack"][0].dim
                        for name, site in trace.nodes.items()
                        if site["type"] == "sample"
                        if any(
                            f.name == self.module.model.list_obs_plate_vars()["name"]
                            for f in site["cond_indep_stack"]
                        )
                    }

                else:
                    args, kwargs = self.module._get_fn_args_from_batch(tensor_dict)
                    args = [a.to(device) for a in args]
                    kwargs = {k: v.to(device) for k, v in kwargs.items()}
                    self.to_device(device)

                    samples_ = self._get_posterior_samples(
                        args, kwargs, **sample_kwargs_obs_plate
                    )
                    samples = {
                        k: np.array(
                            [
                                np.concatenate(
                                    [samples[k][i], samples_[k][i]],
                                    axis=list(obs_plate.values())[0],
                                )
                                for i in range(len(samples[k]))
                            ]
                        )
                        for k in samples.keys()
                    }
                i += 1

        # sample global parameters
        i = 0
        for tensor_dict in train_dl:
            if i == 0:
                args, kwargs = self.module._get_fn_args_from_batch(tensor_dict)
                args = [a.to(device) for a in args]
                kwargs = {k: v.to(device) for k, v in kwargs.items()}
                self.to_device(device)

                global_samples = self._get_posterior_samples(
                    args, kwargs, **sample_kwargs
                )
                global_samples = {
                    k: global_samples[k]
                    for k in global_samples.keys()
                    if k not in self.module.model.list_obs_plate_vars()["sites"]
                }
            i += 1

        for k in global_samples.keys():
            samples[k] = global_samples[k]

        self.module.to(device)

        return samples

    def sample_posterior(
        self,
        num_samples: int = 1000,
        return_sites: Optional[list] = None,
        use_gpu: bool = False,
        sample_kwargs=None,
        return_samples: bool = False,
    ):
        """
        Generate samples from posterior distribution for each parameter

        Parameters
        ----------
        num_samples
            number of posterior samples to generate.
        return_sites
            get samples for pyro model variable, default is all variables, otherwise list variable names).
        use_gpu
            Use gpu?
        sample_kwargs
            dictionary with arguments to _get_posterior_samples (see below):
            return_observed
                return observed sites/variables?
        return_samples
            return samples in addition to sample mean, 5%/95% quantile and SD?

        Returns
        -------
        Posterior distribution samples, a dictionary for each of (mean, 5% quantile, SD, optionally all samples),
            containing dictionaries for each variable with numpy arrays.
        Dictionary of all samples contains samples for each variable as numpy arrays of shape ``(n_samples, ...)``

        """

        sample_kwargs = sample_kwargs if isinstance(sample_kwargs, dict) else dict()
        sample_kwargs["num_samples"] = num_samples
        sample_kwargs["return_sites"] = return_sites

        if self.module.model.batch_size is None:
            # sample using full data
            samples = self._posterior_samples_full_data(
                use_gpu=use_gpu, **sample_kwargs
            )
        else:
            # sample using minibatches
            samples = self._posterior_samples_minibatch(
                use_gpu=use_gpu, **sample_kwargs
            )

        param_names = list(samples.keys())
        results = dict()
        if return_samples:
            results["posterior_samples"] = samples

        results["post_sample_means"] = {v: samples[v].mean(axis=0) for v in param_names}
        results["post_sample_q05"] = self.posterior_quantile(q=0.05, use_gpu=use_gpu)
        results["post_sample_q95"] = self.posterior_quantile(q=0.95, use_gpu=use_gpu)
        results["post_sample_sds"] = {v: samples[v].std(axis=0) for v in param_names}

        return results
