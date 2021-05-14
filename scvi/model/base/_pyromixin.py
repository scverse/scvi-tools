import logging
from typing import Optional, Union

import numpy as np
import pandas as pd
import pyro
import torch
from pyro import poutine
from pyro.infer import SVI
from pytorch_lightning.callbacks import Callback

from scvi.dataloaders import AnnDataLoader, DataSplitter
from scvi.model._utils import parse_use_gpu_arg
from scvi.train import PyroTrainingPlan, Trainer, TrainRunner
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
            Minibatch size to use during training.
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

        data_splitter = DataSplitter(
            self.adata,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            use_gpu=use_gpu,
        )
        training_plan = PyroTrainingPlan(
            pyro_module=self.module, n_obs=len(data_splitter.train_idx), **plan_kwargs
        )

        es = "early_stopping"
        trainer_kwargs[es] = (
            early_stopping if es not in trainer_kwargs.keys() else trainer_kwargs[es]
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

    def _train_full_data(
        self,
        max_epochs: Optional[int] = None,
        use_gpu: bool = False,
        plan_kwargs: Optional[dict] = None,
        lr: float = 0.01,
        optim_kwargs: Optional[dict] = None,
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
        optim_kwargs = optim_kwargs if isinstance(optim_kwargs, dict) else dict()

        args = [a.to(device) for a in args]
        kwargs = {k: v.to(device) for k, v in kwargs.items()}
        self.to_device(device)

        if hasattr(self.module.model, "n_obs"):
            setattr(self.module.model, "n_obs", self.adata.n_obs)
        if hasattr(self.module.guide, "n_obs"):
            setattr(self.module.guide, "n_obs", self.adata.n_obs)

        if not continue_training or not self.is_trained_:
            # models share param store, make sure it is cleared before training
            pyro.clear_param_store()
        # initialise guide params (warmup)
        self.module.guide(*args, **kwargs)

        svi = SVI(
            self.module.model,
            self.module.guide,
            # select optimiser, optionally choosing different lr for autoencoding guide
            pyro.optim.ClippedAdam(self._optim_param(lr, **optim_kwargs)),
            loss=plan_kwargs["loss_fn"],
        )

        iter_iterator = track(
            range(max_epochs),
            style="tqdm",
        )
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
        plan_kwargs["n_obs"] = self.adata.n_obs
        trainer_kwargs = trainer_kwargs if isinstance(trainer_kwargs, dict) else dict()
        optim_kwargs = optim_kwargs if isinstance(optim_kwargs, dict) else dict()

        batch_size = self.batch_size
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
            # callbacks=[PyroJitGuideWarmup(train_dl)],
            **trainer_kwargs,
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

    def train_v2(self, **kwargs):
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

        batch_size = self.batch_size

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

    def _check_obs_plate_return_sites(self, sample_kwargs):
        # check whether any variable requested in return_sites are in obs_plate
        obs_plate_sites = list(self.module.model.list_obs_plate_vars()["sites"].keys())
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

    def _find_plate_dimension(self, args, kwargs):

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
        return obs_plate

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

        train_dl = AnnDataLoader(self.adata, shuffle=False, batch_size=self.batch_size)
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
                sample_kwargs_obs_plate = sample_kwargs.copy()
                sample_kwargs_obs_plate[
                    "return_sites"
                ] = self._check_obs_plate_return_sites(sample_kwargs)
                sample_kwargs_obs_plate["show_progress"] = False
                obs_plate = self._find_plate_dimension(args, kwargs)
                obs_plate_dim = list(obs_plate.values())[0]
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
        for tensor_dict in train_dl:
            args, kwargs = self.module._get_fn_args_from_batch(tensor_dict)
            args = [a.to(device) for a in args]
            kwargs = {k: v.to(device) for k, v in kwargs.items()}
            self.to_device(device)

            global_samples = self._get_posterior_samples(args, kwargs, **sample_kwargs)
            global_samples = {
                k: v
                for k, v in global_samples.items()
                if k
                not in list(self.module.model.list_obs_plate_vars()["sites"].keys())
            }
            break

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

        if self.batch_size is None:
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
        results["post_sample_q05"] = {
            v: np.quantile(samples[v], 0.05, axis=0) for v in param_names
        }
        results["post_sample_q95"] = {
            v: np.quantile(samples[v], 0.95, axis=0) for v in param_names
        }
        results["post_sample_sds"] = {v: samples[v].std(axis=0) for v in param_names}

        return results
