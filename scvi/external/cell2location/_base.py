from datetime import date
from typing import Optional

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyro
import torch
from pyro import poutine
from pyro.infer import SVI
from pyro.infer.autoguide import AutoNormal, init_to_mean
from tqdm.auto import tqdm

from scvi.dataloaders import AnnDataLoader
from scvi.model._utils import parse_use_gpu_arg
from scvi.module.base import PyroBaseModuleClass
from scvi.train import PyroTrainingPlan, Trainer

from .autoguide import AutoGuideList, AutoNormalEncoder


class Cell2locationBaseModule(PyroBaseModuleClass):
    def __init__(
        self,
        model,
        amortised: bool = False,
        single_encoder: bool = True,
        encoder_kwargs=None,
        data_transform="log1p",
        **kwargs
    ):
        """
        Module class which defines AutoGuide given model. Supports multiple model architectures.

        Parameters
        ----------
        amortised
            boolean, use a Neural Network to approximate posterior distribution of location-specific (local) parameters?
        encoder_kwargs
            arguments for Neural Network construction (scvi.nn.FCLayers)
        kwargs
            arguments for specific model class - e.g. number of genes, values of the prior distribution
        """
        super().__init__()
        self.hist = []

        self._model = model(**kwargs)
        self._amortised = amortised

        if not amortised:
            self._guide = AutoNormal(
                self.model,
                init_loc_fn=init_to_mean,
                create_plates=self.model.create_plates,
            )
        else:
            encoder_kwargs = (
                encoder_kwargs if isinstance(encoder_kwargs, dict) else dict()
            )
            n_hidden = (
                encoder_kwargs["n_hidden"]
                if "n_hidden" in encoder_kwargs.keys()
                else 200
            )
            amortised_vars = self.model.list_obs_plate_vars()
            self._guide = AutoGuideList(
                self.model, create_plates=self.model.create_plates
            )
            self._guide.append(
                AutoNormal(
                    pyro.poutine.block(
                        self.model, hide=list(amortised_vars["sites"].keys())
                    ),
                    init_loc_fn=init_to_mean,
                )
            )
            if isinstance(data_transform, np.ndarray):
                # add extra info about gene clusters to the network
                self.register_buffer(
                    "gene_clusters", torch.tensor(data_transform.astype("float32"))
                )
                n_in = self.model.n_vars + data_transform.shape[1]
                data_transform = self.data_transform_clusters()
            elif data_transform == "log1p":
                # use simple log1p transform
                data_transform = torch.log1p
                n_in = self.model.n_vars
            elif (
                isinstance(data_transform, dict)
                and "var_std" in list(data_transform.keys())
                and "var_mean" in list(data_transform.keys())
            ):
                # use data transform by scaling
                n_in = self.model.n_vars
                self.register_buffer(
                    "var_mean",
                    torch.tensor(
                        data_transform["var_mean"].astype("float32").reshape((1, n_in))
                    ),
                )
                self.register_buffer(
                    "var_std",
                    torch.tensor(
                        data_transform["var_std"].astype("float32").reshape((1, n_in))
                    ),
                )
                data_transform = self.data_transform_scale()
            else:
                # use custom data transform
                data_transform = data_transform
                n_in = self.model.n_vars

            self._guide.append(
                AutoNormalEncoder(
                    pyro.poutine.block(
                        self.model, expose=list(amortised_vars["sites"].keys())
                    ),
                    amortised_plate_sites=amortised_vars,
                    n_in=n_in,
                    n_hidden=n_hidden,
                    data_transform=data_transform,
                    encoder_kwargs=encoder_kwargs,
                    single_encoder=single_encoder,
                )
            )

        self._get_fn_args_from_batch = self._model._get_fn_args_from_batch

    @property
    def model(self):
        return self._model

    @property
    def guide(self):
        return self._guide

    @property
    def is_amortised(self):
        return self._amortised

    def data_transform_clusters(self):
        def _data_transform(x):
            return torch.log1p(torch.cat([x, x @ self.gene_clusters], dim=1))

        return _data_transform

    def data_transform_scale(self):
        def _data_transform(x):
            # return (x - self.var_mean) / self.var_std
            return x / self.var_std

        return _data_transform


class TrainSampleMixin:
    """
    Reimplementation of cell2location [Kleshchevnikov20]_ model. This mixin class provides methods for:

    - training models using minibatches (standard scVI interface) and using full data (copies data to GPU only once).
    - computing median and quantiles of the posterior distribution using both direct and amortised inference
    - generating samples from posterior distribution

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
        autoencoding_lr: Optional[float] = None,
        clip_norm: float = 200,
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

        batch_size = self.module.model.batch_size
        # initialise guide params (warmup)
        train_dl = AnnDataLoader(self.adata, shuffle=True, batch_size=batch_size)
        for tensor_dict in train_dl:
            args, kwargs = self.module._get_fn_args_from_batch(tensor_dict)
            args = [a.to(device) for a in args]
            kwargs = {k: v.to(device) for k, v in kwargs.items()}
            self.to_device(device)
            self.module.guide(*args, **kwargs)
            break
        # select optimiser, optionally choosing different lr for autoencoding guide (requires initialised guide)
        plan_kwargs["optim"] = pyro.optim.ClippedAdam(
            self._optim_param(lr, autoencoding_lr, clip_norm)
        )

        # create data loader for training
        train_dl = AnnDataLoader(self.adata, shuffle=True, batch_size=batch_size)
        plan = PyroTrainingPlan(self.module, **plan_kwargs)
        es = "early_stopping"
        trainer_kwargs[es] = (
            early_stopping if es not in trainer_kwargs.keys() else trainer_kwargs[es]
        )
        trainer = Trainer(
            gpus=gpus, max_epochs=max_epochs, max_steps=max_steps, **trainer_kwargs
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

    def _optim_param(
        self,
        lr,
        autoencoding_lr,
        clip_norm,
        module_names=["encoder", "hidden2locs", "hidden2scales"],
    ):
        # create function which fetches different lr for autoencoding guide
        def optim_param(module_name, param_name):
            # detect variables in autoencoding guide
            if autoencoding_lr is not None and np.any(
                [n in module_name + "." + param_name for n in module_names]
            ):
                return {
                    "lr": autoencoding_lr,
                    # limit the gradient step from becoming too large
                    "clip_norm": clip_norm,
                }
            else:
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

    def _posterior_quantile_amortised(self, q: float = 0.5, use_gpu: bool = True):
        """
        Compute median of the posterior distribution of each parameter, separating local (minibatch) variable
        and global variables, which is necessary when performing amortised inference.

        Note for developers: requires model class method which lists observation/minibatch plate
        variables (self.module.model.list_obs_plate_vars()).

        Parameters
        ----------
        q
            quantile to compute
        use_gpu
            Bool, use gpu?

        Returns
        -------
        dictionary {variable_name: posterior median}

        """

        gpus, device = parse_use_gpu_arg(use_gpu)

        self.module.eval()

        train_dl = AnnDataLoader(
            self.adata, shuffle=False, batch_size=self.module.model.batch_size
        )

        # sample local parameters
        i = 0
        for tensor_dict in train_dl:
            if i == 0:
                args, kwargs = self.module._get_fn_args_from_batch(tensor_dict)
                args = [a.to(device) for a in args]
                kwargs = {k: v.to(device) for k, v in kwargs.items()}
                self.to_device(device)

                means = self.module.guide.quantiles([q], *args, **kwargs)
                means = {
                    k: means[k].cpu().detach().numpy()
                    for k in means.keys()
                    if k in self.module.model.list_obs_plate_vars()["sites"]
                }

                # find plate dimension
                trace = poutine.trace(self.module.model).get_trace(*args, **kwargs)
                # print(trace.nodes[self.module.model.list_obs_plate_vars()['name']])
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

                means_ = self.module.guide.quantiles([q], *args, **kwargs)
                means_ = {
                    k: means_[k].cpu().detach().numpy()
                    for k in means_.keys()
                    if k in self.module.model.list_obs_plate_vars()["sites"]
                }
                means = {
                    k: np.concatenate(
                        [means[k], means_[k]], axis=list(obs_plate.values())[0]
                    )
                    for k in means.keys()
                }
            i += 1

        # sample global parameters
        for tensor_dict in train_dl:
            args, kwargs = self.module._get_fn_args_from_batch(tensor_dict)
            args = [a.to(device) for a in args]
            kwargs = {k: v.to(device) for k, v in kwargs.items()}
            self.to_device(device)

            global_means = self.module.guide.quantiles([q], *args, **kwargs)
            global_means = {
                k: global_means[k].cpu().detach().numpy()
                for k in global_means.keys()
                if k not in self.module.model.list_obs_plate_vars()["sites"]
            }
            break

        for k in global_means.keys():
            means[k] = global_means[k]

        self.module.to(device)

        return means

    def _posterior_quantile(self, q: float = 0.5, use_gpu: bool = True):
        """
        Compute median of the posterior distribution of each parameter pyro models trained without amortised inference.

        Parameters
        ----------
        q
            quantile to compute
        use_gpu
            Bool, use gpu?

        Returns
        -------
        dictionary {variable_name: posterior median}

        """

        self.module.eval()
        gpus, device = parse_use_gpu_arg(use_gpu)

        if self.module.model.batch_size is None:
            args, kwargs = self.module.model._get_fn_args_full_data(self.adata)
            args = [a.to(device) for a in args]
            kwargs = {k: v.to(device) for k, v in kwargs.items()}
            self.to_device(device)
        else:
            train_dl = AnnDataLoader(
                self.adata, shuffle=False, batch_size=self.module.model.batch_size
            )
            # sample global parameters
            for tensor_dict in train_dl:
                args, kwargs = self.module._get_fn_args_from_batch(tensor_dict)
                args = [a.to(device) for a in args]
                kwargs = {k: v.to(device) for k, v in kwargs.items()}
                self.to_device(device)
                break

        means = self.module.guide.quantiles([q], *args, **kwargs)
        means = {k: means[k].cpu().detach().numpy() for k in means.keys()}

        return means

    def posterior_quantile(self, q: float = 0.5, use_gpu: bool = True):
        """
        Compute median of the posterior distribution of each parameter.

        Parameters
        ----------
        q
            quantile to compute
        use_gpu

        Returns
        -------

        """

        if self.module.is_amortised:
            return self._posterior_quantile_amortised(q=q, use_gpu=use_gpu)
        else:
            return self._posterior_quantile(q=q, use_gpu=use_gpu)

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


class PltExportMixin:
    @staticmethod
    def plot_posterior_mu_vs_data(mu, data):
        r"""Plot expected value of the model (e.g. mean of NB distribution) vs observed data

        :param mu: expected value
        :param data: data value
        """

        plt.hist2d(
            np.log10(data.flatten() + 1),
            np.log10(mu.flatten() + 1),
            bins=50,
            norm=matplotlib.colors.LogNorm(),
        )
        plt.gca().set_aspect("equal", adjustable="box")
        plt.xlabel("Data, log10")
        plt.ylabel("Posterior expected value, log10")
        plt.tight_layout()

    def plot_history(self, iter_start=0, iter_end=-1, ax=None):
        r"""Plot training history

        Parameters
        ----------
        iter_start
            omit initial iterations from the plot
        iter_end
            omit last iterations from the plot
        ax
            matplotlib axis

        """
        if ax is None:
            ax = plt
            ax.set_xlabel = plt.xlabel
            ax.set_ylabel = plt.ylabel
        if iter_end == -1:
            iter_end = len(self.history_)

        ax.plot(
            self.history_.index[iter_start:iter_end],
            np.array(self.history_.values.flatten())[iter_start:iter_end],
            label="train",
        )
        ax.legend()
        ax.xlim(0, len(self.history_))
        ax.set_xlabel("Training epochs")
        ax.set_ylabel("-ELBO loss")
        plt.tight_layout()

    def export2adata(self, samples):
        r"""
        Export key model variable and samples

        Parameters
        ----------
        samples
            dictionary with posterior mean, 5%/95% quantiles, SD, samples, generated by `.sample_posterior()`

        Returns
        -------
        updated dictionary with additional details which can be saved to `adata.uns['mod']`.
        """
        # add factor filter and samples of all parameters to unstructured data
        results = {
            "model_name": str(self.module.__class__.__name__),
            "date": str(date.today()),
            "factor_filter": list(getattr(self, "factor_filter", [])),
            "factor_names": list(self.factor_names_),
            "var_names": self.adata.var_names.tolist(),
            "obs_names": self.adata.obs_names.tolist(),
            "post_sample_means": samples["post_sample_means"],
            "post_sample_sds": samples["post_sample_sds"],
            "post_sample_q05": samples["post_sample_q05"],
            "post_sample_q95": samples["post_sample_q95"],
        }

        return results

    def sample2df_obs(
        self,
        samples: dict,
        site_name: str = "w_sf",
        name_prefix: str = "cell_abundance_",
    ):
        """Export posterior distribution summary for observation-specific parameters
        (e.g. spatial cell abundance) as Pandas data frames
        (means, 5%/95% quantiles and sd of posterior distribution).

        Parameters
        ----------
        samples
            dictionary with posterior mean, 5%/95% quantiles, SD, samples, generated by `.sample_posterior()`
        site_name
            name of the model parameter to be exported

        Returns
        -------
        list with 4 Pandas data frames corresponding to means, 5%/95% quantiles and sd of posterior distribution

        """

        results = dict()

        results["mean_" + name_prefix + site_name] = pd.DataFrame.from_records(
            samples["post_sample_means"][site_name],
            index=self.adata.obs_names,
            columns=["mean_" + name_prefix + site_name + i for i in self.factor_names_],
        )

        results["sd_" + name_prefix + site_name] = pd.DataFrame.from_records(
            samples["post_sample_sds"][site_name],
            index=self.adata.obs_names,
            columns=["sd_" + name_prefix + site_name + i for i in self.factor_names_],
        )

        results["q05_" + name_prefix + site_name] = pd.DataFrame.from_records(
            samples["post_sample_q05"][site_name],
            index=self.adata.obs_names,
            columns=["q05_" + name_prefix + site_name + i for i in self.factor_names_],
        )

        results["q95_" + name_prefix + site_name] = pd.DataFrame.from_records(
            samples["post_sample_q95"][site_name],
            index=self.adata.obs_names,
            columns=["q95_" + name_prefix + site_name + i for i in self.factor_names_],
        )

        return results

    def sample2df_vars(
        self, samples: dict, site_name: str = "gene_factors", name_prefix: str = ""
    ):
        """Export posterior distribution summary for variable-specific parameters as Pandas data frames
        (means, 5%/95% quantiles and sd of posterior distribution).

        Parameters
        ----------
        samples
            dictionary with posterior mean, 5%/95% quantiles, SD, samples, generated by `.sample_posterior()`
        site_name
            name of the model parameter to be exported

        Returns
        -------
        list with 4 Pandas data frames corresponding to means, 5%/95% quantiles and sd of posterior distribution

        """

        results = dict()

        results["mean_" + name_prefix + site_name] = pd.DataFrame.from_records(
            samples["post_sample_means"][site_name].T,
            index=self.adata.var_names,
            columns=["mean_" + name_prefix + site_name + i for i in self.factor_names_],
        )

        results["sd_" + name_prefix + site_name] = pd.DataFrame.from_records(
            samples["post_sample_sds"][site_name].T,
            index=self.adata.var_names,
            columns=["sd_" + name_prefix + site_name + i for i in self.factor_names_],
        )

        results["q05_" + name_prefix + site_name] = pd.DataFrame.from_records(
            samples["post_sample_q05"][site_name].T,
            index=self.adata.var_names,
            columns=["q05_" + name_prefix + site_name + i for i in self.factor_names_],
        )

        results["q95_" + name_prefix + site_name] = pd.DataFrame.from_records(
            samples["post_sample_q95"][site_name].T,
            index=self.adata.var_names,
            columns=["q95_" + name_prefix + site_name + i for i in self.factor_names_],
        )

        return results
