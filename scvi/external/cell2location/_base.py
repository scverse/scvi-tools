from typing import Optional

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyro
import torch
from pyro import poutine
from pyro.infer import SVI
from tqdm.auto import tqdm

from scvi.dataloaders import AnnDataLoader
from scvi.model._utils import parse_use_gpu_arg
from scvi.train import PyroTrainingPlan, Trainer


class Cell2locationTrainSampleMixin:
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
    ):

        samples = self._get_one_posterior_sample(
            args, kwargs, return_sites=return_sites, sample_observed=return_observed
        )
        samples = {k: [v] for k, v in samples.items()}

        for _ in tqdm(range(1, num_samples)):
            # generate new sample
            samples_ = self._get_one_posterior_sample(
                args, kwargs, return_sites=return_sites, sample_observed=return_observed
            )

            # add new sample
            samples = {k: samples[k] + [samples_[k]] for k in samples.keys()}

        return {k: np.array(v) for k, v in samples.items()}

    def _posterior_samples_full_data(self, use_gpu: bool = True, **sample_kwargs):

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
        sample_kwargs
            arguments to _get_posterior_samples (see below)
        num_samples
            number of samples to generate
        return_sites
            list of site/variable names to be sampled (all by default)
        return_observed
            return observed sites/variables?

        Returns
        -------
        dictionary {variable_name: array with samples in 0 dimension}

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

                sample_kwargs_obs_plate = sample_kwargs.copy()
                if "return_sites" in sample_kwargs.keys():
                    # check whether any
                    return_sites = np.array(sample_kwargs["return_sites"])
                    return_sites = return_sites[
                        np.isin(
                            return_sites,
                            list(
                                self.module.model.list_obs_plate_vars()["sites"].keys()
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
        site="all",
        num_samples: int = 1000,
        use_gpu: bool = False,
        sample_kwargs=None,
        save_samples: bool = False,
    ):
        r"""Sample posterior distribution of parameters - either all or single parameter
        :param node: pyro parameter to sample (e.g. default "all", self.spot_factors)
        :param n_samples: number of posterior samples to generate (1000 is recommended, reduce if you get GPU memory error)
        :param save_samples: save samples in addition to sample mean, 5% quantile, SD.
        :param return_samples: return summarised samples in addition to saving them in `self.samples`
        :param mean_field_slot: string, which mean_field slot to sample? 'init_1' by default
        :return: nothing, a dictionary of dictionaries (mean, 5% quantile, SD, optionally all samples) with numpy arrays for each variables is added to self.adata.uns['mod'].
        Optional dictionary of all samples contains parameters as numpy arrays of shape ``(n_samples, ...)``
        """

        self.num_samples = num_samples
        sample_kwargs = sample_kwargs if isinstance(sample_kwargs, dict) else dict()
        sample_kwargs["num_samples"] = num_samples

        if site == "all":
            # Sample all parameters
            if self.module.model.batch_size is None:
                samples = self._posterior_samples_full_data(
                    use_gpu=use_gpu, **sample_kwargs
                )
                return samples
            else:
                samples = self._posterior_samples_minibatch(
                    use_gpu=use_gpu, **sample_kwargs
                )
                return samples

            self.param_names = list(self.adata.uns["mod"]["post_samples"].keys())

            self.adata.uns["mod"]["post_sample_means"] = {
                v: self.adata.uns["mod"]["post_samples"][v].mean(axis=0)
                for v in self.param_names
            }
            self.adata.uns["mod"]["post_sample_q05"] = {
                v: np.quantile(self.adata.uns["mod"]["post_samples"][v], 0.05, axis=0)
                for v in self.param_names
            }
            self.adata.uns["mod"]["post_sample_q95"] = {
                v: np.quantile(self.adata.uns["mod"]["post_samples"][v], 0.95, axis=0)
                for v in self.param_names
            }
            self.adata.uns["mod"]["post_sample_sds"] = {
                v: self.adata.uns["mod"]["post_samples"][v].std(axis=0)
                for v in self.param_names
            }

            if not save_samples:
                del self.adata.uns["mod"]["post_samples"]

        else:
            self.sample_node(
                site, self.n_sampl_batches, batch_size=self.num_samples_batch, suff=""
            )


class PltExportMixin:
    def plot_posterior_mu_vs_data(self, mu_node_name="mu", data_node="X_data"):
        r"""Plot expected value of the model (e.g. mean of poisson distribution)

        :param mu_node_name: name of the object slot containing expected value
        :param data_node: name of the object slot containing data
        """

        if type(mu_node_name) is str:
            mu = getattr(self, mu_node_name)
        else:
            mu = mu_node_name

        if type(data_node) is str:
            data_node = getattr(self, data_node)

        plt.hist2d(
            np.log10(data_node.flatten() + 1),
            np.log10(mu.flatten() + 1),
            bins=50,
            norm=matplotlib.colors.LogNorm(),
        )
        plt.gca().set_aspect("equal", adjustable="box")
        plt.xlabel("Data, log10(nUMI)")
        plt.ylabel("Posterior sample, log10(nUMI)")
        plt.title("UMI counts (all cell, all genes)")
        plt.tight_layout()

    def plot_history(
        self, iter_start=0, iter_end=-1, history_key=None, log_y=True, ax=None
    ):
        r"""Plot training history

        :param iter_start: omit initial iterations from the plot
        :param iter_end: omit last iterations from the plot
        """

        if ax is None:
            ax = plt
            ax.set_xlabel = plt.xlabel
            ax.set_ylabel = plt.ylabel

        if history_key is None:
            history_key = self.history.keys()

        if type(history_key) == str:
            history_key = [history_key]

        for i in history_key:

            if iter_end == -1:
                iter_end = np.array(self.history[i]).flatten().shape[0]

            y = np.array(self.history[i]).flatten()[iter_start:iter_end]
            if log_y:
                y = np.log10(y)
            ax.plot(np.arange(iter_start, iter_end), y, label="train")
            ax.set_xlabel("Training epochs")
            ax.set_ylabel("Reconstruction accuracy (-ELBO loss)")
            ax.legend()
            plt.tight_layout()

    def export2adata(self, adata, slot_name="mod"):
        r"""Add posterior mean and sd for all parameters to unstructured data `adata.uns['mod']`.

        :param adata: anndata object
        """
        # add factor filter and samples of all parameters to unstructured data
        adata.uns[slot_name] = {}

        adata.uns[slot_name]["mod_name"] = str(self.module.__class__.__name__)
        adata.uns[slot_name]["fact_filt"] = self.fact_filt
        adata.uns[slot_name]["fact_names"] = self.fact_names.tolist()
        adata.uns[slot_name]["var_names"] = self.var_names.tolist()
        adata.uns[slot_name]["obs_names"] = self.obs_names.tolist()
        adata.uns[slot_name]["post_sample_means"] = self.samples["post_sample_means"]
        adata.uns[slot_name]["post_sample_sds"] = self.samples["post_sample_sds"]
        adata.uns[slot_name]["post_sample_q05"] = self.samples["post_sample_q05"]
        adata.uns[slot_name]["post_sample_q95"] = self.samples["post_sample_q95"]

        return adata
