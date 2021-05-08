from typing import Optional

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pyro
import torch
from pyro import poutine
from pyro.infer import SVI, Predictive
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

        Returns
        -------
        ELBO history in self.module.history_

        """

        args, kwargs = self.module.model._get_fn_args_full_data(self.adata)
        gpus, device = parse_use_gpu_arg(use_gpu)

        args = [a.to(device) for a in args]
        kwargs = {k: v.to(device) for k, v in kwargs.items()}
        self.to_device(device)

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

        Returns
        -------
        ELBO history in self.module.history_

        """

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
            self.module.history_ = trainer.logger.history
            self.history_ = trainer.logger.history
        except AttributeError:
            self.history_ = None

        self.module.is_trained_ = True
        self.is_trained_ = True

    def _optim_param(self, lr, autoencoding_lr, clip_norm, module_names=["encoder"]):
        # detect variables in autoencoding guide
        if autoencoding_lr is not None:
            all_param_list = np.array(list(self.module.state_dict().keys()))
            ind = np.zeros(len(all_param_list)).astype(bool)
            for n in module_names:
                ind = ind | np.array([n in i for i in all_param_list]).astype(bool)
            param_list = all_param_list[ind]
            print(all_param_list)
        else:
            param_list = []
        print(param_list)

        # create function which fetches different lr for autoencoding guide
        def optim_param(module_name, param_name):
            if module_name + "." + param_name in param_list:
                print(module_name + "." + param_name)
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
        i = 0
        for tensor_dict in train_dl:
            if i == 0:
                args, kwargs = self.module._get_fn_args_from_batch(tensor_dict)
                args = [a.to(device) for a in args]
                kwargs = {k: v.to(device) for k, v in kwargs.items()}

                global_means = self.module.guide.quantiles([q], *args, **kwargs)
                global_means = {
                    k: global_means[k].cpu().detach().numpy()
                    for k in global_means.keys()
                    if k not in self.module.model.list_obs_plate_vars()["sites"]
                }
            i += 1

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

        args, kwargs = self.module.model._get_fn_args_full_data(self.adata)
        gpus, device = parse_use_gpu_arg(use_gpu)
        args = [a.to(device) for a in args]
        kwargs = {k: v.to(device) for k, v in kwargs.items()}

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

    def _sample_node(self, node, num_samples_batch: int = 10):

        self.module.batch_size = self.adata.n_obs

        args, kwargs = self.module._get_fn_args_for_predictive(self.adata)

        if self.use_gpu is True:
            self.module.cuda()

        predictive = Predictive(
            self.module.model, guide=self.module.guide, num_samples=num_samples_batch
        )

        post_samples = {
            k: v.detach().cpu().numpy()
            for k, v in predictive(*args, **kwargs).items()
            if k == node
        }

        return post_samples[node]

    def sample_node(self, node, n_sampl_batches, num_samples_batch: int = 10, suff=""):

        # sample first batch
        self.samples[node + suff] = self._sample_node(
            node, num_samples_batch=num_samples_batch
        )

        for it in tqdm(range(n_sampl_batches - 1)):
            # sample remaining batches
            post_node = self._sample_node(node, num_samples_batch=num_samples_batch)

            # concatenate batches
            self.samples[node + suff] = np.concatenate(
                (self.samples[node + suff], post_node), axis=0
            )

        # compute mean across samples
        self.samples[node + suff] = self.samples[node + suff].mean(0)

    def _sample_all(self, num_samples_batch: int = 10):

        self.module.batch_size = self.adata.n_obs

        args, kwargs = self.module._get_fn_args_for_predictive(self.adata)

        if self.use_gpu is True:
            self.module.cuda()

        predictive = Predictive(
            self.module.model, guide=self.module.guide, num_samples=num_samples_batch
        )

        post_samples = {
            k: v.detach().cpu().numpy() for k, v in predictive(*args, **kwargs).items()
        }

        return post_samples

    def sample_all(self, n_sampl_batches, num_samples_batch: int = 10):

        self.adata.uns["mod"] = {}

        # sample first batch
        self.adata.uns["mod"]["post_samples"] = self._sample_all(
            num_samples_batch=num_samples_batch
        )

        for it in tqdm(range(n_sampl_batches - 1)):
            # sample remaining batches
            post_samples = self._sample_all(num_samples_batch=num_samples_batch)

            # concatenate batches
            self.adata.uns["mod"]["post_samples"] = {
                k: np.concatenate(
                    (self.adata.uns["mod"]["post_samples"][k], post_samples[k]), axis=0
                )
                for k in post_samples.keys()
            }

    def sample_posterior(
        self,
        node="all",
        n_samples: int = 1000,
        num_samples_batch: int = 10,
        save_samples=False,
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

        self.n_samples = n_samples
        self.n_sampl_batches = int(np.ceil(n_samples / num_samples_batch))
        self.num_samples_batch = num_samples_batch

        if node == "all":
            # Sample all parameters - might use a lot of GPU memory

            self.sample_all(
                self.n_sampl_batches, num_samples_batch=self.num_samples_batch
            )

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
                node, self.n_sampl_batches, batch_size=self.num_samples_batch, suff=""
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
