from typing import Optional

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyro
import torch
from anndata import AnnData
from pyro.infer import Predictive
from scipy.optimize import linear_sum_assignment
from tqdm.auto import tqdm

from scvi.data import register_tensor_from_anndata
from scvi.dataloaders import AnnDataLoader
from scvi.external.cell2location._module import Cell2locationModule
from scvi.lightning import PyroTrainingPlan, Trainer
from scvi.model.base import BaseModelClass


class Cell2locationBaseModelClass(BaseModelClass):
    """
    Reimplementation of cell2location [Kleshchevnikov20]_ model. BaseModelClass.

    https://github.com/BayraktarLab/cell2location

    Parameters
    ----------
    sc_adata
        single-cell AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    use_gpu
        Use the GPU or not.
    **model_kwargs
        Keyword args for :class:`~scvi.external.cell2location...`

    Examples
    --------
    >>>
    """

    def __init__(
        self,
        adata: AnnData,
        n_fact: int = 50,
        var_names_read=None,
        fact_names=None,
        sample_id=None,
        use_gpu: bool = True,
        batch_size: int = 1024,
        **model_kwargs,
    ):
        super(Cell2locationBaseModelClass, self).__init__(adata, use_gpu=use_gpu)

        adata.obs["_indices"] = np.arange(adata.n_obs).astype("int64")
        register_tensor_from_anndata(adata, "ind_x", "obs", "_indices")

        self.adata = adata

        self.n_var = self.summary_stats["n_vars"]
        self.n_obs = self.summary_stats["n_cells"]
        self.n_fact = n_fact
        self.batch_size = batch_size

        self.var_names = pd.Series(adata.var_names, index=adata.var_names)
        if var_names_read is None:
            self.var_names_read = pd.Series(self.var_names, index=self.var_names)
        else:
            self.var_names_read = pd.Series(
                adata.var[var_names_read], index=self.var_names
            )
        self.obs_names = pd.Series(adata.obs_names, index=adata.obs_names)

        if fact_names is None:
            self.fact_names = pd.Series(["fact_" + str(i) for i in range(self.n_fact)])
        else:
            self.fact_names = pd.Series(fact_names)

        if sample_id is None:
            self.sample_id = pd.Series(
                ["sample" for i in range(self.n_obs)], index=self.obs_names
            )
        else:
            self.sample_id = pd.Series(adata.obs[sample_id], index=self.obs_names)

        self.fact_filt = None

        # generate one-hot encoded matrix telling which obs belong to whic samples
        self.obs2sample_df = pd.get_dummies(self.sample_id)
        # convert to np.ndarray and register
        adata.obsm["obs2sample_obsm"] = self.obs2sample_df.values
        self.n_exper = self.obs2sample_df.shape[1]
        register_tensor_from_anndata(adata, "obs2sample", "obsm", "obs2sample_obsm")

    def train_pyro_v2(
        self,
        max_epochs: Optional[int] = None,
        use_gpu: Optional[bool] = None,
        train_size: float = 1,
        validation_size: Optional[float] = None,
        batch_size: Optional[int] = None,
        lr: float = 0.001,
        total_grad_norm_constraint: float = 200,
        **kwargs,
    ):

        if use_gpu is True:
            self.module.cuda()

        plan_kwargs = {
            "loss_fn": pyro.infer.Trace_ELBO(),  # for some reason JitTrace_ELBO does not work
            "optim": pyro.optim.ClippedAdam(
                {
                    "lr": lr,
                    # limit the gradient step from becoming too large
                    "clip_norm": total_grad_norm_constraint,
                }
            ),
        }

        if batch_size is None:
            batch_size = self.batch_size

        self.train(
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            plan_kwargs=plan_kwargs,
            plan_class=PyroTrainingPlan,
        )

    def train_pyro(
        self,
        max_epochs: int = 30000,
        lr: float = 0.001,
        total_grad_norm_constraint: float = 200,
        use_gpu: Optional[bool] = None,
        train_size: float = 1,
        validation_size: Optional[float] = None,
        batch_size=None,
        plan_kwargs: Optional[dict] = None,
        **kwargs,
    ):
        """
        Train the model using variational inference (Pyro).

        Parameters
        ----------
        max_epochs
            Number of epochs to train for
        lr
            Learning rate for optimization.
        use_gpu
            If `True`, use the GPU if available.
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        **kwargs
            Other keyword args for :class:`~scvi.lightning.Trainer`.
        """

        if use_gpu is None:
            use_gpu = self.use_gpu
        else:
            use_gpu = use_gpu and torch.cuda.is_available()
        gpus = 1 if use_gpu else None

        if self.use_gpu:
            torch.set_default_tensor_type(torch.cuda.FloatTensor)
        else:
            torch.set_default_tensor_type(torch.FloatTensor)

        if use_gpu:
            self.module.cuda()

        train_dl = AnnDataLoader(
            self.adata, shuffle=True, batch_size=self.batch_size, drop_last=True
        )
        pyro.clear_param_store()
        module = self.module
        # warmup guide for JIT
        for tensors in train_dl:
            args, kwargs = module._get_fn_args_from_batch(tensors)
            module.guide(*args, **kwargs)
            break

        train_dl = AnnDataLoader(
            self.adata, shuffle=True, batch_size=self.batch_size, drop_last=True
        )
        plan = PyroTrainingPlan(
            module,
            loss_fn=pyro.infer.Trace_ELBO(),
            # for some reason JitTrace_ELBO does not work (AssertionError)
            optim=pyro.optim.ClippedAdam(
                {
                    "lr": lr,
                    # limit the gradient step from becoming too large
                    "clip_norm": total_grad_norm_constraint,
                }
            ),
        )
        self.trainer = Trainer(
            gpus=gpus,
            max_epochs=max_epochs,
        )
        self.trainer.fit(plan, train_dl)

        # try:
        #    self.history = self.trainer.logger.history
        # except AttributeError:
        #    self.history = None

        self.is_trained_ = True

    def sample_node1(self, node, num_samples_batch: int = 10):

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
        self.samples[node + suff] = self.sample_node1(
            node, num_samples_batch=num_samples_batch
        )

        for it in tqdm(range(n_sampl_batches - 1)):
            # sample remaining batches
            post_node = self.sample_node1(node, num_samples_batch=num_samples_batch)

            # concatenate batches
            self.samples[node + suff] = np.concatenate(
                (self.samples[node + suff], post_node), axis=0
            )

        # compute mean across samples
        self.samples[node + suff] = self.samples[node + suff].mean(0)

    def sample_all1(self, num_samples_batch: int = 10):

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
        self.adata.uns["mod"]["post_samples"] = self.sample_all1(
            num_samples_batch=num_samples_batch
        )

        for it in tqdm(range(n_sampl_batches - 1)):
            # sample remaining batches
            post_samples = self.sample_all1(num_samples_batch=num_samples_batch)

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

    @property
    def _plan_class(self):
        return PyroTrainingPlan

    @property
    def _data_loader_cls(self):
        return AnnDataLoader

    @staticmethod
    def align_plot_stability(
        fac1, fac2, name1, name2, align=True, return_aligned=False
    ):
        r"""Align columns between two np.ndarrays using scipy.optimize.linear_sum_assignment,
            then plot correlations between columns in fac1 and fac2, ordering fac2 according to alignment

        :param fac1: np.ndarray 1, factors in columns
        :param fac2: np.ndarray 2, factors in columns
        :param name1: axis x name
        :param name2: axis y name
        :param align: boolean, match columns in fac1 and fac2 using linear_sum_assignment?
        """

        corr12 = np.corrcoef(fac1, fac2, False)
        ind_top = np.arange(0, fac1.shape[1])
        ind_right = np.arange(0, fac1.shape[1]) + fac1.shape[1]
        corr12 = corr12[ind_top, :][:, ind_right]
        corr12[np.isnan(corr12)] = -1

        if align:
            img = corr12[:, linear_sum_assignment(2 - corr12)[1]]
        else:
            img = corr12

        plt.imshow(img)

        plt.title(f"Training initialisation \n {name1} vs {name2}")
        plt.xlabel(name2)
        plt.ylabel(name1)

        plt.tight_layout()

        if return_aligned:
            return linear_sum_assignment(2 - corr12)[1]

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
            ax.set_ylabel("Reconstruction accuracy (ELBO loss)")
            ax.legend()
            plt.tight_layout()

    def export2adata(self, adata, slot_name="mod"):
        r"""Add posterior mean and sd for all parameters to unstructured data `adata.uns['mod']`.

        :param adata: anndata object
        """
        # add factor filter and samples of all parameters to unstructured data
        adata.uns[slot_name] = {}

        adata.uns[slot_name]["mod_name"] = str(self.__class__.__name__)
        adata.uns[slot_name]["fact_filt"] = self.fact_filt
        adata.uns[slot_name]["fact_names"] = self.fact_names.tolist()
        adata.uns[slot_name]["var_names"] = self.var_names.tolist()
        adata.uns[slot_name]["obs_names"] = self.obs_names.tolist()
        adata.uns[slot_name]["post_sample_means"] = self.samples["post_sample_means"]
        adata.uns[slot_name]["post_sample_sds"] = self.samples["post_sample_sds"]
        adata.uns[slot_name]["post_sample_q05"] = self.samples["post_sample_q05"]
        adata.uns[slot_name]["post_sample_q95"] = self.samples["post_sample_q95"]

        return adata


class Cell2locationModelClass(Cell2locationBaseModelClass):
    """
    Reimplementation of cell2location [Kleshchevnikov20]_ model. ModelClass.

    https://github.com/BayraktarLab/cell2location

    Parameters
    ----------
    sc_adata
        single-cell AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    use_gpu
        Use the GPU or not.
    **model_kwargs
        Keyword args for :class:`~scvi.external.cell2location...`

    Examples
    --------
    >>>
    """

    def __init__(
        self,
        adata: AnnData,
        cell_state_df: pd.DataFrame,
        var_names_read=None,
        sample_id=None,
        use_gpu: bool = True,
        batch_size: int = 1024,
        **model_kwargs,
    ):

        intersect = np.intersect1d(cell_state_df.index, adata.var_names)
        cell_state_df = cell_state_df.loc[intersect, :]
        adata = adata[:, intersect]
        adata.varm["cell_state_varm"] = cell_state_df.values

        super(Cell2locationModelClass, self).__init__(
            adata,
            n_fact=cell_state_df.shape[1],
            var_names_read=var_names_read,
            fact_names=cell_state_df.columns,
            sample_id=sample_id,
            use_gpu=use_gpu,
            batch_size=batch_size,
        )

        # register_tensor_from_anndata(adata, "cell_state", "varm", "cell_state_varm")
        self.cell_state_df = cell_state_df

    def sample2df(self, node_name="w_sf"):
        r"""Export cell locations as Pandas data frames.

        :param node_name: name of the model parameter to be exported
        :return: 4 Pandas dataframes added to model object:
            .spot_factors_df, .spot_factors_sd, .spot_factors_q05, .spot_factors_q95
        """

        if len(self.samples) == 0:
            raise ValueError(
                "Please run `.sample_posterior()` first to generate samples & summarise posterior of each parameter"
            )

        self.w_sf_df = pd.DataFrame.from_records(
            self.samples["post_sample_means"][node_name],
            index=self.obs_names,
            columns=["mean_" + node_name + i for i in self.fact_names],
        )

        self.w_sf_sd = pd.DataFrame.from_records(
            self.samples["post_sample_sds"][node_name],
            index=self.obs_names,
            columns=["sd_" + node_name + i for i in self.fact_names],
        )

        self.w_sf_q05 = pd.DataFrame.from_records(
            self.samples["post_sample_q05"][node_name],
            index=self.obs_names,
            columns=["q05_" + node_name + i for i in self.fact_names],
        )

        self.w_sf_q95 = pd.DataFrame.from_records(
            self.samples["post_sample_q95"][node_name],
            index=self.obs_names,
            columns=["q95_" + node_name + i for i in self.fact_names],
        )

    def annotate_spot_adata(self, adata):
        r"""Add cell locations to adata.obs

        :param adata: anndata object to annotate
        :return: updated anndata object
        """

        if self.spot_factors_df is None:
            self.sample2df()

        # add cell factors to adata
        adata.obs[self.w_sf_df.columns] = self.w_sf_df.loc[adata.obs.index, :]

        # add cell factor sd to adata
        adata.obs[self.w_sf_sd.columns] = self.w_sf_sd.loc[adata.obs.index, :]

        # add cell factor 5% and 95% quantiles to adata
        adata.obs[self.w_sf_q05.columns] = self.w_sf_q05.loc[adata.obs.index, :]
        adata.obs[self.w_sf_q95.columns] = self.w_sf_q95.loc[adata.obs.index, :]

        return adata


class Cell2location(Cell2locationModelClass):
    """
    Reimplementation of cell2location [Kleshchevnikov20]_ model. ModelClass.

    https://github.com/BayraktarLab/cell2location

    Parameters
    ----------
    sc_adata
        single-cell AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    use_gpu
        Use the GPU or not.
    **model_kwargs
        Keyword args for :class:`~scvi.external.cell2location...`

    Examples
    --------
    >>>
    """

    def __init__(
        self,
        adata: AnnData,
        cell_state_df: pd.DataFrame,
        var_names_read=None,
        sample_id=None,
        module=None,
        use_gpu: bool = True,
        batch_size: int = 2048,
        **model_kwargs,
    ):

        super(Cell2location, self).__init__(
            adata=adata,
            cell_state_df=cell_state_df,
            var_names_read=var_names_read,
            sample_id=sample_id,
            use_gpu=use_gpu,
            batch_size=batch_size,
        )

        if module is None:
            module = Cell2locationModule

        self.module = module(
            n_obs=self.n_obs,
            n_var=self.n_var,
            n_fact=self.n_fact,
            n_exper=self.n_exper,
            batch_size=self.batch_size,
            cell_state_mat=self.cell_state_df.values,
            **model_kwargs,
        )
