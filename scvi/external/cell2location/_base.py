from datetime import date

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyro
import torch
from pyro import poutine
from pyro.infer.autoguide import AutoNormal, init_to_mean
from scipy.sparse import issparse

from scvi import _CONSTANTS
from scvi.data._anndata import get_from_registry
from scvi.dataloaders import AnnDataLoader
from scvi.model._utils import parse_use_gpu_arg

from .autoguide import AutoGuideList, AutoNormalEncoder


class AutoGuideMixinModule:
    def _create_autoguide(
        self,
        model,
        amortised,
        encoder_kwargs,
        data_transform,
        single_encoder,
        init_loc_fn=init_to_mean,
    ):

        if not amortised:
            _guide = AutoNormal(
                model,
                init_loc_fn=init_to_mean,
                create_plates=model.create_plates,
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
            amortised_vars = model.list_obs_plate_vars()
            _guide = AutoGuideList(model, create_plates=model.create_plates)
            _guide.append(
                AutoNormal(
                    pyro.poutine.block(
                        model, hide=list(amortised_vars["sites"].keys())
                    ),
                    init_loc_fn=init_loc_fn,
                )
            )
            if isinstance(data_transform, np.ndarray):
                # add extra info about gene clusters to the network
                self.register_buffer(
                    "gene_clusters", torch.tensor(data_transform.astype("float32"))
                )
                n_in = model.n_vars + data_transform.shape[1]
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
                n_in = model.n_vars
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
                n_in = model.n_vars

            _guide.append(
                AutoNormalEncoder(
                    pyro.poutine.block(
                        model, expose=list(amortised_vars["sites"].keys())
                    ),
                    amortised_plate_sites=amortised_vars,
                    n_in=n_in,
                    n_hidden=n_hidden,
                    data_transform=data_transform,
                    encoder_kwargs=encoder_kwargs,
                    single_encoder=single_encoder,
                )
            )
            return _guide

    def _data_transform_clusters(self):
        def _data_transform(x):
            return torch.log1p(torch.cat([x, x @ self.gene_clusters], dim=1))

        return _data_transform

    def _data_transform_scale(self):
        def _data_transform(x):
            # return (x - self.var_mean) / self.var_std
            return x / self.var_std

        return _data_transform


class QuantileMixin:
    """
    Reimplementation of cell2location [Kleshchevnikov20]_ model. This mixin class provides methods for:

    - computing median and quantiles of the posterior distribution using both direct and amortised inference

    """

    def _optim_param(
        self,
        lr: float = 0.01,
        autoencoding_lr: float = None,
        clip_norm: float = 200,
        module_names: list = ["encoder", "hidden2locs", "hidden2scales"],
    ):
        # TODO implement custom training method that can use this function.
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

    @torch.no_grad()
    def _posterior_quantile_amortised(
        self, q: float = 0.5, batch_size: int = 2048, use_gpu: bool = True
    ):
        """
        Compute median of the posterior distribution of each parameter, separating local (minibatch) variable
        and global variables, which is necessary when performing amortised inference.

        Note for developers: requires model class method which lists observation/minibatch plate
        variables (self.module.model.list_obs_plate_vars()).

        Parameters
        ----------
        q
            quantile to compute
        batch_size
            number of observations per batch
        use_gpu
            Bool, use gpu?

        Returns
        -------
        dictionary {variable_name: posterior median}

        """

        gpus, device = parse_use_gpu_arg(use_gpu)

        self.module.eval()

        train_dl = AnnDataLoader(self.adata, shuffle=False, batch_size=batch_size)

        # sample local parameters
        i = 0
        for tensor_dict in train_dl:

            args, kwargs = self.module._get_fn_args_from_batch(tensor_dict)
            args = [a.to(device) for a in args]
            kwargs = {k: v.to(device) for k, v in kwargs.items()}
            self.to_device(device)

            if i == 0:

                means = self.module.guide.quantiles([q], *args, **kwargs)
                means = {
                    k: means[k].cpu().numpy()
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

                means_ = self.module.guide.quantiles([q], *args, **kwargs)
                means_ = {
                    k: means_[k].cpu().numpy()
                    for k in means_.keys()
                    if k
                    in list(self.module.model.list_obs_plate_vars()["sites"].keys())
                }
                means = {
                    k: np.concatenate(
                        [means[k], means_[k]], axis=list(obs_plate.values())[0]
                    )
                    for k in means.keys()
                }
            i += 1

        # sample global parameters
        tensor_dict = next(iter(train_dl))
        args, kwargs = self.module._get_fn_args_from_batch(tensor_dict)
        args = [a.to(device) for a in args]
        kwargs = {k: v.to(device) for k, v in kwargs.items()}
        self.to_device(device)

        global_means = self.module.guide.quantiles([q], *args, **kwargs)
        global_means = {
            k: global_means[k].cpu().numpy()
            for k in global_means.keys()
            if k not in list(self.module.model.list_obs_plate_vars()["sites"].keys())
        }

        for k in global_means.keys():
            means[k] = global_means[k]

        self.module.to(device)

        return means

    @torch.no_grad()
    def _posterior_quantile(
        self, q: float = 0.5, batch_size: int = 2048, use_gpu: bool = True
    ):
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

        train_dl = AnnDataLoader(self.adata, shuffle=False, batch_size=batch_size)
        # sample global parameters
        tensor_dict = next(iter(train_dl))
        args, kwargs = self.module._get_fn_args_from_batch(tensor_dict)
        args = [a.to(device) for a in args]
        kwargs = {k: v.to(device) for k, v in kwargs.items()}
        self.to_device(device)

        means = self.module.guide.quantiles([q], *args, **kwargs)
        means = {k: means[k].cpu().detach().numpy() for k in means.keys()}

        return means

    def posterior_quantile(
        self, q: float = 0.5, batch_size: int = 2048, use_gpu: bool = True
    ):
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
            return self._posterior_quantile_amortised(
                q=q, batch_size=batch_size, use_gpu=use_gpu
            )
        else:
            return self._posterior_quantile(q=q, batch_size=batch_size, use_gpu=use_gpu)


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
        plt.title("Reconstruction accuracy")
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
            iter_end = len(self.history_["train_loss_epoch"])

        ax.plot(
            self.history_["train_loss_epoch"].index[iter_start:iter_end],
            np.array(self.history_["train_loss_epoch"].values.flatten())[
                iter_start:iter_end
            ],
            label="train",
        )
        ax.legend()
        ax.xlim(0, len(self.history_["train_loss_epoch"]))
        ax.set_xlabel("Training epochs")
        ax.set_ylabel("-ELBO loss")
        plt.tight_layout()

    def _export2adata(self, samples):
        r"""
        Export key model variables and samples

        Parameters
        ----------
        samples
            dictionary with posterior mean, 5%/95% quantiles, SD, samples, generated by `.sample_posterior()`

        Returns
        -------
        updated dictionary with additional details is saved to `adata.uns['mod']`.
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
        summary_name: str = "means",
        name_prefix: str = "cell_abundance",
    ):
        """Export posterior distribution summary for observation-specific parameters
        (e.g. spatial cell abundance) as Pandas data frame
        (means, 5%/95% quantiles or sd of posterior distribution).

        Parameters
        ----------
        samples
            dictionary with posterior mean, 5%/95% quantiles, SD, samples, generated by `.sample_posterior()`
        site_name
            name of the model parameter to be exported
        summary_name
            posterior distribution summary to return ['means', 'sds', 'q05', 'q95']
        name_prefix
            prefix to add to column names (f'{summary_name}{name_prefix}_{site_name}_{self.factor_names_}')

        Returns
        -------
        Pandas data frame corresponding to either means, 5%/95% quantiles or sd of the posterior distribution

        """

        return pd.DataFrame(
            samples[f"post_sample_{summary_name}"].get(site_name, None),
            index=self.adata.obs_names,
            columns=[
                f"{summary_name}{name_prefix}_{site_name}_{i}"
                for i in self.factor_names_
            ],
        )

    def sample2df_vars(
        self,
        samples: dict,
        site_name: str = "gene_factors",
        summary_name: str = "means",
        name_prefix: str = "",
    ):
        """Export posterior distribution summary for variable-specific parameters as Pandas data frame
        (means, 5%/95% quantiles or sd of posterior distribution).

        Parameters
        ----------
        samples
            dictionary with posterior mean, 5%/95% quantiles, SD, samples, generated by `.sample_posterior()`
        site_name
            name of the model parameter to be exported
        summary_name
            posterior distribution summary to return ('means', 'sds', 'q05', 'q95')
        name_prefix
            prefix to add to column names (f'{summary_name}{name_prefix}_{site_name}_{self.factor_names_}')

        Returns
        -------
        Pandas data frame corresponding to either means, 5%/95% quantiles or sd of the posterior distribution

        """

        return pd.DataFrame(
            samples[f"post_sample_{summary_name}"].get(site_name, None),
            columns=self.adata.var_names,
            index=[
                f"{summary_name}{name_prefix}_{site_name}_{i}"
                for i in self.factor_names_
            ],
        ).T

    def plot_QC(self, summary_name: str = "means", use_n_obs: int = 1000):
        """
        Show quality control plots:
        1. Reconstruction accuracy to assess if there are any issues with model training.
            The plot should be roughly diagonal, strong deviations signal problems that need to be investigated.
            Plotting is slow because expected value of mRNA count needs to be computed from model parameters. Random
            observations are used to speed up computation.

        Parameters
        ----------
        summary_name
            posterior distribution summary to use ('means', 'sds', 'q05', 'q95')

        Returns
        -------

        """

        if getattr(self, "samples", False) is False:
            raise RuntimeError(
                "self.samples is missing, please run self.export_posterior() first"
            )
        if use_n_obs is not None:
            ind_x = np.random.choice(self.adata.n_obs, use_n_obs, replace=False)
        else:
            ind_x = None

        self.expected_nb_param = self.module.model.compute_expected(
            self.samples[f"post_sample_{summary_name}"], self.adata, ind_x=ind_x
        )
        x_data = get_from_registry(self.adata, _CONSTANTS.X_KEY)[ind_x, :]
        if issparse(x_data):
            x_data = np.asarray(x_data.toarray())
        self.plot_posterior_mu_vs_data(self.expected_nb_param["mu"], x_data)
