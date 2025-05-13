import logging
import warnings
from collections.abc import Sequence
from functools import partial

import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from pyro import infer

from scvi import settings
from scvi.model._utils import _get_batch_code_from_category, parse_device_args
from scvi.utils import track

logger = logging.getLogger(__name__)


class ResolVIPredictiveMixin:
    """Mixin class for generating samples from posterior distribution using infer.predictive."""

    @torch.inference_mode()
    def get_latent_representation(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        give_mean: bool = True,
        mc_samples: int = 1,  # consistency, noqa, pylint: disable=unused-argument
        batch_size: int | None = None,
        return_dist: bool = False,
    ) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
        """
        Return the latent representation for each cell.

        This is denoted as :math:`z` in RESOLVI.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        give_mean
            Give mean of distribution or sample from it.
        mc_samples
            For consistency with scVI, this parameter is ignored.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_dist
            Return the distribution parameters of the latent variables rather than their sampled
            values. If `True`, ignores `give_mean` and `mc_samples`.

        Returns
        -------
        Low-dimensional representation for each cell or a tuple containing its mean and variance.
        """
        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        latent = []
        latent_qzm = []
        latent_qzv = []

        _, _, device = parse_device_args(
            accelerator="auto",
            devices="auto",
            return_device="torch",
            validate_single_device=True,
        )

        for tensors in scdl:
            _, kwargs = self.module._get_fn_args_from_batch(tensors)
            kwargs = {k: v.to(device) if v is not None else v for k, v in kwargs.items()}

            if kwargs["cat_covs"] is not None and self.module.encode_covariates:
                categorical_input = list(torch.split(kwargs["cat_covs"], 1, dim=1))
            else:
                categorical_input = ()

            qz_m, qz_v, z = self.module.z_encoder(
                torch.log1p(kwargs["x"] / torch.mean(kwargs["x"], dim=1, keepdim=True)),
                kwargs["batch_index"],
                *categorical_input,
            )
            qz = torch.distributions.Normal(qz_m, qz_v.sqrt())
            if give_mean:
                z = qz.loc

            latent += [z.cpu()]
            latent_qzm += [qz.loc.cpu()]
            latent_qzv += [qz.scale.square().cpu()]
        return (
            (torch.cat(latent_qzm).numpy(), torch.cat(latent_qzv).numpy())
            if return_dist
            else torch.cat(latent).numpy()
        )

    @torch.inference_mode()
    def get_normalized_expression_importance(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        transform_batch: Sequence[int | str] | None = None,
        gene_list: Sequence[str] | None = None,
        library_size: float | None = 1,
        n_samples: int = 30,
        n_samples_overall: int = None,
        batch_size: int | None = None,
        weights: str | np.ndarray | None = None,
        return_mean: bool = True,
        return_numpy: bool | None = None,
    ) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
        adata = self._validate_anndata(adata)

        if indices is None:
            indices = np.arange(adata.n_obs)
        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        transform_batch = _get_batch_code_from_category(
            self.get_anndata_manager(adata, required=True), transform_batch
        )

        gene_mask = slice(None) if gene_list is None else adata.var_names.isin(gene_list)

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                warnings.warn(
                    "`return_numpy` must be `True` if `n_samples > 1` and `return_mean` "
                    "is`False`, returning an `np.ndarray`.",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )
            return_numpy = True

        exprs = []
        weighting = []

        _, _, device = parse_device_args(
            accelerator="auto",
            devices="auto",
            return_device="torch",
            validate_single_device=True,
        )

        for tensors in scdl:
            args, kwargs = self.module._get_fn_args_from_batch(tensors)
            kwargs = {k: v.to(device) if v is not None else v for k, v in kwargs.items()}
            model_now = partial(self.module.model_simplified, corrected_rate=True)
            importance_dist = infer.Importance(
                model_now, guide=self.module.guide.guide_simplified, num_samples=10 * n_samples
            )
            posterior = importance_dist.run(*args, **kwargs)
            marginal = infer.EmpiricalMarginal(posterior, sites=["mean_poisson", "px_scale"])
            samples = torch.cat([marginal().unsqueeze(1) for i in range(n_samples)], 1)
            log_weights = (
                torch.distributions.Poisson(samples[0, ...] + 1e-3)
                .log_prob(kwargs["x"].to(samples.device))
                .sum(-1)
            )
            log_weights = log_weights / kwargs["x"].to(samples.device).sum(-1)
            weighting.append(log_weights.reshape(-1).cpu())
            exprs.append(samples[1, ...].cpu())
        exprs = torch.cat(exprs, axis=1).numpy()
        if return_mean:
            exprs = exprs.mean(0)
        weighting = torch.cat(weighting, axis=0).numpy()
        if library_size is not None:
            exprs = library_size * exprs

        if n_samples_overall is not None:
            # Converts the 3d tensor to a 2d tensor
            exprs = exprs.reshape(-1, exprs.shape[-1])
            n_samples_ = exprs.shape[0]
            if (weights is None) or weights == "uniform":
                p = None
            else:
                weighting -= weighting.max()
                weighting = np.exp(weighting)
                p = weighting / weighting.sum(axis=0, keepdims=True)

            ind_ = np.random.choice(n_samples_, n_samples_overall, p=p, replace=True)
            exprs = exprs[ind_]

        if return_numpy is None or return_numpy is False:
            return pd.DataFrame(
                exprs,
                columns=adata.var_names[gene_mask],
                index=adata.obs_names[indices],
            )
        else:
            return exprs

    @torch.inference_mode()
    def get_normalized_expression(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        transform_batch: Sequence[int | str] | None = None,
        gene_list: Sequence[str] | None = None,
        library_size: float | None = 1,
        n_samples: int = 1,
        n_samples_overall: int = None,
        batch_size: int | None = None,
        return_mean: bool = True,
        return_numpy: bool | None = None,
        silent: bool = True,
    ) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
        r"""Returns the normalized (decoded) gene expression.

        This is denoted as :math:`\rho_n` in the scVI paper.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        transform_batch
            Batch to condition on.
            If transform_batch is:

            - None, then real observed batch is used.
            - int, then batch transform_batch is used.
        gene_list
            Return frequencies of expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
        library_size
            Scale the expression frequencies to a common library size.
            This allows gene expression levels to be interpreted on a common scale of relevant
            magnitude.
        n_samples
            Number of posterior samples to use for estimation.
        n_samples_overall
            Number of posterior samples to use for estimation. Overrides `n_samples`.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`. DataFrame
            includes gene names as columns. If either `n_samples=1` or `return_mean=True`, defaults
             to `False`. Otherwise, it defaults to `True`.
        %(de_silent)s

        Returns
        -------
        If `n_samples` is provided and `return_mean` is False,
        this method returns a 3d tensor of shape (n_samples, n_cells, n_genes).
        If `n_samples` is provided and `return_mean` is True, it returns a 2d tensor
        of shape (n_cells, n_genes).
        In this case, return type is :class:`~pandas.DataFrame` unless `return_numpy` is True.
        Otherwise, the method expects `n_samples_overall` to be provided and returns a 2d tensor
        of shape (n_samples_overall, n_genes).
        """
        adata = self._validate_anndata(adata)

        if indices is None:
            indices = np.arange(adata.n_obs)
        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        transform_batch = _get_batch_code_from_category(
            self.get_anndata_manager(adata, required=True), transform_batch
        )

        gene_mask = slice(None) if gene_list is None else adata.var_names.isin(gene_list)

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                warnings.warn(
                    "`return_numpy` must be `True` if `n_samples > 1` and `return_mean` "
                    "is`False`, returning an `np.ndarray`.",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )
            return_numpy = True

        exprs = []

        _, _, device = parse_device_args(
            accelerator="auto",
            devices="auto",
            return_device="torch",
            validate_single_device=True,
        )

        for tensors in scdl:
            per_batch_exprs = []
            for batch in track(transform_batch, disable=silent):
                _, kwargs = self.module._get_fn_args_from_batch(tensors)
                kwargs = {k: v.to(device) if v is not None else v for k, v in kwargs.items()}

                if kwargs["cat_covs"] is not None and self.module.encode_covariates:
                    categorical_input = list(torch.split(kwargs["cat_covs"], 1, dim=1))
                else:
                    categorical_input = ()

                qz_m, qz_v, _ = self.module.z_encoder(
                    torch.log1p(kwargs["x"] / torch.mean(kwargs["x"], dim=1, keepdim=True)),
                    kwargs["batch_index"],
                    *categorical_input,
                )
                z = torch.distributions.Normal(qz_m, qz_v.sqrt()).sample([n_samples])

                if kwargs["cat_covs"] is not None:
                    categorical_input = list(torch.split(kwargs["cat_covs"], 1, dim=1))
                else:
                    categorical_input = ()
                if batch is not None:
                    batch = torch.full_like(kwargs["batch_index"], batch)
                else:
                    batch = kwargs["batch_index"]

                px_scale, _, px_rate, _ = self.module.model.decoder(
                    self.module.model.dispersion, z, kwargs["library"], batch, *categorical_input
                )
                if library_size is not None:
                    exp_ = library_size * px_scale
                else:
                    exp_ = px_rate

                exp_ = exp_[..., gene_mask]
                per_batch_exprs.append(exp_[None].cpu())
            per_batch_exprs = torch.cat(per_batch_exprs, dim=0).mean(0).numpy()
            exprs.append(per_batch_exprs)

        exprs = np.concatenate(exprs, axis=1)
        if return_mean:
            exprs = exprs.mean(0)

        if n_samples_overall is not None:
            # Converts the 3d tensor to a 2d tensor
            exprs = exprs.reshape(-1, exprs.shape[-1])
            n_samples_ = exprs.shape[0]
            ind_ = np.random.choice(n_samples_, n_samples_overall, replace=True)
            exprs = exprs[ind_]

        if return_numpy is None or return_numpy is False:
            return pd.DataFrame(
                exprs,
                columns=adata.var_names[gene_mask],
                index=adata.obs_names[indices],
            )
        else:
            return exprs

    @torch.inference_mode()
    def get_neighbor_abundance(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        neighbor_key: str | None = None,
        n_samples: int = 1,
        n_samples_overall: int = None,
        batch_size: int | None = None,
        summary_frequency: int = 2,
        weights: str | None = None,
        return_mean: bool = True,
        return_numpy: bool | None = None,
        **kwargs,
    ) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
        r"""Returns the normalized (decoded) gene expression.

        This is denoted as :math:`\rho_n` in the scVI paper.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        neighbor_key
            Obsm key containing the spatial neighbors of each cell.
        n_samples
            Number of posterior samples to use for estimation.
        n_samples_overall
            Number of posterior samples to use for estimation. Overrides `n_samples`.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        summary_frequency
            Compute summary_fn after summary_frequency batches. Reduces memory footprint.
        weights
            Spatial weights for each neighbor.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`. DataFrame
            includes gene names as columns. If either `n_samples=1` or `return_mean=True`, defaults
             to `False`. Otherwise, it defaults to `True`.
        kwargs
            Additional keyword arguments that have no effect and only serve for compatibility.

        Returns
        -------
        If `n_samples` is provided and `return_mean` is False,
        this method returns a 3d tensor of shape (n_samples, n_cells, n_celltypes).
        If `n_samples` is provided and `return_mean` is True, it returns a 2d tensor
        of shape (n_cells, n_celltypes).
        In this case, return type is :class:`~pandas.DataFrame` unless `return_numpy` is True.
        Otherwise, the method expects `n_samples_overall` to be provided and returns a 2d tensor
        of shape (n_samples_overall, n_celltypes).
        """
        if adata:
            assert neighbor_key is not None, "Must provide `neighbor_key` if `adata` is provided."
        adata = self._validate_anndata(adata)
        if indices is None:
            indices = np.arange(adata.n_obs)
        if neighbor_key is None:
            neighbor_key = self.adata_manager.registry["field_registries"]["index_neighbor"][
                "data_registry"
            ]["attr_key"]
            neighbor_obsm = adata.obsm[neighbor_key]
        else:
            neighbor_obsm = adata.obsm[neighbor_key]
        n_neighbors = neighbor_obsm.shape[-1]

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                warnings.warn(
                    "`return_numpy` must be `True` if `n_samples > 1` and `return_mean` "
                    "is `False`, returning an `np.ndarray`.",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )
            return_numpy = True

        if batch_size is not None:
            if batch_size % n_neighbors != 0:
                raise ValueError("Batch size must be divisible by the number of neighbors.")
        batch_size = batch_size if batch_size is not None else n_neighbors * settings.batch_size
        indices_ = neighbor_obsm[indices].reshape(-1)
        dl = self._make_data_loader(
            adata=adata, indices=indices_, shuffle=False, batch_size=batch_size
        )

        sampled_prediction = self.sample_posterior(
            input_dl=dl,
            model=self.module.model_corrected,
            return_sites=["probs_prediction"],
            summary_frequency=summary_frequency,
            num_samples=n_samples,
            return_samples=True,
        )
        flat_neighbor_abundance_ = sampled_prediction["posterior_samples"]["probs_prediction"]
        neighbor_abundance_ = flat_neighbor_abundance_.reshape(
            n_samples, len(indices), n_neighbors, -1
        )
        neighbor_abundance = np.average(neighbor_abundance_, axis=-2, weights=weights)

        if return_mean:
            neighbor_abundance = np.mean(neighbor_abundance, axis=0)

        if n_samples_overall is not None:
            # Converts the 3d tensor to a 2d tensor
            neighbor_abundance = neighbor_abundance.reshape(-1, neighbor_abundance.shape[-1])
            n_samples_ = neighbor_abundance.shape[0]
            ind_ = np.random.choice(n_samples_, n_samples_overall, replace=True)
            neighbor_abundance = neighbor_abundance[ind_]

        if return_numpy is None or return_numpy is False:
            assert return_mean, "Only numpy output is supported when `return_mean` is False."
            n_labels = len(neighbor_abundance[-1])
            return pd.DataFrame(
                neighbor_abundance,
                columns=self._label_mapping[:n_labels],
                index=adata.obs_names[indices],
            )
        else:
            return neighbor_abundance
