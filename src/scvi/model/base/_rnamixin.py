from __future__ import annotations

import inspect
import logging
import warnings
from functools import partial
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import torch
import torch.distributions as db
from pyro.distributions.util import deep_to

from scvi import REGISTRY_KEYS, settings
from scvi.data._utils import _validate_adata_dataloader_input
from scvi.distributions._utils import DistributionConcatenator, subset_distribution
from scvi.model._utils import _get_batch_code_from_category, scrna_raw_counts_properties
from scvi.model.base._de_core import _de_core
from scvi.module.base._decorators import _move_data_to_device
from scvi.utils import de_dsp, dependencies, track, unsupported_if_adata_minified

if TYPE_CHECKING:
    from collections.abc import Iterator
    from typing import Literal

    from anndata import AnnData
    from torch import Tensor

    from scvi._types import Number

try:
    from sparse import GCXS
except ImportError:
    GCXS = type(None)

logger = logging.getLogger(__name__)


class RNASeqMixin:
    """General purpose methods for RNA-seq analysis."""

    def _get_transform_batch_gen_kwargs(self, batch):
        if "transform_batch" in inspect.signature(self.module.generative).parameters:
            return {"transform_batch": batch}
        else:
            raise NotImplementedError("Transforming batches is not implemented for this model.")

    def _get_importance_weights(
        self,
        adata: AnnData | None,
        indices: list[int] | None,
        qz: db.Distribution,
        px: db.Distribution,
        zs: torch.Tensor,
        max_cells: int = 1024,
        truncation: bool = False,
        n_mc_samples: int = 500,
        n_mc_samples_per_pass: int = 250,
    ) -> np.ndarray:
        """Computes importance weights for the given samples.

        This method computes importance weights for every latent code in `zs` as a way to
        encourage latent codes providing high likelihoods across many cells in the considered
        subpopulation.

        Parameters
        ----------
        adata
            Data to use for computing importance weights.
        indices
            Indices of cells in adata to use.
        distributions
            Dictionary of distributions associated with `indices`.
        qz
            Variational posterior distributions of the cells, aligned with `indices`.
        px
            Count distributions of the cells, aligned with `indices`.
        zs
            Samples associated with `indices`.
        max_cells
            Maximum number of cells used to estimated the importance weights
        truncation
            Whether importance weights should be truncated. If True, the importance weights are
            truncated as described in :cite:p:`Ionides2008`. In particular, the provided value
            is used to threshold importance weights as a way to reduce the variance of the
            estimator.
        n_mc_samples
            Number of Monte Carlo samples to use for estimating the importance weights, by default
            500
        n_mc_samples_per_pass
            Number of Monte Carlo samples to use for each pass, by default 250

        Returns
        -------
        importance_weights
            Numpy array containing importance weights aligned with the provided `indices`.

        Notes
        -----
        This method assumes a normal prior on the latent space.
        """
        device = self.device
        log_pz = db.Normal(0, 1).log_prob(zs).sum(dim=-1)
        all_cell_indices = np.arange(len(indices))
        anchor_cells = (
            np.random.choice(all_cell_indices, size=max_cells, replace=False)
            if len(indices) > max_cells
            else all_cell_indices
        )

        log_px = self.get_marginal_ll(
            adata,
            indices=indices[anchor_cells],
            return_mean=False,
            n_mc_samples=n_mc_samples,
            n_mc_samples_per_pass=n_mc_samples_per_pass,
        )  # n_anchors
        mask = torch.tensor(anchor_cells)
        qz_anchor = subset_distribution(qz, mask, 0)  # n_anchors, n_latent
        log_qz = qz_anchor.log_prob(zs.unsqueeze(-2)).sum(dim=-1)  # n_samples, n_cells, n_anchors

        log_px_z = []
        distributions_px = deep_to(px, device=device)
        scdl_anchor = self._make_data_loader(
            adata=adata, indices=indices[anchor_cells], batch_size=1
        )
        for tensors_anchor in scdl_anchor:
            tensors_anchor = _move_data_to_device(tensors_anchor, device)
            x_anchor = tensors_anchor[REGISTRY_KEYS.X_KEY]  # 1, n_genes
            distributions_px.mu = distributions_px.scale * x_anchor.sum(-1)
            log_px_z.append(
                distributions_px.log_prob(x_anchor).sum(dim=-1)[..., None].cpu()
            )  # n_samples, n_cells, 1
        log_px_z = torch.cat(log_px_z, dim=-1)  # n_samples, n_cells, n_anchors

        log_pz = log_pz.reshape(-1, 1)
        log_px_z = log_px_z.reshape(-1, len(anchor_cells))
        log_qz = log_qz.reshape(-1, len(anchor_cells))
        log_px = log_px.reshape(1, len(anchor_cells))

        importance_weight = torch.logsumexp(
            log_pz + log_px_z - log_px - torch.logsumexp(log_qz, 1, keepdims=True),
            dim=1,
        )
        if truncation:
            tau = torch.logsumexp(importance_weight, 0) - np.log(importance_weight.shape[0])
            importance_weight = torch.clamp(importance_weight, min=tau)

        log_probs = importance_weight - torch.logsumexp(importance_weight, 0)
        return log_probs.exp().numpy()

    @torch.inference_mode()
    def get_normalized_expression(
        self,
        adata: AnnData | None = None,
        indices: list[int] | None = None,
        transform_batch: list[Number | str] | None = None,
        gene_list: list[str] | None = None,
        library_size: float | Literal["latent"] = 1,
        n_samples: int = 1,
        n_samples_overall: int = None,
        weights: Literal["uniform", "importance"] | None = None,
        batch_size: int | None = None,
        return_mean: bool = True,
        return_numpy: bool | None = None,
        silent: bool = True,
        dataloader: Iterator[dict[str, Tensor | None]] | None = None,
        **importance_weighting_kwargs,
    ) -> np.ndarray | pd.DataFrame:
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
            magnitude. If set to `"latent"`, use the latent library size.
        n_samples
            Number of posterior samples to use for estimation.
        n_samples_overall
            Number of posterior samples to use for estimation. Overrides `n_samples`.
        weights
            Weights to use for sampling. If `None`, defaults to `"uniform"`.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`. DataFrame
            includes gene names as columns. If either `n_samples=1` or `return_mean=True`, defaults
            to `False`. Otherwise, it defaults to `True`.
        %(de_silent)s
        dataloader
            An iterator over minibatches of data on which to compute the metric. The minibatches
            should be formatted as a dictionary of :class:`~torch.Tensor` with keys as expected by
            the model. If ``None``, a dataloader is created from ``adata``.
        importance_weighting_kwargs
            Keyword arguments passed into
            :meth:`~scvi.model.base.RNASeqMixin._get_importance_weights`.

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
        _validate_adata_dataloader_input(self, adata, dataloader)

        if dataloader is None:
            adata = self._validate_anndata(adata)

            if indices is None:
                indices = np.arange(adata.n_obs)
            if n_samples_overall is not None:
                assert n_samples == 1  # default value
                n_samples = n_samples_overall // len(indices) + 1
            scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

            transform_batch = _get_batch_code_from_category(
                self.get_anndata_manager(adata, required=True), transform_batch
            )

            gene_mask = slice(None) if gene_list is None else adata.var_names.isin(gene_list)

        else:
            scdl = dataloader
            for param in [indices, batch_size, n_samples]:
                if param is not None:
                    Warning(
                        f"Using {param} after custom Dataloader was initialize is redundant, "
                        f"please re-initialize with selected {param}",
                    )
            gene_mask = slice(None)
            transform_batch = [None]

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                warnings.warn(
                    "`return_numpy` must be `True` if `n_samples > 1` and `return_mean` "
                    "is`False`, returning an `np.ndarray`.",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )
            return_numpy = True
        if library_size == "latent":
            generative_output_key = "mu"
            scaling = 1
        else:
            generative_output_key = "scale"
            scaling = library_size

        store_distributions = weights == "importance"
        if store_distributions and len(transform_batch) > 1:
            raise NotImplementedError(
                "Importance weights cannot be computed when expression levels are averaged across "
                "batches."
            )

        exprs = []
        zs = []
        qz_store = DistributionConcatenator()
        px_store = DistributionConcatenator()
        for tensors in scdl:
            per_batch_exprs = []
            for batch in track(transform_batch, disable=silent):
                generative_kwargs = self._get_transform_batch_gen_kwargs(batch)
                inference_kwargs = {"n_samples": n_samples}
                inference_outputs, generative_outputs = self.module.forward(
                    tensors=tensors,
                    inference_kwargs=inference_kwargs,
                    generative_kwargs=generative_kwargs,
                    compute_loss=False,
                )
                px_generative = generative_outputs["px"]
                if isinstance(px_generative, torch.Tensor):
                    exp_ = px_generative
                else:
                    exp_ = px_generative.get_normalized(generative_output_key)
                exp_ = exp_[..., gene_mask]
                exp_ *= scaling
                per_batch_exprs.append(exp_[None].cpu())
                if store_distributions:
                    qz_store.store_distribution(inference_outputs["qz"])
                    px_store.store_distribution(generative_outputs["px"])

            zs.append(inference_outputs["z"].cpu())
            per_batch_exprs = torch.cat(per_batch_exprs, dim=0).mean(0).numpy()
            exprs.append(per_batch_exprs)

        cell_axis = 1 if n_samples > 1 else 0
        exprs = np.concatenate(exprs, axis=cell_axis)
        zs = torch.concat(zs, dim=cell_axis)

        if n_samples_overall is not None:
            # Converts the 3d tensor to a 2d tensor
            exprs = exprs.reshape(-1, exprs.shape[-1])
            n_samples_ = exprs.shape[0]
            if (weights is None) or weights == "uniform":
                p = None
            else:
                qz = qz_store.get_concatenated_distributions(axis=0)
                x_axis = 0 if n_samples == 1 else 1
                px = px_store.get_concatenated_distributions(axis=x_axis)
                p = self._get_importance_weights(
                    adata,
                    indices,
                    qz=qz,
                    px=px,
                    zs=zs,
                    **importance_weighting_kwargs,
                )
            ind_ = np.random.choice(n_samples_, n_samples_overall, p=p, replace=True)
            exprs = exprs[ind_]
        elif n_samples > 1 and return_mean:
            exprs = exprs.mean(0)

        if (return_numpy is None or return_numpy is False) and dataloader is None:
            return pd.DataFrame(
                exprs,
                columns=adata.var_names[gene_mask],
                index=adata.obs_names[indices],
            )
        else:
            return exprs

    @de_dsp.dedent
    def differential_expression(
        self,
        adata: AnnData | None = None,
        groupby: str | None = None,
        group1: list[str] | None = None,
        group2: str | None = None,
        idx1: list[int] | list[bool] | str | None = None,
        idx2: list[int] | list[bool] | str | None = None,
        mode: Literal["vanilla", "change"] = "vanilla",
        delta: float = 0.25,
        batch_size: int | None = None,
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: list[str] | None = None,
        batchid2: list[str] | None = None,
        fdr_target: float = 0.05,
        silent: bool = False,
        weights: Literal["uniform", "importance"] | None = "uniform",
        filter_outlier_cells: bool = False,
        importance_weighting_kwargs: dict | None = None,
        **kwargs,
    ) -> pd.DataFrame:
        r"""A unified method for differential expression analysis.

        Implements ``'vanilla'`` DE :cite:p:`Lopez18` and ``'change'`` mode DE :cite:p:`Boyeau19`.

        Parameters
        ----------
        %(de_adata)s
        %(de_groupby)s
        %(de_group1)s
        %(de_group2)s
        %(de_idx1)s
        %(de_idx2)s
        %(de_mode)s
        %(de_delta)s
        %(de_batch_size)s
        %(de_all_stats)s
        %(de_batch_correction)s
        %(de_batchid1)s
        %(de_batchid2)s
        %(de_fdr_target)s
        %(de_silent)s
        weights
            Weights to use for sampling. If `None`, defaults to `"uniform"`.
        filter_outlier_cells
            Whether to filter outlier cells with
            :meth:`~scvi.model.base.DifferentialComputation.filter_outlier_cells`.
        importance_weighting_kwargs
            Keyword arguments passed into
            :meth:`~scvi.model.base.RNASeqMixin._get_importance_weights`.
        **kwargs
            Keyword args for :meth:`scvi.model.base.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential expression DataFrame.
        """
        adata = self._validate_anndata(adata)
        col_names = adata.var_names
        importance_weighting_kwargs = importance_weighting_kwargs or {}
        model_fn = partial(
            self.get_normalized_expression,
            return_numpy=True,
            n_samples=1,
            batch_size=batch_size,
            weights=weights,
            **importance_weighting_kwargs,
        )
        representation_fn = self.get_latent_representation if filter_outlier_cells else None

        result = _de_core(
            self.get_anndata_manager(adata, required=True),
            model_fn,
            representation_fn,
            groupby,
            group1,
            group2,
            idx1,
            idx2,
            all_stats,
            scrna_raw_counts_properties,
            col_names,
            mode,
            batchid1,
            batchid2,
            delta,
            batch_correction,
            fdr_target,
            silent,
            **kwargs,
        )

        return result

    @dependencies("sparse")
    def posterior_predictive_sample(
        self,
        adata: AnnData | None = None,
        indices: list[int] | None = None,
        n_samples: int = 1,
        gene_list: list[str] | None = None,
        batch_size: int | None = None,
        dataloader: Iterator[dict[str, Tensor | None]] | None = None,
    ) -> GCXS:
        r"""Generate predictive samples from the posterior predictive distribution.

        The posterior predictive distribution is denoted as :math:`p(\hat{x} \mid x)`, where
        :math:`x` is the input data and :math:`\hat{x}` is the sampled data.

        We sample from this distribution by first sampling ``n_samples`` times from the posterior
        distribution :math:`q(z \mid x)` for a given observation, and then sampling from the
        likelihood :math:`p(\hat{x} \mid z)` for each of these.

        Parameters
        ----------
        adata
            :class:`~anndata.AnnData` object with an equivalent structure to the model's dataset.
            If ``None``, defaults to the :class:`~anndata.AnnData` object used to initialize the
            model.
        indices
            Indices of the observations in ``adata`` to use. If ``None``, defaults to all the
            observations.
        n_samples
            Number of Monte Carlo samples to draw from the posterior predictive distribution for
            each observation.
        gene_list
            Names of the genes to which to subset. If ``None``, defaults to all genes.
        batch_size
            Minibatch size to use for data loading and model inference. Defaults to
            ``scvi.settings.batch_size``. Passed into
            :meth:`~scvi.model.base.BaseModelClass._make_data_loader`.
        dataloader
            An iterator over minibatches of data on which to compute the metric. The minibatches
            should be formatted as a dictionary of :class:`~torch.Tensor` with keys as expected by
            the model. If ``None``, a dataloader is created from ``adata``.

        Returns
        -------
        Sparse multidimensional array of shape ``(n_obs, n_vars)`` if ``n_samples == 1``, else
        ``(n_obs, n_vars, n_samples)``.
        """
        import sparse

        _validate_adata_dataloader_input(self, adata, dataloader)

        if dataloader is None:
            adata = self._validate_anndata(adata)
            dataloader = self._make_data_loader(
                adata=adata, indices=indices, batch_size=batch_size
            )

            if gene_list is None:
                gene_mask = slice(None)
            else:
                gene_mask = [gene in gene_list for gene in adata.var_names]
                if not np.any(gene_mask):
                    raise ValueError(
                        "None of the provided genes in ``gene_list`` were detected in the data."
                    )
        else:
            for param in [indices, batch_size, gene_list]:
                if param is not None:
                    Warning(
                        f"Using {param} after custom Dataloader was initialize is redundant, "
                        f"please re-initialize with selected {param}",
                    )
            gene_mask = slice(None)

        x_hat = []
        for tensors in dataloader:
            # (batch_size, n_vars) if n_samples == 1, else (batch_size, n_vars, n_samples)
            samples = self.module.sample(tensors, n_samples=n_samples)[:, gene_mask]
            x_hat.append(sparse.GCXS.from_numpy(samples.numpy()))

        # (n_minibatches, batch_size, n_vars) -> (n_obs, n_vars) if n_samples == 1, else
        # (n_minibatches, batch_size, n_vars, n_samples) -> (n_obs, n_vars, n_samples)
        return sparse.concatenate(x_hat, axis=0)

    @torch.inference_mode()
    def _get_denoised_samples(
        self,
        adata: AnnData | None = None,
        indices: list[int] | None = None,
        n_samples: int = 25,
        batch_size: int = 64,
        rna_size_factor: int = 1000,
        transform_batch: list[int] | None = None,
        dataloader: Iterator[dict[str, Tensor | None]] | None = None,
    ) -> np.ndarray:
        """Return samples from an adjusted posterior predictive.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        n_samples
            Number of posterior samples to use for estimation.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        rna_size_factor
            size factor for RNA prior to sampling gamma distribution.
        transform_batch
            int of which batch to condition on for all cells.
        dataloader
            An iterator over minibatches of data on which to compute the metric. The minibatches
            should be formatted as a dictionary of :class:`~torch.Tensor` with keys as expected by
            the model. If ``None``, a dataloader is created from ``adata``.

        Returns
        -------
        denoised_samples
        """
        _validate_adata_dataloader_input(self, adata, dataloader)

        if dataloader is None:
            adata = self._validate_anndata(adata)
            scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        else:
            scdl = dataloader
            for param in [indices, batch_size, n_samples]:
                if param is not None:
                    Warning(
                        f"Using {param} after custom Dataloader was initialize is redundant, "
                        f"please re-initialize with selected {param}",
                    )
            transform_batch = None

        data_loader_list = []
        for tensors in scdl:
            x = tensors[REGISTRY_KEYS.X_KEY]
            generative_kwargs = self._get_transform_batch_gen_kwargs(transform_batch)
            inference_kwargs = {"n_samples": n_samples}
            _, generative_outputs = self.module.forward(
                tensors=tensors,
                inference_kwargs=inference_kwargs,
                generative_kwargs=generative_kwargs,
                compute_loss=False,
            )
            if "px" in generative_outputs:
                px_scale = generative_outputs["px"].scale
                px_r = generative_outputs["px"].theta
            else:
                px_scale = generative_outputs["px_scale"]
                px_r = generative_outputs["px_r"]
            device = px_r.device

            rate = rna_size_factor * px_scale
            if len(px_r.size()) == 2:
                px_dispersion = px_r
            else:
                px_dispersion = torch.ones_like(x).to(device) * px_r

            # This gamma is using scVI manuscript notation
            p = rate / (rate + px_dispersion)
            r = px_dispersion
            # TODO: NEED TORCH MPS FIX for 'aten::_standard_gamma'
            l_train = (
                torch.distributions.Gamma(r.to("cpu"), ((1 - p) / p).to("cpu")).sample()
                if device.type == "mps"
                else torch.distributions.Gamma(r, (1 - p) / p).sample().cpu()
            )
            data = l_train.numpy()
            data_loader_list += [data]

            if n_samples > 1:
                data_loader_list[-1] = np.transpose(data_loader_list[-1], (1, 2, 0))

        return np.concatenate(data_loader_list, axis=0)

    @torch.inference_mode()
    def get_feature_correlation_matrix(
        self,
        adata: AnnData | None = None,
        indices: list[int] | None = None,
        n_samples: int = 10,
        batch_size: int = 64,
        rna_size_factor: int = 1000,
        transform_batch: list[Number | str] | None = None,
        correlation_type: Literal["spearman", "pearson"] = "spearman",
        silent: bool = True,
    ) -> pd.DataFrame:
        """Generate gene-gene correlation matrix using scvi uncertainty and expression.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        n_samples
            Number of posterior samples to use for estimation.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        rna_size_factor
            size factor for RNA prior to sampling gamma distribution.
        transform_batch
            Batches to condition on.
            If transform_batch is:

            - None, then real observed batch is used.
            - int, then batch transform_batch is used.
            - list of int, then values are averaged over provided batches.
        correlation_type
            One of "pearson", "spearman".
        %(de_silent)s

        Returns
        -------
        Gene-gene correlation matrix
        """
        from scipy.stats import spearmanr

        adata = self._validate_anndata(adata)

        transform_batch = _get_batch_code_from_category(
            self.get_anndata_manager(adata, required=True), transform_batch
        )

        corr_mats = []
        for b in track(transform_batch, disable=silent):
            denoised_data = self._get_denoised_samples(
                adata=adata,
                indices=indices,
                n_samples=n_samples,
                batch_size=batch_size,
                rna_size_factor=rna_size_factor,
                transform_batch=b,
            )
            flattened = np.zeros((denoised_data.shape[0] * n_samples, denoised_data.shape[1]))
            for i in range(n_samples):
                if n_samples == 1:
                    flattened[denoised_data.shape[0] * (i) : denoised_data.shape[0] * (i + 1)] = (
                        denoised_data[:, :]
                    )
                else:
                    flattened[denoised_data.shape[0] * (i) : denoised_data.shape[0] * (i + 1)] = (
                        denoised_data[:, :, i]
                    )
            if correlation_type == "pearson":
                corr_matrix = np.corrcoef(flattened, rowvar=False)
            elif correlation_type == "spearman":
                corr_matrix, _ = spearmanr(flattened)
            else:
                raise ValueError("Unknown correlation type. Choose one of 'spearman', 'pearson'.")
            corr_mats.append(corr_matrix)
        corr_matrix = np.mean(np.stack(corr_mats), axis=0)
        var_names = adata.var_names
        return pd.DataFrame(corr_matrix, index=var_names, columns=var_names)

    @torch.inference_mode()
    def get_likelihood_parameters(
        self,
        adata: AnnData | None = None,
        indices: list[int] | None = None,
        n_samples: int | None = 1,
        give_mean: bool | None = False,
        batch_size: int | None = None,
        dataloader: Iterator[dict[str, Tensor | None]] | None = None,
    ) -> dict[str, np.ndarray]:
        r"""Estimates for the parameters of the likelihood :math:`p(x \mid z)`.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        n_samples
            Number of posterior samples to use for estimation.
        give_mean
            Return expected value of parameters or a samples
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        dataloader
            An iterator over minibatches of data on which to compute the metric. The minibatches
            should be formatted as a dictionary of :class:`~torch.Tensor` with keys as expected by
            the model. If ``None``, a dataloader is created from ``adata``.
        """
        _validate_adata_dataloader_input(self, adata, dataloader)

        if dataloader is None:
            adata = self._validate_anndata(adata)
            scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        else:
            scdl = dataloader
            for param in [indices, batch_size, n_samples]:
                if param is not None:
                    Warning(
                        f"Using {param} after custom Dataloader was initialize is redundant, "
                        f"please re-initialize with selected {param}",
                    )

        dropout_list = []
        mean_list = []
        dispersion_list = []
        for tensors in scdl:
            inference_kwargs = {"n_samples": n_samples}
            _, generative_outputs = self.module.forward(
                tensors=tensors,
                inference_kwargs=inference_kwargs,
                compute_loss=False,
            )
            px = generative_outputs["px"]
            px_r = px.theta
            px_rate = px.mu
            if self.module.gene_likelihood == "zinb":
                px_dropout = px.zi_probs
                dropout_list += [px_dropout.cpu().numpy()]
                dropout = np.concatenate(dropout_list, axis=-2)

            n_batch = px_rate.size(0) if n_samples == 1 else px_rate.size(1)
            if self.module.gene_likelihood != "poisson":
                px_r = px_r.cpu().numpy()
                if len(px_r.shape) == 1:
                    dispersion_list += [np.repeat(px_r[np.newaxis, :], n_batch, axis=0)]
                else:
                    dispersion_list += [px_r]
            mean_list += [px_rate.cpu().numpy()]

        means = np.concatenate(mean_list, axis=-2)
        dispersions = np.concatenate(dispersion_list, axis=-2)

        if give_mean and n_samples > 1:
            if self.module.gene_likelihood == "zinb":
                dropout = dropout.mean(0)
            if self.module.gene_likelihood != "poisson":
                dispersions = dispersions.mean(0)
            means = means.mean(0)

        return_dict = {}
        return_dict["mean"] = means

        if self.module.gene_likelihood == "zinb":
            return_dict["dropout"] = dropout
        if self.module.gene_likelihood != "poisson":
            return_dict["dispersions"] = dispersions

        return return_dict

    @torch.inference_mode()
    @unsupported_if_adata_minified
    def get_latent_library_size(
        self,
        adata: AnnData | None = None,
        indices: list[int] | None = None,
        give_mean: bool = True,
        batch_size: int | None = None,
        dataloader: Iterator[dict[str, Tensor | None]] | None = None,
    ) -> np.ndarray:
        r"""Returns the latent library size for each cell.

        This is denoted as :math:`\ell_n` in the scVI paper.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        give_mean
            Return the mean or a sample from the posterior distribution.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        dataloader
            An iterator over minibatches of data on which to compute the metric. The minibatches
            should be formatted as a dictionary of :class:`~torch.Tensor` with keys as expected by
            the model. If ``None``, a dataloader is created from ``adata``.
        """
        _validate_adata_dataloader_input(self, adata, dataloader)

        if dataloader is None:
            self._check_if_trained(warn=False)
            adata = self._validate_anndata(adata)
            scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        else:
            scdl = dataloader
            for param in [indices, batch_size]:
                if param is not None:
                    Warning(
                        f"Using {param} after custom Dataloader was initialize is redundant, "
                        f"please re-initialize with selected {param}",
                    )

        libraries = []
        for tensors in scdl:
            inference_inputs = self.module._get_inference_input(tensors)
            outputs = self.module.inference(**inference_inputs)

            library = outputs["library"]
            if not give_mean:
                library = torch.exp(library)
            else:
                ql = outputs["ql"]
                if ql is None:
                    raise RuntimeError(
                        "The module for this model does not compute the posterior distribution "
                        "for the library size. Set `give_mean` to False to use the observed "
                        "library size instead."
                    )
                library = torch.distributions.LogNormal(ql.loc, ql.scale).mean
            libraries += [library.cpu()]
        return torch.cat(libraries).numpy()
