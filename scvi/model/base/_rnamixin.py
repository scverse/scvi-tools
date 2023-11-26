import inspect
import logging
import warnings
from collections.abc import Iterable, Sequence
from functools import partial
from typing import Literal, Optional, Union

import numpy as np
import pandas as pd
import sparse
import torch
import torch.distributions as db
from anndata import AnnData
from pyro.distributions.util import deep_to

from scvi import REGISTRY_KEYS, settings
from scvi._types import Number
from scvi.distributions._utils import DistributionConcatenator, subset_distribution
from scvi.model._utils import _get_batch_code_from_category, scrna_raw_counts_properties
from scvi.module.base._decorators import _move_data_to_device
from scvi.utils import de_dsp, unsupported_if_adata_minified

from ._utils import (
    _de_core,
)

logger = logging.getLogger(__name__)


class RNASeqMixin:
    """General purpose methods for RNA-seq analysis."""

    def _get_transform_batch_gen_kwargs(self, batch):
        if "transform_batch" in inspect.signature(self.module.generative).parameters:
            return {"transform_batch": batch}
        else:
            raise NotImplementedError(
                "Transforming batches is not implemented for this model."
            )

    def _get_importance_weights(
        self,
        adata: Optional[AnnData],
        indices: Optional[Sequence[int]],
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
            truncated as described in :cite:p:`ionides2008`. In particular, the provided value
            is used to threshold importance weights as a way to reduce the variance of the estimator.
        n_mc_samples
            Number of Monte Carlo samples to use for estimating the importance weights, by default 500
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
        )
        mask = torch.tensor(anchor_cells)
        qz_anchor = subset_distribution(qz, mask, 0)  # n_anchors, n_latent
        log_qz = qz_anchor.log_prob(zs.unsqueeze(-2)).sum(
            dim=-1
        )  # n_samples, n_cells, n_anchors

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
        log_px_z = torch.cat(log_px_z, dim=-1)

        log_pz = log_pz.reshape(-1, 1)
        log_px_z = log_px_z.reshape(-1, len(anchor_cells))
        log_qz = log_qz.reshape(-1, len(anchor_cells))
        log_px = log_px.reshape(1, len(anchor_cells))

        importance_weight = torch.logsumexp(
            log_pz + log_px_z - log_px - torch.logsumexp(log_qz, 1, keepdims=True),
            dim=1,
        )
        if truncation:
            tau = torch.logsumexp(importance_weight, 0) - np.log(
                importance_weight.shape[0]
            )
            importance_weight = torch.clamp(importance_weight, min=tau)

        log_probs = importance_weight - torch.logsumexp(importance_weight, 0)
        return log_probs.exp().numpy()

    @torch.inference_mode()
    def get_normalized_expression(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        transform_batch: Optional[Sequence[Union[Number, str]]] = None,
        gene_list: Optional[Sequence[str]] = None,
        library_size: Union[float, Literal["latent"]] = 1,
        n_samples: int = 1,
        n_samples_overall: int = None,
        weights: Optional[Literal["uniform", "importance"]] = None,
        batch_size: Optional[int] = None,
        return_mean: bool = True,
        return_numpy: Optional[bool] = None,
        **importance_weighting_kwargs,
    ) -> Union[np.ndarray, pd.DataFrame]:
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
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`. DataFrame includes
            gene names as columns. If either `n_samples=1` or `return_mean=True`, defaults to `False`.
            Otherwise, it defaults to `True`.
        importance_weighting_kwargs
            Keyword arguments passed into :meth:`~scvi.model.base.RNASeqMixin._get_importance_weights`.

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
        if n_samples_overall is not None:
            assert n_samples == 1  # default value
            n_samples = n_samples_overall // len(indices) + 1
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )

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
        if library_size == "latent":
            generative_output_key = "mu"
            scaling = 1
        else:
            generative_output_key = "scale"
            scaling = library_size

        store_distributions = weights == "importance"
        if store_distributions and len(transform_batch) > 1:
            raise NotImplementedError(
                "Importance weights cannot be computed when expression levels are averaged across batches."
            )

        exprs = []
        zs = []
        qz_store = DistributionConcatenator()
        px_store = DistributionConcatenator()
        for tensors in scdl:
            per_batch_exprs = []
            for batch in transform_batch:
                generative_kwargs = self._get_transform_batch_gen_kwargs(batch)
                inference_kwargs = {"n_samples": n_samples}
                inference_outputs, generative_outputs = self.module.forward(
                    tensors=tensors,
                    inference_kwargs=inference_kwargs,
                    generative_kwargs=generative_kwargs,
                    compute_loss=False,
                )
                exp_ = getattr(generative_outputs["px"], generative_output_key)
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

        if return_numpy is None or return_numpy is False:
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
        adata: Optional[AnnData] = None,
        groupby: Optional[str] = None,
        group1: Optional[Iterable[str]] = None,
        group2: Optional[str] = None,
        idx1: Optional[Union[Sequence[int], Sequence[bool], str]] = None,
        idx2: Optional[Union[Sequence[int], Sequence[bool], str]] = None,
        mode: Literal["vanilla", "change"] = "change",
        delta: float = 0.25,
        batch_size: Optional[int] = None,
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: Optional[Iterable[str]] = None,
        batchid2: Optional[Iterable[str]] = None,
        fdr_target: float = 0.05,
        silent: bool = False,
        weights: Optional[Literal["uniform", "importance"]] = "uniform",
        filter_outlier_cells: bool = False,
        importance_weighting_kwargs: Optional[dict] = None,
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
            Whether to filter outlier cells with :meth:`~scvi.model.base.DifferentialComputation.filter_outlier_cells`.
        importance_weighting_kwargs
            Keyword arguments passed into :meth:`~scvi.model.base.RNASeqMixin._get_importance_weights`.
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
        representation_fn = (
            self.get_latent_representation if filter_outlier_cells else None
        )

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

    @torch.inference_mode()
    def posterior_predictive_sample(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        n_samples: int = 1,
        gene_list: Optional[Sequence[str]] = None,
        batch_size: Optional[int] = None,
    ) -> sparse.GCXS:
        r"""
        Generate observation samples from the posterior predictive distribution.

        The posterior predictive distribution is written as :math:`p(\hat{x} \mid x)`.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        n_samples
            Number of samples for each cell.
        gene_list
            Names of genes of interest.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        x_new : :py:class:`torch.Tensor`
            tensor with shape (n_cells, n_genes, n_samples)
        """
        if self.module.gene_likelihood not in ["zinb", "nb", "poisson"]:
            raise ValueError("Invalid gene_likelihood.")

        adata = self._validate_anndata(adata)

        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )

        if indices is None:
            indices = np.arange(adata.n_obs)

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = adata.var_names
            gene_mask = [True if gene in gene_list else False for gene in all_genes]

        x_new = []
        for tensors in scdl:
            samples = self.module.sample(
                tensors,
                n_samples=n_samples,
            )
            if gene_list is not None:
                samples = samples[:, gene_mask, ...]
            x_new.append(sparse.GCXS.from_numpy(samples.numpy()))

        x_new = sparse.concatenate(x_new)  # Shape (n_cells, n_genes, n_samples)

        return x_new

    @torch.inference_mode()
    def _get_denoised_samples(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        n_samples: int = 25,
        batch_size: int = 64,
        rna_size_factor: int = 1000,
        transform_batch: Optional[Sequence[int]] = None,
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

        Returns
        -------
        denoised_samples
        """
        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )

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

            # This gamma is really l*w using scVI manuscript notation
            p = rate / (rate + px_dispersion)
            r = px_dispersion
            l_train = torch.distributions.Gamma(r, (1 - p) / p).sample()
            data = l_train.cpu().numpy()
            # """
            # In numpy (shape, scale) => (concentration, rate), with scale = p /(1 - p)
            # rate = (1 - p) / p  # = 1/scale # used in pytorch
            # """
            data_loader_list += [data]

            if n_samples > 1:
                data_loader_list[-1] = np.transpose(data_loader_list[-1], (1, 2, 0))

        return np.concatenate(data_loader_list, axis=0)

    @torch.inference_mode()
    def get_feature_correlation_matrix(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        n_samples: int = 10,
        batch_size: int = 64,
        rna_size_factor: int = 1000,
        transform_batch: Optional[Sequence[Union[Number, str]]] = None,
        correlation_type: Literal["spearman", "pearson"] = "spearman",
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
        for b in transform_batch:
            denoised_data = self._get_denoised_samples(
                adata=adata,
                indices=indices,
                n_samples=n_samples,
                batch_size=batch_size,
                rna_size_factor=rna_size_factor,
                transform_batch=b,
            )
            flattened = np.zeros(
                (denoised_data.shape[0] * n_samples, denoised_data.shape[1])
            )
            for i in range(n_samples):
                if n_samples == 1:
                    flattened[
                        denoised_data.shape[0] * (i) : denoised_data.shape[0] * (i + 1)
                    ] = denoised_data[:, :]
                else:
                    flattened[
                        denoised_data.shape[0] * (i) : denoised_data.shape[0] * (i + 1)
                    ] = denoised_data[:, :, i]
            if correlation_type == "pearson":
                corr_matrix = np.corrcoef(flattened, rowvar=False)
            elif correlation_type == "spearman":
                corr_matrix, _ = spearmanr(flattened)
            else:
                raise ValueError(
                    "Unknown correlation type. Choose one of 'spearman', 'pearson'."
                )
            corr_mats.append(corr_matrix)
        corr_matrix = np.mean(np.stack(corr_mats), axis=0)
        var_names = adata.var_names
        return pd.DataFrame(corr_matrix, index=var_names, columns=var_names)

    @torch.inference_mode()
    def get_likelihood_parameters(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        n_samples: Optional[int] = 1,
        give_mean: Optional[bool] = False,
        batch_size: Optional[int] = None,
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
        """
        adata = self._validate_anndata(adata)

        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
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

            n_batch = px_rate.size(0) if n_samples == 1 else px_rate.size(1)

            px_r = px_r.cpu().numpy()
            if len(px_r.shape) == 1:
                dispersion_list += [np.repeat(px_r[np.newaxis, :], n_batch, axis=0)]
            else:
                dispersion_list += [px_r]
            mean_list += [px_rate.cpu().numpy()]
            if self.module.gene_likelihood == "zinb":
                dropout_list += [px_dropout.cpu().numpy()]
                dropout = np.concatenate(dropout_list, axis=-2)
        means = np.concatenate(mean_list, axis=-2)
        dispersions = np.concatenate(dispersion_list, axis=-2)

        if give_mean and n_samples > 1:
            if self.module.gene_likelihood == "zinb":
                dropout = dropout.mean(0)
            means = means.mean(0)
            dispersions = dispersions.mean(0)

        return_dict = {}
        return_dict["mean"] = means

        if self.module.gene_likelihood == "zinb":
            return_dict["dropout"] = dropout
            return_dict["dispersions"] = dispersions
        if self.module.gene_likelihood == "nb":
            return_dict["dispersions"] = dispersions

        return return_dict

    @torch.inference_mode()
    @unsupported_if_adata_minified
    def get_latent_library_size(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        give_mean: bool = True,
        batch_size: Optional[int] = None,
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
        """
        self._check_if_trained(warn=False)

        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
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
                        "for the library size. Set `give_mean` to False to use the observed library size instead."
                    )
                library = torch.distributions.LogNormal(ql.loc, ql.scale).mean
            libraries += [library.cpu()]
        return torch.cat(libraries).numpy()
