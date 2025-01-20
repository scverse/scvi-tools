from __future__ import annotations

import logging
import warnings
from collections import defaultdict
from functools import partial
from typing import TYPE_CHECKING

from torch import nn
from torch.distributions import Binomial

from scvi.distributions import BetaBinomial
from scvi.external.methylvi._utils import METHYLVI_REGISTRY_KEYS
from scvi.nn import FCLayers

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from typing import Literal

    from mudata import MuData

    from scvi._types import Number

import numpy as np
import pandas as pd
import torch

from scvi import settings
from scvi.model.base._de_core import (
    _de_core,
)

from ._utils import scmc_raw_counts_properties

logger = logging.getLogger(__name__)


class BSSeqMixin:
    """General purpose methods for BS-seq analysis."""

    @torch.inference_mode()
    def get_normalized_methylation(
        self,
        mdata: MuData | None = None,
        indices: Sequence[int] | None = None,
        region_list: Sequence[str] | None = None,
        n_samples: int = 1,
        n_samples_overall: int = None,
        batch_size: int | None = None,
        return_mean: bool = True,
        return_numpy: bool | None = None,
        context: str | None = None,
        **importance_weighting_kwargs,
    ) -> (np.ndarray | pd.DataFrame) | dict[str, np.ndarray | pd.DataFrame]:
        r"""Returns the normalized (decoded) methylation.

        This is denoted as :math:`\mu_n` in the methylVI paper.

        Parameters
        ----------
        mdata
            MuData object with equivalent structure to initial Mudata.
            If `None`, defaults to the MuData object used to initialize the model.
        indices
            Indices of cells in mdata to use. If `None`, all cells are used.
        region_list
            Return frequencies of expression for a subset of regions.
            This can save memory when working with large datasets and few regions are
            of interest.
        n_samples
            Number of posterior samples to use for estimation.
        n_samples_overall
            Number of posterior samples to use for estimation. Overrides `n_samples`.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`.
            DataFrame includes region names as columns. If either `n_samples=1` or
            `return_mean=True`, defaults to `False`. Otherwise, it defaults to `True`.
        context
            If not `None`, returns normalized methylation levels for the specified
            methylation context. Otherwise, a dictionary with contexts as keys and normalized
            methylation levels as values is returned.

        Returns
        -------
        If `n_samples` is provided and `return_mean` is False,
        this method returns a 3d tensor of shape (n_samples, n_cells, n_regions).
        If `n_samples` is provided and `return_mean` is True, it returns a 2d tensor
        of shape (n_cells, n_regions).
        In this case, return type is :class:`~pandas.DataFrame` unless `return_numpy` is True.
        Otherwise, the method expects `n_samples_overall` to be provided and returns a 2d tensor
        of shape (n_samples_overall, n_regions).

        If model was set up using a MuData object, a dictionary is returned with keys
        corresponding to individual methylation contexts with values determined as
        described above.
        """
        mdata = self._validate_anndata(mdata)

        if context is not None and context not in self.contexts:
            raise ValueError(
                f"{context} is not a valid methylation context for this model. "
                f"Valid contexts are {self.contexts}."
            )

        if indices is None:
            indices = np.arange(mdata.n_obs)
        if n_samples_overall is not None:
            assert n_samples == 1  # default value
            n_samples = n_samples_overall // len(indices) + 1
        scdl = self._make_data_loader(adata=mdata, indices=indices, batch_size=batch_size)

        region_mask = slice(None) if region_list is None else mdata.var_names.isin(region_list)

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                warnings.warn(
                    "`return_numpy` must be `True` if `n_samples > 1` and `return_mean` "
                    "is`False`, returning an `np.ndarray`.",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )
            return_numpy = True

        exprs = defaultdict(list)

        for tensors in scdl:
            inference_kwargs = {"n_samples": n_samples}
            inference_outputs, generative_outputs = self.module.forward(
                tensors=tensors,
                inference_kwargs=inference_kwargs,
                generative_kwargs={},
                compute_loss=False,
            )

            for ctxt in self.contexts:
                exp_ = generative_outputs["px_mu"][ctxt]
                exp_ = exp_[..., region_mask]
                exprs[ctxt].append(exp_.cpu())

        cell_axis = 1 if n_samples > 1 else 0

        for ctxt in self.contexts:
            exprs[ctxt] = np.concatenate(exprs[ctxt], axis=cell_axis)

        if n_samples_overall is not None:
            # Converts the 3d tensor to a 2d tensor
            for ctxt in self.contexts:
                exprs[ctxt] = exprs[ctxt].reshape(-1, exprs[ctxt].shape[-1])
                n_samples_ = exprs[ctxt].shape[0]
                ind_ = np.random.choice(n_samples_, n_samples_overall, replace=True)
                exprs[ctxt] = exprs[ctxt][ind_]
                return_numpy = True

        elif n_samples > 1 and return_mean:
            for ctxt in self.contexts:
                exprs[ctxt] = exprs[ctxt].mean(0)

        if return_numpy is None or return_numpy is False:
            exprs_dfs = {}
            for ctxt in self.contexts:
                exprs_dfs[ctxt] = pd.DataFrame(
                    exprs[ctxt],
                    columns=mdata[ctxt].var_names[region_mask],
                    index=mdata[ctxt].obs_names[indices],
                )
            exprs_ = exprs_dfs
        else:
            exprs_ = exprs

        if context is not None:
            return exprs_[context]
        else:
            return exprs_

    @torch.inference_mode()
    def get_specific_normalized_methylation(
        self,
        mdata: MuData | None = None,
        context: str = None,
        indices: Sequence[int] | None = None,
        transform_batch: Sequence[Number | str] | None = None,
        region_list: Sequence[str] | None = None,
        n_samples: int = 1,
        n_samples_overall: int = None,
        weights: Literal["uniform", "importance"] | None = None,
        batch_size: int | None = None,
        return_mean: bool = True,
        return_numpy: bool | None = None,
        **importance_weighting_kwargs,
    ) -> (np.ndarray | pd.DataFrame) | dict[str, np.ndarray | pd.DataFrame]:
        r"""Convenience function to obtain normalized methylation values for a single context.

        Parameters
        ----------
        mdata
            MuData object with equivalent structure to initial MuData. If `None`, defaults to the
            MuData object used to initialize the model.
        context
            Methylation context for which to obtain normalized methylation levels.
        indices
            Indices of cells in mdata to use. If `None`, all cells are used.
        transform_batch
            Batch to condition on.
            If transform_batch is:

            - None, then real observed batch is used.
            - int, then batch transform_batch is used.
        region_list
            Return frequencies of expression for a subset of regions.
            This can save memory when working with large datasets and few regions are
            of interest.
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
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`.
            DataFrame includes region names as columns. If either `n_samples=1` or
            `return_mean=True`, defaults to `False`. Otherwise, it defaults to `True`.
        importance_weighting_kwargs
            Keyword arguments passed into
            :meth:`~scvi.model.base.RNASeqMixin._get_importance_weights`.

        Returns
        -------
        If `n_samples` is provided and `return_mean` is False,
        this method returns a 3d tensor of shape (n_samples, n_cells, n_regions).
        If `n_samples` is provided and `return_mean` is True, it returns a 2d tensor
        of shape (n_cells, n_regions).
        In this case, return type is :class:`~pandas.DataFrame` unless `return_numpy` is True.
        Otherwise, the method expects `n_samples_overall` to be provided and returns a 2d tensor
        of shape (n_samples_overall, n_regions).
        """
        exprs = self.get_normalized_methylation(
            mdata=mdata,
            indices=indices,
            transform_batch=transform_batch,
            region_list=region_list,
            n_samples=n_samples,
            n_samples_overall=n_samples_overall,
            weights=weights,
            batch_size=batch_size,
            return_mean=return_mean,
            return_numpy=return_numpy,
            **importance_weighting_kwargs,
        )
        return exprs[context]

    def differential_methylation(
        self,
        mdata: MuData | None = None,
        groupby: str | None = None,
        group1: Iterable[str] | None = None,
        group2: str | None = None,
        idx1: Sequence[int] | Sequence[bool] | str | None = None,
        idx2: Sequence[int] | Sequence[bool] | str | None = None,
        mode: Literal["vanilla", "change"] = "vanilla",
        delta: float = 0.05,
        batch_size: int | None = None,
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: Iterable[str] | None = None,
        batchid2: Iterable[str] | None = None,
        fdr_target: float = 0.05,
        silent: bool = False,
        two_sided: bool = True,
        **kwargs,
    ) -> dict[str, pd.DataFrame] | pd.DataFrame:
        r"""\.

        A unified method for differential methylation analysis.

        Implements `"vanilla"` DE :cite:p:`Lopez18`. and `"change"` mode DE :cite:p:`Boyeau19`.

        Parameters
        ----------
        %(de_mdata)s
        %(de_modality)s
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
        two_sided
            Whether to perform a two-sided test, or a one-sided test.
        **kwargs
            Keyword args for :meth:`scvi.model.base.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential methylation DataFrame with the following columns:
        proba_de
            the probability of the region being differentially methylated
        is_de_fdr
            whether the region passes a multiple hypothesis correction procedure
            with the target_fdr threshold
        bayes_factor
            Bayes Factor indicating the level of significance of the analysis
        effect_size
            the effect size, computed as (accessibility in population 2) -
            (accessibility in population 1)
        emp_effect
            the empirical effect, based on observed detection rates instead of the estimated
            accessibility scores from the methylVI model
        scale1
            the estimated methylation level in population 1
        scale2
            the estimated methylation level in population 2
        emp_mean1
            the empirical (observed) methylation level in population 1
        emp_mean2
            the empirical (observed) methylation level in population 2

        """
        mdata = self._validate_anndata(mdata)

        def change_fn(a, b):
            return a - b

        if two_sided:

            def m1_domain_fn(samples):
                return np.abs(samples) >= delta

        else:

            def m1_domain_fn(samples):
                return samples >= delta

        result = {}
        for context in self.contexts:
            col_names = mdata[context].var_names
            model_fn = partial(
                self.get_specific_normalized_methylation,
                batch_size=batch_size,
                context=context,
            )
            all_stats_fn = partial(scmc_raw_counts_properties, context=context)

            result[context] = _de_core(
                adata_manager=self.get_anndata_manager(mdata, required=True),
                model_fn=model_fn,
                representation_fn=None,
                groupby=groupby,
                group1=group1,
                group2=group2,
                idx1=idx1,
                idx2=idx2,
                all_stats=all_stats,
                all_stats_fn=all_stats_fn,
                col_names=col_names,
                mode=mode,
                batchid1=batchid1,
                batchid2=batchid2,
                delta=delta,
                batch_correction=batch_correction,
                fdr=fdr_target,
                silent=silent,
                change_fn=change_fn,
                m1_domain_fn=m1_domain_fn,
                **kwargs,
            )

        return result


class BSSeqModuleMixin:
    """Shared methods for BS-seq VAE modules."""

    data_input_keys = [METHYLVI_REGISTRY_KEYS.MC_KEY, METHYLVI_REGISTRY_KEYS.COV_KEY]

    def _compute_minibatch_reconstruction_loss(self, minibatch_size, tensors, generative_outputs):
        reconst_loss = torch.zeros(minibatch_size).to(self.device)

        for context in self.contexts:
            px_mu = generative_outputs["px_mu"][context]
            px_gamma = generative_outputs["px_gamma"][context]
            mc = tensors[f"{context}_{METHYLVI_REGISTRY_KEYS.MC_KEY}"]
            cov = tensors[f"{context}_{METHYLVI_REGISTRY_KEYS.COV_KEY}"]

            if self.dispersion == "region":
                px_gamma = torch.sigmoid(self.px_gamma[context])

            if self.likelihood == "binomial":
                dist = Binomial(probs=px_mu, total_count=cov)
            elif self.likelihood == "betabinomial":
                dist = BetaBinomial(mu=px_mu, gamma=px_gamma, total_count=cov)

            reconst_loss += -dist.log_prob(mc).sum(dim=-1)

        return reconst_loss


class DecoderMETHYLVI(nn.Module):
    """Decodes data from latent space of ``n_input`` dimensions into ``n_output`` dimensions.

    Uses a fully-connected neural network of ``n_hidden`` layers.

    Parameters
    ----------
    n_input
        The dimensionality of the input (latent space)
    n_output
        The dimensionality of the output (data space)
    n_cat_list
        A list containing the number of categories
        for each category of interest. Each category will be
        included using a one-hot encoding
    n_layers
        The number of fully-connected hidden layers
    n_hidden
        The number of nodes per hidden layer
    dropout_rate
        Dropout rate to apply to each of the hidden layers
    inject_covariates
        Whether to inject covariates in each layer, or just the first (default).
    use_batch_norm
        Whether to use batch norm in layers
    use_layer_norm
        Whether to use layer norm in layers
    scale_activation
        Activation layer to use for px_scale_decoder
    **kwargs
        Keyword args for :class:`~scvi.nn.FCLayers`.
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_cat_list: Iterable[int] = None,
        n_layers: int = 1,
        n_hidden: int = 128,
        inject_covariates: bool = True,
        use_batch_norm: bool = False,
        use_layer_norm: bool = False,
        **kwargs,
    ):
        super().__init__()
        self.px_decoder = FCLayers(
            n_in=n_input,
            n_out=n_hidden,
            n_cat_list=n_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=0,
            inject_covariates=inject_covariates,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            **kwargs,
        )

        self.px_mu_decoder = nn.Sequential(
            nn.Linear(n_hidden, n_output),
            nn.Sigmoid(),
        )
        self.px_gamma_decoder = nn.Sequential(
            nn.Linear(n_hidden, n_output),
            nn.Sigmoid(),
        )

    def forward(
        self,
        dispersion: str,
        z: torch.Tensor,
        *cat_list: int,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """The forward computation for a single sample.

         #. Decodes the data from the latent space using the decoder network
         #. Returns parameters for the beta-binomial distribution of methylation
         #. If ``dispersion != 'region-cell'`` then value for that param will be ``None``

        Parameters
        ----------
        dispersion
            One of the following

            * ``'region'`` - dispersion parameter of NB is constant per region across cells
            * ``'region-cell'`` - dispersion can differ for every region in every cell
        z :
            tensor with shape ``(n_input,)``
        library_size
            library size
        cat_list
            list of category membership(s) for this sample

        Returns
        -------
        2-tuple of :py:class:`torch.Tensor`
            parameters for the Beta distribution of mean methylation values

        """
        px = self.px_decoder(z, *cat_list)
        px_mu = self.px_mu_decoder(px)
        px_gamma = self.px_gamma_decoder(px) if dispersion == "region-cell" else None

        return px_mu, px_gamma
