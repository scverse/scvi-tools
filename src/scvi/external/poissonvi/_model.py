from __future__ import annotations

import logging
from functools import partial
from typing import TYPE_CHECKING

import numpy as np
import torch

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
)
from scvi.model import PEAKVI
from scvi.model._utils import _init_library_size, scatac_raw_counts_properties
from scvi.model.base import (
    RNASeqMixin,
)
from scvi.model.base._de_core import _de_core
from scvi.module import VAE
from scvi.utils import setup_anndata_dsp
from scvi.utils._docstrings import de_dsp

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from typing import Literal

    import pandas as pd
    from anndata import AnnData

logger = logging.getLogger(__name__)


class POISSONVI(PEAKVI, RNASeqMixin):
    """
    Peak Variational Inference using a Poisson distribution :cite:p:`Martens2023`.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.POISSONVI.setup_anndata`.
    n_hidden
        Number of nodes per hidden layer. If `None`, defaults to square root
        of number of regions.
    n_latent
        Dimensionality of the latent space. If `None`, defaults to square root
        of `n_hidden`.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    dropout_rate
        Dropout rate for neural networks
    latent_distribution
        One of

        * ``'normal'`` - Normal distribution (Default)
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)
    **model_kwargs
        Keyword args for :class:`~scvi.module.VAE`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.external.POISSINVI.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.external.POISSONVI(adata)
    >>> vae.train()

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/atac/PoissonVI`
    """

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int | None = None,
        n_latent: int | None = None,
        n_layers: int = 2,
        dropout_rate: float = 0.1,
        latent_distribution: Literal["normal", "ln"] = "normal",
        **model_kwargs,
    ):
        # need to pass these in to get the correct defaults for peakvi
        super().__init__(adata, n_hidden=n_hidden, n_latent=n_latent)

        n_batch = self.summary_stats.n_batch
        use_size_factor_key = REGISTRY_KEYS.SIZE_FACTOR_KEY in self.adata_manager.data_registry
        library_log_means, library_log_vars = None, None
        if use_size_factor_key is not None:
            library_log_means, library_log_vars = _init_library_size(self.adata_manager, n_batch)

        self._module_cls = VAE
        self.get_normalized_function_name = "get_normalized_accessibility"

        self.module = self._module_cls(
            n_input=self.summary_stats.n_vars,
            n_batch=n_batch,
            n_labels=self.summary_stats.n_labels,
            n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
            n_cats_per_cov=self.module.n_cats_per_cov,
            n_hidden=self.module.n_hidden,
            n_latent=self.module.n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion="gene",  # not needed here
            gene_likelihood="poisson",  # fixed value for now, but we could think of allowing nb
            latent_distribution=latent_distribution,
            use_size_factor_key=use_size_factor_key,
            library_log_means=library_log_means,
            library_log_vars=library_log_vars,
            use_batch_norm="none",
            use_layer_norm="both",
            # to be consitent with PEAKVI architecture
            extra_encoder_kwargs={"activation_fn": torch.nn.LeakyReLU},
            extra_decoder_kwargs={"activation_fn": torch.nn.LeakyReLU},
            **model_kwargs,
        )

        self._model_summary_string = (
            "PoissonVI Model with the following params: \nn_hidden: {}, n_latent: {}, "
            "n_layers: {}, dropout_rate: {}, peak_likelihood: {}, latent_distribution: {}"
        ).format(
            n_hidden,
            n_latent,
            n_layers,
            dropout_rate,
            "poisson",
            latent_distribution,
        )
        self.init_params_ = self._get_init_params(locals())

    @torch.inference_mode()
    def get_normalized_accessibility(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] = None,
        transform_batch: str | int | None = None,
        region_list: Sequence[str] | None = None,
        library_size: float | Literal["latent"] = 1,
        normalize_regions: bool = False,
        n_samples: int = 1,
        n_samples_overall: int = None,
        weights: Literal["uniform", "importance"] | None = None,
        batch_size: int = 128,
        return_mean: bool = True,
        return_numpy: bool = False,
        **importance_weighting_kwargs,
    ) -> pd.DataFrame | np.ndarray:
        """Returns the normalized accessibility matrix.

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
        region_list
            Return frequencies of accessibility for a subset of regions.
            This can save memory when working with large datasets and few regions are
            of interest.
        library_size
            Scale the accessibility frequencies to a common library size.
            This allows accessibility counts to be interpreted on a common scale of relevant
            magnitude. If set to `"latent"`, use the latent library size.
        normalize_regions
            Whether to reintroduce region factors to scale the normalized accessibility. This makes
            the estimates closer to the input, but removes the region-level bias correction. False
            by default.
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
            includes region names as columns. If either `n_samples=1` or `return_mean=True`,
            defaults to `False`. Otherwise, it defaults to `True`.
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
        # this is similar to PeakVI's region normalization where we ignore the factor that is
        # learnt per region
        if not normalize_regions:
            region_factors = self.module.decoder.px_scale_decoder[-2].bias
            # set region_factors (bias) to 0
            self.module.decoder.px_scale_decoder[-2].bias = torch.nn.Parameter(
                torch.zeros_like(region_factors)
            )
        accs = super().get_normalized_expression(
            adata=adata,
            indices=indices,
            transform_batch=transform_batch,
            gene_list=region_list,
            library_size=library_size,
            n_samples=n_samples,
            n_samples_overall=n_samples_overall,
            weights=weights,
            batch_size=batch_size,
            return_mean=return_mean,
            return_numpy=return_numpy,
            **importance_weighting_kwargs,
        )
        if not normalize_regions:
            # reset region_factors (bias)
            self.module.decoder.px_scale_decoder[-2].bias = torch.nn.Parameter(region_factors)
        return accs

    @torch.inference_mode()
    def get_region_factors(self):
        """Return region-specific factors. CPU/GPU dependent"""
        if self.device.type == "cpu":
            region_factors = self.module.decoder.px_scale_decoder[-2].bias.numpy()
        else:
            region_factors = self.module.decoder.px_scale_decoder[-2].bias.cpu().numpy()  # gpu
        if region_factors is None:
            raise RuntimeError("region factors were not included in this model")
        return region_factors

    @de_dsp.dedent
    def differential_accessibility(
        self,
        adata: AnnData | None = None,
        groupby: str | None = None,
        group1: Iterable[str] | None = None,
        group2: str | None = None,
        idx1: Sequence[int] | Sequence[bool] | str | None = None,
        idx2: Sequence[int] | Sequence[bool] | str | None = None,
        mode: Literal["vanilla", "change"] = "change",
        delta: float = 0.05,
        batch_size: int | None = None,
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: Iterable[str] | None = None,
        batchid2: Iterable[str] | None = None,
        fdr_target: float = 0.05,
        silent: bool = False,
        two_sided: bool = True,
        weights: Literal["uniform", "importance"] | None = "uniform",
        filter_outlier_cells: bool = False,
        importance_weighting_kwargs: dict | None = None,
        **kwargs,
    ) -> pd.DataFrame:
        r"""\.

        A unified method for differential accessibility analysis.

        Implements `"vanilla"` DE :cite:p:`Lopez18`. and `"change"` mode DE :cite:p:`Boyeau19`.

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
        two_sided
            Whether to perform a two-sided test, or a one-sided test.
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
        Differential accessibility DataFrame with the following columns:
        prob_da
            the probability of the region being differentially accessible
        is_da_fdr
            whether the region passes a multiple hypothesis correction procedure with the
            target_fdr threshold
        bayes_factor
            Bayes Factor indicating the level of significance of the analysis
        effect_size
            the effect size, computed as (accessibility in population 2) - (accessibility in
            population 1)
        emp_effect
            the empirical effect, based on observed detection rates instead of the estimated
            accessibility scores from the PeakVI model
        est_prob1
            the estimated probability of accessibility in population 1
        est_prob2
            the estimated probability of accessibility in population 2
        emp_prob1
            the empirical (observed) probability of accessibility in population 1
        emp_prob2
            the empirical (observed) probability of accessibility in population 2

        """
        adata = self._validate_anndata(adata)
        col_names = adata.var_names
        importance_weighting_kwargs = importance_weighting_kwargs or {}
        model_fn = partial(
            self.get_normalized_accessibility,
            return_numpy=True,
            n_samples=1,
            batch_size=batch_size,
            weights=weights,
            **importance_weighting_kwargs,
        )
        representation_fn = self.get_latent_representation if filter_outlier_cells else None

        if two_sided:

            def m1_domain_fn(samples):
                return np.abs(samples) >= delta

        else:

            def m1_domain_fn(samples):
                return samples >= delta

        result = _de_core(
            adata_manager=self.get_anndata_manager(adata, required=True),
            model_fn=model_fn,
            representation_fn=representation_fn,
            groupby=groupby,
            group1=group1,
            group2=group2,
            idx1=idx1,
            idx2=idx2,
            all_stats=all_stats,
            all_stats_fn=scatac_raw_counts_properties,
            col_names=col_names,
            mode=mode,
            batchid1=batchid1,
            batchid2=batchid2,
            delta=delta,
            batch_correction=batch_correction,
            fdr=fdr_target,
            m1_domain_fn=m1_domain_fn,
            silent=silent,
            **kwargs,
        )

        # change the column names to prob as done in PeakVI
        result = result.rename(
            columns={
                "emp_mean1": "emp_prob1",
                "emp_mean2": "emp_prob2",
            }
        )
        return result

    def differential_expression(
        self,
    ):
        # Refer to function differential_accessibility
        msg = (
            f"differential_expression is not implemented for {self.__class__.__name__}, please "
            f"use {self.__class__.__name__}.differential_accessibility"
        )
        raise NotImplementedError(msg)

        return None

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        batch_key: str | None = None,
        labels_key: str | None = None,
        size_factor_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        %(param_batch_key)s
        %(param_labels_key)s
        %(param_size_factor_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, check_fragment_counts=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
            NumericalObsField(REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
