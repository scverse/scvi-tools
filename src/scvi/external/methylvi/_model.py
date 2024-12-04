from __future__ import annotations

import logging
import warnings
from collections import defaultdict
from functools import partial
from typing import TYPE_CHECKING

from scvi.data._utils import get_anndata_attribute
from scvi.external.methylvi import METHYLVI_REGISTRY_KEYS

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from typing import Literal

    from anndata import AnnData
    from mudata import MuData

    from scvi._types import Number

import numpy as np
import pandas as pd
import sparse
import torch

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager, fields
from scvi.data._constants import _SETUP_ARGS_KEY
from scvi.dataloaders import SemiSupervisedDataSplitter
from scvi.external.methylvi._utils import _context_cov_key, _context_mc_key
from scvi.model._utils import (
    get_max_epochs_heuristic,
)
from scvi.model.base import (
    ArchesMixin,
    BaseModelClass,
    UnsupervisedTrainingMixin,
    VAEMixin,
)
from scvi.model.base._de_core import (
    _de_core,
)
from scvi.train import SemiSupervisedTrainingPlan, TrainRunner
from scvi.train._callbacks import SubSampleLabels
from scvi.utils import setup_anndata_dsp
from scvi.utils._docstrings import devices_dsp

from ._module import METHYLANVAE, METHYLVAE
from ._utils import scmc_raw_counts_properties

logger = logging.getLogger(__name__)


class METHYLVI(VAEMixin, UnsupervisedTrainingMixin, ArchesMixin, BaseModelClass):
    """
    Model class for methylVI :cite:p:`Weinberger2023a`

    Parameters
    ----------
    mdata
        MuData object that has been registered via :meth:`~scvi.external.METHYLVI.setup_mudata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    **model_kwargs
        Keyword args for :class:`~scvi.external.methylvi.METHYLVAE`

    Examples
    --------
    >>> mdata = mudata.read_h5mu(path_to_mudata)
    >>> MethylVI.setup_mudata(mdata, batch_key="batch")
    >>> vae = MethylVI(mdata)
    >>> vae.train()
    >>> mdata.obsm["X_methylVI"] = vae.get_latent_representation()
    """

    def __init__(
        self,
        mdata: MuData,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        **model_kwargs,
    ):
        super().__init__(mdata)

        n_batch = self.summary_stats.n_batch
        n_cats_per_cov = (
            self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY)[
                fields.CategoricalJointObsField.N_CATS_PER_KEY
            ]
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )

        self.contexts = self.get_anndata_manager(mdata, required=True).registry[_SETUP_ARGS_KEY][
            "methylation_contexts"
        ]
        self.num_features_per_context = [mdata[context].shape[1] for context in self.contexts]

        n_input = np.sum(self.num_features_per_context)

        self.module = METHYLVAE(
            n_input=n_input,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            n_batch=n_batch,
            n_cats_per_cov=n_cats_per_cov,
            contexts=self.contexts,
            num_features_per_context=self.num_features_per_context,
            **model_kwargs,
        )
        self._model_summary_string = (
            "Overwrite this attribute to get an informative representation for your model"
        )
        # necessary line to get params that will be used for saving/loading
        self.init_params_ = self._get_init_params(locals())

        logger.info("The model has been initialized")

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        **kwargs,
    ) -> AnnData | None:
        """
        %(summary)s.

        Parameters
        ----------
        %(param_adata)s

        Returns
        -------
        %(returns)s
        """
        raise NotImplementedError("METHYLVI must be used with a MuData object.")

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_mudata(
        cls,
        mdata: MuData,
        mc_layer: str,
        cov_layer: str,
        methylation_contexts: Iterable[str],
        batch_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        modalities=None,
        **kwargs,
    ):
        """%(summary_mdata)s.

        Parameters
        ----------
        %(param_mdata)s
        mc_layer
            Layer containing methylated cytosine counts for each set of methylation features.
        cov_layer
            Layer containing total coverage counts for each set of methylation features.
        methylation_contexts
            List of modality fields in `mdata` object representing different methylation contexts.
            Each context must be equipped with a layer containing the number of methylated counts
            (specified by `mc_layer`) and total number of counts (specified by `cov_layer`) for
            each genomic region feature.
        %(param_batch_key)s
        %(param_categorical_covariate_keys)s
        %(param_modalities)s

        Examples
        --------
        MethylVI.setup_mudata(
            mdata,
            mc_layer="mc",
            cov_layer="cov",
            batch_key="Platform",
            methylation_modalities=['mCG', 'mCH'],
            modalities={
                "batch_key": "mCG"
            },
        )

        """
        if modalities is None:
            modalities = {}
        setup_method_args = METHYLVI._get_setup_method_args(**locals())

        if methylation_contexts is None:
            raise ValueError("Methylation contexts cannot be None.")

        modalities_ = cls._create_modalities_attr_dict(modalities, setup_method_args)

        batch_field = fields.MuDataCategoricalObsField(
            REGISTRY_KEYS.BATCH_KEY,
            batch_key,
            mod_key=modalities_.batch_key,
        )

        cat_cov_field = fields.MuDataCategoricalJointObsField(
            REGISTRY_KEYS.CAT_COVS_KEY,
            categorical_covariate_keys,
            mod_key=modalities_.categorical_covariate_keys,
        )

        mc_fields = []
        cov_fields = []

        for context in methylation_contexts:
            mc_fields.append(
                fields.MuDataLayerField(
                    _context_mc_key(context),
                    mc_layer,
                    mod_key=context,
                    is_count_data=True,
                    mod_required=True,
                )
            )

            cov_fields.append(
                fields.MuDataLayerField(
                    _context_cov_key(context),
                    cov_layer,
                    mod_key=context,
                    is_count_data=True,
                    mod_required=True,
                )
            )

        mudata_fields = mc_fields + cov_fields + [batch_field] + [cat_cov_field]
        adata_manager = AnnDataManager(fields=mudata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(mdata, **kwargs)

        cls.register_manager(adata_manager)

    @torch.inference_mode()
    def posterior_predictive_sample(
        self,
        mdata: MuData | None = None,
        n_samples: int = 1,
        batch_size: int | None = None,
    ) -> dict[str, sparse.GCXS] | sparse.GCXS:
        r"""
        Generate observation samples from the posterior predictive distribution.

        The posterior predictive distribution is written as :math:`p(\hat{x} \mid x)`.

        Parameters
        ----------
        mdata
            MuData object with equivalent structure to initial MuData. If `None`, defaults to the
            MuData object used to initialize the model.
        n_samples
            Number of samples for each cell.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        x_new : :py:class:`torch.Tensor`
            tensor with shape (n_cells, n_regions, n_samples)
        """
        mdata = self._validate_anndata(mdata)

        scdl = self._make_data_loader(adata=mdata, batch_size=batch_size)

        x_new = defaultdict(list)
        for tensors in scdl:
            samples = self.module.sample(
                tensors,
                n_samples=n_samples,
            )

            for context in self.contexts:
                x_new[context].append(sparse.GCXS.from_numpy(samples[context].numpy()))

        for context in self.contexts:
            x_new[context] = sparse.concatenate(
                x_new[context]
            )  # Shape (n_cells, n_regions, n_samples)

        return x_new

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

        Only applicable to MuData models.

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


class METHYLANVI(VAEMixin, ArchesMixin, BaseModelClass):
    """Methylation annotation using variational inference :cite:p:`Weinberger23`.

    Inspired from M1 + M2 model, as described in (https://arxiv.org/pdf/1406.5298.pdf).

    Parameters
    ----------
    mdata
        MuData object registered via :meth:`~scvi.external.methylvi.METHYLVI.setup_mudata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    dropout_rate
        Dropout rate for neural networks.
    likelihood
        One of
        * ``'betabinomial'`` - BetaBinomial distribution
        * ``'binomial'`` - Binomial distribution
    dispersion
        One of the following
        * ``'region'`` - dispersion parameter of BetaBinomial is constant per region across cells
        * ``'region-cell'`` - dispersion can differ for every region in every cell
    linear_classifier
        If ``True``, uses a single linear layer for classification instead of a
        multi-layer perceptron.
    **model_kwargs
        Keyword args for :class:`~scvi.module.SCANVAE`

    Examples
    --------
    >>> mdata = mudata.read_h5mu(path_to_mudata)
    >>> scvi.external.methylvi.METHYLANVI.setup_mudata(
    ...     mdata, labels_key="labels", unlabeled_category="Unknown"
    ... )
    >>> vae = scvi.external.methylvi.METHYLANVI(mdata)
    >>> vae.train()
    >>> mdata.obsm["X_scVI"] = vae.get_latent_representation()
    >>> mdata.obs["pred_label"] = vae.predict()

    """

    _module_cls = METHYLANVAE
    _training_plan_cls = SemiSupervisedTrainingPlan

    def __init__(
        self,
        mdata: MuData,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        likelihood: Literal["betabinomial", "binomial"] = "betabinomial",
        dispersion: Literal["region", "region-cell"] = "region",
        linear_classifier: bool = False,
        **model_kwargs,
    ):
        super().__init__(mdata)
        methylanvae_model_kwargs = dict(model_kwargs)

        self._set_indices_and_labels()

        # ignores unlabeled category
        n_labels = self.summary_stats.n_labels - 1
        n_cats_per_cov = (
            self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )

        n_batch = self.summary_stats.n_batch

        self.contexts = self.get_anndata_manager(mdata, required=True).registry[_SETUP_ARGS_KEY][
            "methylation_contexts"
        ]
        self.num_features_per_context = [mdata[context].shape[1] for context in self.contexts]

        n_input = np.sum(self.num_features_per_context)

        self.module = self._module_cls(
            n_input=n_input,
            n_batch=n_batch,
            n_cats_per_cov=n_cats_per_cov,
            n_labels=n_labels,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            likelihood=likelihood,
            linear_classifier=linear_classifier,
            contexts=self.contexts,
            num_features_per_context=self.num_features_per_context,
            **methylanvae_model_kwargs,
        )

        self.unsupervised_history_ = None
        self.semisupervised_history_ = None

        self._model_summary_string = (
            f"MethylANVI Model with the following params: \nunlabeled_category: "
            f"{self.unlabeled_category_}, n_hidden: {n_hidden}, n_latent: {n_latent}"
            f", n_layers: {n_layers}, dropout_rate: {dropout_rate}, dispersion: "
            f"{dispersion}, likelihood: {likelihood}"
        )
        self.init_params_ = self._get_init_params(locals())
        self.was_pretrained = False
        self.n_labels = n_labels

    def _set_indices_and_labels(self):
        """Set indices for labeled and unlabeled cells."""
        labels_state_registry = self.adata_manager.get_state_registry(REGISTRY_KEYS.LABELS_KEY)
        self.original_label_key = labels_state_registry.original_key
        self.unlabeled_category_ = labels_state_registry.unlabeled_category

        labels = get_anndata_attribute(
            self.adata,
            self.adata_manager.data_registry.labels.attr_name,
            self.original_label_key,
        ).ravel()
        self._label_mapping = labels_state_registry.categorical_mapping

        # set unlabeled and labeled indices
        self._unlabeled_indices = np.argwhere(labels == self.unlabeled_category_).ravel()
        self._labeled_indices = np.argwhere(labels != self.unlabeled_category_).ravel()
        self._code_to_label = dict(enumerate(self._label_mapping))

    def predict(
        self,
        mdata: MuData | None = None,
        indices: Sequence[int] | None = None,
        soft: bool = False,
        batch_size: int | None = None,
        use_posterior_mean: bool = True,
    ) -> np.ndarray | pd.DataFrame:
        """Return cell label predictions.

        Parameters
        ----------
        mdata
            MuData object registered via :meth:`~scvi.external.methylvi.MethylANVI.setup_mudata`.
        indices
            Return probabilities for each class label.
        soft
            If True, returns per class probabilities
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        use_posterior_mean
            If ``True``, uses the mean of the posterior distribution to predict celltype
            labels. Otherwise, uses a sample from the posterior distribution - this
            means that the predictions will be stochastic.
        """
        mdata = self._validate_anndata(mdata)

        if indices is None:
            indices = np.arange(mdata.n_obs)

        scdl = self._make_data_loader(
            adata=mdata,
            indices=indices,
            batch_size=batch_size,
        )
        y_pred = []
        for _, tensors in enumerate(scdl):
            inference_inputs = self.module._get_inference_input(tensors)  # (n_obs, n_vars)

            mc = inference_inputs[METHYLVI_REGISTRY_KEYS.MC_KEY]
            cov = inference_inputs[METHYLVI_REGISTRY_KEYS.COV_KEY]
            batch = tensors[REGISTRY_KEYS.BATCH_KEY]

            pred = self.module.classify(
                mc,
                cov,
                batch_index=batch,
                use_posterior_mean=use_posterior_mean,
            )
            if self.module.classifier.logits:
                pred = torch.nn.functional.softmax(pred, dim=-1)
            if not soft:
                pred = pred.argmax(dim=1)
            y_pred.append(pred.detach().cpu())

        y_pred = torch.cat(y_pred).numpy()
        if not soft:
            predictions = []
            for p in y_pred:
                predictions.append(self._code_to_label[p])

            return np.array(predictions)
        else:
            n_labels = len(pred[0])
            pred = pd.DataFrame(
                y_pred,
                columns=self._label_mapping[:n_labels],
                index=mdata.obs_names[indices],
            )
            return pred

    @devices_dsp.dedent
    def train(
        self,
        max_epochs: int | None = None,
        n_samples_per_label: float | None = None,
        check_val_every_n_epoch: int | None = None,
        train_size: float = 0.9,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        batch_size: int = 128,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        **trainer_kwargs,
    ):
        """Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset for semisupervised training.
        n_samples_per_label
            Number of subsamples for each label class to sample per epoch. By default, there
            is no label subsampling.
        check_val_every_n_epoch
            Frequency with which metrics are computed on the data for validation set for both
            the unsupervised and semisupervised trainers. If you'd like a different frequency for
            the semisupervised trainer, set check_val_every_n_epoch in semisupervised_train_kwargs.
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        shuffle_set_split
            Whether to shuffle indices before splitting. If `False`, the val, train,
            and test set are split in the sequential order of the data according to
            `validation_size` and `train_size` percentages.
        batch_size
            Minibatch size to use during training.
        %(param_accelerator)s
        %(param_devices)s
        datasplitter_kwargs
            Additional keyword arguments passed into
            :class:`~scvi.dataloaders.SemiSupervisedDataSplitter`.
        plan_kwargs
            Keyword args for :class:`~scvi.train.SemiSupervisedTrainingPlan`. Keyword
            arguments passed to `train()` will overwrite values present in `plan_kwargs`,
            when appropriate.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        if max_epochs is None:
            max_epochs = get_max_epochs_heuristic(self.adata.n_obs)

            if self.was_pretrained:
                max_epochs = int(np.min([10, np.max([2, round(max_epochs / 3.0)])]))

        logger.info(f"Training for {max_epochs} epochs.")

        plan_kwargs = {} if plan_kwargs is None else plan_kwargs
        datasplitter_kwargs = datasplitter_kwargs or {}

        # if we have labeled cells, we want to subsample labels each epoch
        sampler_callback = [SubSampleLabels()] if len(self._labeled_indices) != 0 else []

        data_splitter = SemiSupervisedDataSplitter(
            adata_manager=self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            n_samples_per_label=n_samples_per_label,
            batch_size=batch_size,
            **datasplitter_kwargs,
        )
        training_plan = self._training_plan_cls(self.module, self.n_labels, **plan_kwargs)

        if "callbacks" in trainer_kwargs.keys():
            trainer_kwargs["callbacks"] + [sampler_callback]
        else:
            trainer_kwargs["callbacks"] = sampler_callback

        runner = TrainRunner(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=devices,
            check_val_every_n_epoch=check_val_every_n_epoch,
            **trainer_kwargs,
        )
        return runner()

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        **kwargs,
    ) -> AnnData | None:
        """
        %(summary)s.

        Parameters
        ----------
        %(param_adata)s

        Returns
        -------
        %(returns)s
        """
        raise NotImplementedError("METHYLANVI must be used with a MuData object.")

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_mudata(
        cls,
        mdata: MuData,
        mc_layer: str,
        cov_layer: str,
        labels_key: str,
        unlabeled_category: str,
        methylation_contexts: Iterable[str],
        batch_key: str | None = None,
        categorical_covariate_keys: Iterable[str] | None = None,
        modalities=None,
        **kwargs,
    ):
        """%(summary_mdata)s.

        Parameters
        ----------
        %(param_mdata)s
        mc_layer
            Layer containing methylated cytosine counts for each set of methylation features.
        cov_layer
            Layer containing total coverage counts for each set of methylation features.
        labels_key
            Obs field in `mdata` object containing cell type labels
        unlabeled_category
            Value of `mdata.obs[labels_key]` representing an unknown cell type label
        methylation_contexts
            List of modality fields in `mdata` object representing different methylation contexts.
            Each context must be equipped with a layer containing the number of methylated counts
            (specified by `mc_layer`) and total number of counts (specified by `cov_layer`) for
            each genomic region feature.
        %(param_batch_key)s
        %(param_categorical_covariate_keys)s
        %(param_modalities)s

        Examples
        --------
        METHYLANVI.setup_mudata(
            mdata,
            mc_layer="mc",
            cov_layer="cov",
            labels_key="CellType",
            unlabeled_category="Unknown",
            methylation_contexts=["mCG", "mCH"],
            categorical_covariate_keys=["Platform"],
            modalities={
                "categorical_covariate_keys": "mCG"
            },
        )

        """
        if modalities is None:
            modalities = {}
        setup_method_args = METHYLANVI._get_setup_method_args(**locals())

        if methylation_contexts is None:
            raise ValueError("Methylation contexts cannot be None.")

        modalities_ = cls._create_modalities_attr_dict(modalities, setup_method_args)

        batch_field = fields.MuDataCategoricalObsField(
            REGISTRY_KEYS.BATCH_KEY,
            batch_key,
            mod_key=modalities_.batch_key,
        )

        cat_cov_field = fields.MuDataCategoricalJointObsField(
            REGISTRY_KEYS.CAT_COVS_KEY,
            categorical_covariate_keys,
            mod_key=modalities_.categorical_covariate_keys,
        )
        cell_type_field = fields.LabelsWithUnlabeledObsField(
            REGISTRY_KEYS.LABELS_KEY, labels_key, unlabeled_category
        )

        mc_fields = []
        cov_fields = []

        for context in methylation_contexts:
            mc_fields.append(
                fields.MuDataLayerField(
                    _context_mc_key(context),
                    mc_layer,
                    mod_key=context,
                    is_count_data=True,
                    mod_required=True,
                )
            )

            cov_fields.append(
                fields.MuDataLayerField(
                    _context_cov_key(context),
                    cov_layer,
                    mod_key=context,
                    is_count_data=True,
                    mod_required=True,
                )
            )

        mudata_fields = (
            mc_fields + cov_fields + [batch_field] + [cat_cov_field] + [cell_type_field]
        )
        adata_manager = AnnDataManager(fields=mudata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(mdata, **kwargs)

        cls.register_manager(adata_manager)
