"""Model class for methylVI for single cell methylation data."""

from __future__ import annotations

import inspect
import logging
import warnings
from collections import defaultdict
from collections.abc import Iterable, Sequence
from functools import partial
from typing import Literal

import numpy as np
import pandas as pd
import sparse
import torch
from anndata import AnnData
from mudata import MuData

from scvi import REGISTRY_KEYS, settings
from scvi._types import Number
from scvi.data import AnnDataManager, _constants, fields
from scvi.data._constants import _MODEL_NAME_KEY, _SETUP_ARGS_KEY
from scvi.data.fields import (
    CategoricalObsField,
    LayerField,
)
from scvi.distributions._utils import DistributionConcatenator
from scvi.external.methylvi import METHYLVAE, METHYLVI_REGISTRY_KEYS
from scvi.model._utils import _get_batch_code_from_category, parse_device_args
from scvi.model.base import (
    ArchesMixin,
    BaseModelClass,
    UnsupervisedTrainingMixin,
    VAEMixin,
)
from scvi.model.base._archesmixin import _get_loaded_data, _set_params_online_update
from scvi.model.base._utils import (
    _de_core,
    _initialize_model,
    _validate_var_names,
)
from scvi.utils import setup_anndata_dsp

from ._utils import scmc_raw_counts_properties

logger = logging.getLogger(__name__)


class MethylVIModel(VAEMixin, UnsupervisedTrainingMixin, ArchesMixin, BaseModelClass):
    """
    Model class for methylVI

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~methyl_vi.MethylVI.setup_anndata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    **model_kwargs
        Keyword args for :class:`~MethylVIModule`

    Examples
    --------
    >>> mdata = anndata.read_h5mu(path_to_mudata)
    >>> MethylVI.setup_mudata(mdata, batch_key="batch")
    >>> vae = MethylVI(adata)
    >>> vae.train()
    >>> mdata['mCG'].obsm["X_methylVI"] = vae.get_latent_representation()
    """

    def __init__(
        self,
        adata: AnnData | MuData,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        **model_kwargs,
    ):
        super().__init__(adata)

        n_batch = self.summary_stats.n_batch

        # We feed in both the number of methylated counts (mc) and the
        # total number of counts (cov) as inputs
        adata_manager = self.get_anndata_manager(adata)

        if type(adata) == AnnData:
            self.modalities = ["mCG"]
            self.num_features_per_modality = [adata.shape[1]]
            n_input = adata.layers["cov"].shape[1]
        else:
            self.modalities = [
                x.split("_")[0] for x in adata_manager.data_registry.keys() if x.endswith("cov")
            ]
            self.num_features_per_modality = [
                adata[modality].shape[1] for modality in self.modalities
            ]
            n_input = np.sum(self.num_features_per_modality)

        self.module = METHYLVAE(
            n_input=n_input,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            n_batch=n_batch,
            modalities=self.modalities,
            num_features_per_modality=self.num_features_per_modality,
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
        mc_layer: str,
        cov_layer: str,
        batch_key: str | None = None,
        labels_key: str | None = None,
        **kwargs,
    ) -> AnnData | None:
        """
        %(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_batch_key)s
        %(param_labels_key)s
        %(param_layer)s

        Returns
        -------
        %(returns)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())

        # If we don't have
        anndata_fields = [
            LayerField(METHYLVI_REGISTRY_KEYS.MC_KEY, mc_layer, is_count_data=True),
            LayerField(METHYLVI_REGISTRY_KEYS.COV_KEY, cov_layer, is_count_data=True),
        ]

        anndata_fields = anndata_fields + [
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_mudata(
        cls,
        mdata: MuData,
        mc_layer: str,
        cov_layer: str,
        batch_key: str | None = None,
        size_factor_key: str | None = None,
        methylation_modalities: dict[str, str] | None = None,
        covariate_modalities=None,
        **kwargs,
    ):
        """%(summary_mdata)s.

        Parameters
        ----------
        %(param_mdata)s
        mc_layer
            Layer containing methylated cytosine counts for each methylation modality.
        cov_layer
            Layer containing total coverage counts for each methylation modality.
        %(param_batch_key)s
        %(param_size_factor_key)s
        %(param_methylation_modalities)s
        %(param_covariate_modalities)s

        Examples
        --------
        MethylVI.setup_mudata(
            mdata,
            mc_layer="mc",
            cov_layer="cov",
            batch_key="Platform",
            methylation_modalities={
                "mCG": "mCG",
                "mCH": "mCH"
            },
            covariate_modalities={
                "batch_key": "mCG"
            },
        )

        """
        if covariate_modalities is None:
            covariate_modalities = {}
        setup_method_args = MethylVIModel._get_setup_method_args(**locals())

        if methylation_modalities is None:
            raise ValueError("Methylation modalities cannot be None.")

        covariate_modalities_ = cls._create_modalities_attr_dict(
            covariate_modalities, setup_method_args
        )

        batch_field = fields.MuDataCategoricalObsField(
            REGISTRY_KEYS.BATCH_KEY,
            batch_key,
            mod_key=covariate_modalities_.batch_key,
        )

        mc_fields = []
        cov_fields = []

        for mod in methylation_modalities:
            mc_fields.append(
                fields.MuDataLayerField(
                    f"{mod}_{METHYLVI_REGISTRY_KEYS.MC_KEY}",
                    mc_layer,
                    mod_key=methylation_modalities[mod],
                    is_count_data=True,
                    mod_required=True,
                )
            )

            cov_fields.append(
                fields.MuDataLayerField(
                    f"{mod}_{METHYLVI_REGISTRY_KEYS.COV_KEY}",
                    cov_layer,
                    mod_key=methylation_modalities[mod],
                    is_count_data=True,
                    mod_required=True,
                )
            )

        mudata_fields = mc_fields + cov_fields + [batch_field]
        adata_manager = AnnDataManager(fields=mudata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(mdata, **kwargs)
        cls.register_manager(adata_manager)

    @torch.inference_mode()
    def posterior_predictive_sample(
        self,
        mdata: MuData | None = None,
        n_samples: int = 1,
        batch_size: int | None = None,
    ) -> dict[str, sparse.GCXS]:
        r"""
        Generate observation samples from the posterior predictive distribution.

        The posterior predictive distribution is written as :math:`p(\hat{x} \mid x)`.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        n_samples
            Number of samples for each cell.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        x_new : :py:class:`torch.Tensor`
            tensor with shape (n_cells, n_genes, n_samples)
        """
        adata = self._validate_anndata(mdata)

        scdl = self._make_data_loader(adata=adata, batch_size=batch_size)

        x_new = defaultdict(list)
        for tensors in scdl:
            samples = self.module.sample(
                tensors,
                n_samples=n_samples,
            )

            for modality in self.modalities:
                x_new[modality].append(sparse.GCXS.from_numpy(samples[modality].numpy()))

        for modality in self.modalities:
            x_new[modality] = sparse.concatenate(
                x_new[modality]
            )  # Shape (n_cells, n_genes, n_samples)

        return x_new

    def _get_transform_batch_gen_kwargs(self, batch):
        if "transform_batch" in inspect.signature(self.module.generative).parameters:
            return {"transform_batch": batch}
        else:
            raise NotImplementedError("Transforming batches is not implemented for this model.")

    @torch.inference_mode()
    def get_normalized_expression(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        transform_batch: Sequence[Number | str] | None = None,
        gene_list: Sequence[str] | None = None,
        n_samples: int = 1,
        n_samples_overall: int = None,
        weights: Literal["uniform", "importance"] | None = None,
        batch_size: int | None = None,
        return_mean: bool = True,
        return_numpy: bool | None = None,
        **importance_weighting_kwargs,
    ) -> (np.ndarray | pd.DataFrame) | dict[str, np.ndarray | pd.DataFrame]:
        r"""Returns the normalized (decoded) methylation expression.

        This is denoted as :math:`\mu_n` in the methylVI paper.

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
            DataFrame includes gene names as columns. If either `n_samples=1` or
            `return_mean=True`, defaults to `False`. Otherwise, it defaults to `True`.
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

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                warnings.warn(
                    "`return_numpy` must be `True` if `n_samples > 1` and `return_mean` "
                    "is`False`, returning an `np.ndarray`.",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )
            return_numpy = True

        store_distributions = weights == "importance"
        if store_distributions and len(transform_batch) > 1:
            raise NotImplementedError(
                "Importance weights cannot be computed when expression "
                "levels are averaged across batches."
            )

        exprs = defaultdict(list)
        zs = []
        qz_store = DistributionConcatenator()
        px_store = DistributionConcatenator()
        for tensors in scdl:
            per_batch_exprs = defaultdict(list)
            for batch in transform_batch:
                generative_kwargs = self._get_transform_batch_gen_kwargs(batch)
                inference_kwargs = {"n_samples": n_samples}
                inference_outputs, generative_outputs = self.module.forward(
                    tensors=tensors,
                    inference_kwargs=inference_kwargs,
                    generative_kwargs=generative_kwargs,
                    compute_loss=False,
                )

                for modality in self.modalities:
                    exp_ = generative_outputs["px_mu"][modality]
                    exp_ = exp_[..., gene_mask]
                    per_batch_exprs[modality].append(exp_[None].cpu())
                if store_distributions:
                    qz_store.store_distribution(inference_outputs["qz"])
                    px_store.store_distribution(generative_outputs["px"])

            zs.append(inference_outputs["z"].cpu())

            for modality in self.modalities:
                per_batch_exprs[modality] = (
                    torch.cat(per_batch_exprs[modality], dim=0).mean(0).numpy()
                )
                exprs[modality].append(per_batch_exprs[modality])

        cell_axis = 1 if n_samples > 1 else 0

        for modality in self.modalities:
            exprs[modality] = np.concatenate(exprs[modality], axis=cell_axis)

        zs = torch.concat(zs, dim=cell_axis)

        if n_samples_overall is not None:
            # Converts the 3d tensor to a 2d tensor
            for modality in self.modalities:
                exprs[modality] = exprs[modality].reshape(-1, exprs[modality].shape[-1])
                n_samples_ = exprs[modality].shape[0]
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
                exprs[modality] = exprs[modality][ind_]
                return_numpy = True

        elif n_samples > 1 and return_mean:
            for modality in self.modalities:
                exprs[modality] = exprs[modality].mean(0)

        if return_numpy is None or return_numpy is False:
            exprs_dfs = {}
            if len(self.modalities) > 1:
                for modality in self.modalities:
                    exprs_dfs[modality] = pd.DataFrame(
                        exprs[modality],
                        columns=adata[modality].var_names[gene_mask],
                        index=adata[modality].obs_names[indices],
                    )
                return exprs_dfs
            else:
                modality = self.modalities[0]
                return pd.DataFrame(
                    exprs,
                    columns=adata[modality].var_names[gene_mask],
                    index=adata[modality].obs_names[indices],
                )
        else:
            if len(self.modalities) > 1:
                return exprs
            else:
                return exprs[self.modalities[0]]

    @torch.inference_mode()
    def get_specific_normalized_expression(
        self,
        adata: AnnData | None = None,
        modality: str = None,
        indices: Sequence[int] | None = None,
        transform_batch: Sequence[Number | str] | None = None,
        gene_list: Sequence[str] | None = None,
        n_samples: int = 1,
        n_samples_overall: int = None,
        weights: Literal["uniform", "importance"] | None = None,
        batch_size: int | None = None,
        return_mean: bool = True,
        return_numpy: bool | None = None,
        **importance_weighting_kwargs,
    ) -> (np.ndarray | pd.DataFrame) | dict[str, np.ndarray | pd.DataFrame]:
        r"""Returns normalized methylation expression for a specific modality (e.g. mCG, mCH).

        This is denoted as :math:`\mu_n` in the methylVI paper.

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
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`.
            DataFrame includes gene names as columns. If either `n_samples=1` or
            `return_mean=True`, defaults to `False`. Otherwise, it defaults to `True`.
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
        exprs = self.get_normalized_expression(
            adata=adata,
            indices=indices,
            transform_batch=transform_batch,
            gene_list=gene_list,
            n_samples=n_samples,
            n_samples_overall=n_samples_overall,
            weights=weights,
            batch_size=batch_size,
            return_mean=return_mean,
            return_numpy=return_numpy,
            **importance_weighting_kwargs,
        )
        return exprs[modality]

    def differential_methylation(
        self,
        modality: str,
        adata: AnnData | None = None,
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
    ) -> pd.DataFrame:
        r"""\.

        A unified method for differential methylation analysis.

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
        **kwargs
            Keyword args for :meth:`scvi.model.base.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential accessibility DataFrame with the following columns:
        prob_da
            the probability of the region being differentially accessible
        is_da_fdr
            whether the region passes a multiple hypothesis correction procedure
            with the target_fdr threshold
        bayes_factor
            Bayes Factor indicating the level of significance of the analysis
        effect_size
            the effect size, computed as (accessibility in population 2) -
            (accessibility in population 1)
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
        col_names = adata[modality].var_names
        model_fn = partial(
            self.get_specific_normalized_expression,
            batch_size=batch_size,
            modality=modality,
        )

        # TODO check if change_fn in kwargs and raise error if so
        def change_fn(a, b):
            return a - b

        if two_sided:

            def m1_domain_fn(samples):
                return np.abs(samples) >= delta

        else:

            def m1_domain_fn(samples):
                return samples >= delta

        result = _de_core(
            adata_manager=self.get_anndata_manager(adata, required=True),
            model_fn=model_fn,
            representation_fn=None,
            groupby=groupby,
            group1=group1,
            group2=group2,
            idx1=idx1,
            idx2=idx2,
            all_stats=all_stats,
            all_stats_fn=partial(scmc_raw_counts_properties, modality=modality),
            col_names=col_names,
            mode=mode,
            batchid1=batchid1,
            batchid2=batchid2,
            delta=delta,
            batch_correction=batch_correction,
            fdr=fdr_target,
            silent=silent,
            **kwargs,
        )

        return result

    @classmethod
    def load_query_mdata(
        cls,
        mdata: MuData,
        reference_model: str | BaseModelClass,
        inplace_subset_query_vars: bool = False,
        accelerator: str = "auto",
        device: int | str = "auto",
        unfrozen: bool = False,
        freeze_dropout: bool = False,
        freeze_expression: bool = True,
        freeze_decoder_first_layer: bool = True,
        freeze_batchnorm_encoder: bool = True,
        freeze_batchnorm_decoder: bool = False,
        freeze_classifier: bool = True,
    ):
        """Online update of a reference model with scArches algorithm :cite:p:`Lotfollahi21`.

        Parameters
        ----------
        adata
            AnnData organized in the same way as data used to train model.
            It is not necessary to run setup_anndata,
            as AnnData is validated against the ``registry``.
        reference_model
            Either an already instantiated model of the same class, or a path to
            saved outputs for reference model.
        inplace_subset_query_vars
            Whether to subset and rearrange query vars inplace based on vars used to
            train reference model.
        %(param_accelerator)s
        %(param_device)s
        unfrozen
            Override all other freeze options for a fully unfrozen model
        freeze_dropout
            Whether to freeze dropout during training
        freeze_expression
            Freeze neurons corersponding to expression in first layer
        freeze_decoder_first_layer
            Freeze neurons corersponding to first layer in decoder
        freeze_batchnorm_encoder
            Whether to freeze batchnorm weight and bias during training for encoder
        freeze_batchnorm_decoder
            Whether to freeze batchnorm weight and bias during training for decoder
        freeze_classifier
            Whether to freeze classifier completely. Only applies to `SCANVI`.
        """
        _, _, device = parse_device_args(
            accelerator=accelerator,
            devices=device,
            return_device="torch",
            validate_single_device=True,
        )

        attr_dict, var_names, load_state_dict = _get_loaded_data(reference_model, device=device)

        if inplace_subset_query_vars:
            logger.debug("Subsetting query vars to reference vars.")
            mdata._inplace_subset_var(var_names)
        _validate_var_names(mdata, var_names)

        registry = attr_dict.pop("registry_")
        if _MODEL_NAME_KEY in registry and registry[_MODEL_NAME_KEY] != cls.__name__:
            raise ValueError("It appears you are loading a model from a different class.")

        if _SETUP_ARGS_KEY not in registry:
            raise ValueError(
                "Saved model does not contain original setup inputs. "
                "Cannot load the original setup."
            )

        cls.setup_mudata(
            mdata,
            source_registry=registry,
            extend_categories=True,
            allow_missing_labels=True,
            **registry[_SETUP_ARGS_KEY],
        )

        model = _initialize_model(cls, mdata, attr_dict)
        adata_manager = model.get_anndata_manager(mdata, required=True)

        if REGISTRY_KEYS.CAT_COVS_KEY in adata_manager.data_registry:
            raise NotImplementedError(
                "scArches currently does not support models with extra categorical covariates."
            )

        version_split = adata_manager.registry[_constants._SCVI_VERSION_KEY].split(".")
        if int(version_split[1]) < 8 and int(version_split[0]) == 0:
            warnings.warn(
                "Query integration should be performed using models trained with "
                "version >= 0.8",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )

        model.to_device(device)

        # model tweaking
        new_state_dict = model.module.state_dict()
        for key, load_ten in load_state_dict.items():
            new_ten = new_state_dict[key]
            if new_ten.size() == load_ten.size():
                continue
            # new categoricals changed size
            else:
                dim_diff = new_ten.size()[-1] - load_ten.size()[-1]
                fixed_ten = torch.cat([load_ten, new_ten[..., -dim_diff:]], dim=-1)
                load_state_dict[key] = fixed_ten

        model.module.load_state_dict(load_state_dict)
        model.module.eval()

        _set_params_online_update(
            model.module,
            unfrozen=unfrozen,
            freeze_decoder_first_layer=freeze_decoder_first_layer,
            freeze_batchnorm_encoder=freeze_batchnorm_encoder,
            freeze_batchnorm_decoder=freeze_batchnorm_decoder,
            freeze_dropout=freeze_dropout,
            freeze_expression=freeze_expression,
            freeze_classifier=freeze_classifier,
        )
        model.is_trained_ = False

        return model
