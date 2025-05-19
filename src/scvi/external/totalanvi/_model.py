from __future__ import annotations

import logging
import warnings
from copy import deepcopy
from typing import TYPE_CHECKING

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager, fields
from scvi.data._constants import _SETUP_ARGS_KEY
from scvi.model import TOTALVI
from scvi.model._utils import _init_library_size
from scvi.model.base import SemisupervisedTrainingMixin
from scvi.utils._docstrings import setup_anndata_dsp

from ._module import TOTALANVAE

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData
    from mudata import MuData

logger = logging.getLogger(__name__)


class TOTALANVI(SemisupervisedTrainingMixin, TOTALVI):
    """total Variational Inference :cite:p:`GayosoSteier21`.

    Parameters
    ----------
    adata
        AnnData object that has been registered via
        :meth:`~scvi.external.model.TOTALANVI.setup_anndata`.
    n_latent
        Dimensionality of the latent space.
    gene_dispersion
        One of the following:

        * ``'gene'`` - genes_dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - genes_dispersion can differ between different batches
        * ``'gene-label'`` - genes_dispersion can differ between different labels
    protein_dispersion
        One of the following:

        * ``'protein'`` - protein_dispersion parameter is constant per protein across cells
        * ``'protein-batch'`` - protein_dispersion can differ between different batches NOT TESTED
        * ``'protein-label'`` - protein_dispersion can differ between different labels NOT TESTED
    gene_likelihood
        One of:

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
    latent_distribution
        One of:

        * ``'normal'`` - Normal distribution
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)
    empirical_protein_background_prior
        Set the initialization of protein background prior empirically. This option fits a GMM for
        each of 100 cells per batch and averages the distributions. Note that even with this option
        set to `True`, this only initializes a parameter that is learned during inference. If
        `False`, randomly initializes. The default (`None`), sets this to `True` if greater than 10
        proteins are used.
    override_missing_proteins
        If `True`, will not treat proteins with all 0 expression in a particular batch as missing.
    linear_classifier
        If ``True``, uses a single linear layer for classification instead of a
        multi-layer perceptron.
    **model_kwargs
        Keyword args for :class:`~scvi.module.TOTALVAE`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.TOTALANVI.setup_anndata(
            adata, batch_key="batch", protein_expression_obsm_key="protein_expression"
        )
    >>> vae = scvi.model.TOTALANVI(adata)
    >>> vae.train()
    >>> adata.obsm["X_totalANVI"] = vae.get_latent_representation()
    >>> adata.obs["pred_label"] = vae.predict()

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/multimodal/totalANVI`
    2. :doc:`/tutorials/notebooks/multimodal/cite_scrna_integration_w_totalANVI`
    3. :doc:`/tutorials/notebooks/scrna/scarches_scvi_tools`
    """

    _module_cls = TOTALANVAE

    def __init__(
        self,
        adata: AnnData,
        n_latent: int = 20,
        gene_dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        protein_dispersion: Literal["protein", "protein-batch", "protein-label"] = "protein",
        gene_likelihood: Literal["zinb", "nb"] = "nb",
        latent_distribution: Literal["normal", "ln"] = "normal",
        empirical_protein_background_prior: bool | None = None,
        override_missing_proteins: bool = False,
        linear_classifier: bool = False,
        **model_kwargs,
    ):
        super().__init__(adata)
        self._set_indices_and_labels()
        self.protein_state_registry = self.adata_manager.get_state_registry(
            REGISTRY_KEYS.PROTEIN_EXP_KEY
        )
        if (
            fields.ProteinObsmField.PROTEIN_BATCH_MASK in self.protein_state_registry
            and not override_missing_proteins
        ):
            batch_mask = self.protein_state_registry.protein_batch_mask
            msg = (
                "Some proteins have all 0 counts in some batches. "
                + "These proteins will be treated as missing measurements; however, "
                + "this can occur due to experimental design/biology. "
                + "Reinitialize the model with `override_missing_proteins=True`,"
                + "to override this behavior."
            )
            warnings.warn(msg, UserWarning, stacklevel=settings.warnings_stacklevel)
            self._use_adversarial_classifier = True
        else:
            batch_mask = None
            self._use_adversarial_classifier = False

        emp_prior = (
            empirical_protein_background_prior
            if empirical_protein_background_prior is not None
            else (self.summary_stats.n_proteins > 10)
        )
        if emp_prior:
            prior_mean, prior_scale = self._get_totalvi_protein_priors(adata)
        else:
            prior_mean, prior_scale = None, None

        n_cats_per_cov = (
            self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY)[
                fields.CategoricalJointObsField.N_CATS_PER_KEY
            ]
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )
        n_labels = self.summary_stats.n_labels - 1
        n_cats_per_cov = (
            self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )

        n_batch = self.summary_stats.n_batch
        if "n_panel" in self.summary_stats:
            n_panel = self.summary_stats.n_panel
            panel_key = "panel"
        else:
            n_panel = self.summary_stats.n_batch
            panel_key = REGISTRY_KEYS.BATCH_KEY

        use_size_factor_key = REGISTRY_KEYS.SIZE_FACTOR_KEY in self.adata_manager.data_registry
        library_log_means, library_log_vars = None, None
        if not use_size_factor_key:
            library_log_means, library_log_vars = _init_library_size(self.adata_manager, n_batch)

        self.module = self._module_cls(
            n_input_genes=self.summary_stats.n_vars,
            n_input_proteins=self.summary_stats.n_proteins,
            n_batch=n_batch,
            n_labels=n_labels,
            n_latent=n_latent,
            n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
            n_cats_per_cov=n_cats_per_cov,
            gene_dispersion=gene_dispersion,
            protein_dispersion=protein_dispersion,
            gene_likelihood=gene_likelihood,
            latent_distribution=latent_distribution,
            protein_batch_mask=batch_mask,
            protein_background_prior_mean=prior_mean,
            protein_background_prior_scale=prior_scale,
            use_size_factor_key=use_size_factor_key,
            library_log_means=library_log_means,
            library_log_vars=library_log_vars,
            n_panel=n_panel,
            panel_key=panel_key,
            linear_classifier=linear_classifier,
            **model_kwargs,
        )
        self._model_summary_string = (
            f"TotalANVI Model with the following params: \nunlabeled_category: "
            f"{self.unlabeled_category_}, n_latent:{n_latent}, gene_dispersion: {gene_dispersion},"
            f" protein_dispersion: {protein_dispersion}, gene_likelihood: {gene_likelihood}, "
            f"latent_distribution: {latent_distribution}"
        )
        self.unsupervised_history_ = None
        self.semisupervised_history_ = None

        self.init_params_ = self._get_init_params(locals())
        self.was_pretrained = False
        self.n_labels = n_labels

    @classmethod
    def from_totalvi_model(
        cls,
        totalvi_model: TOTALVI,
        unlabeled_category: str,
        labels_key: str | None = None,
        adata: AnnData | None = None,
        **totalanvi_kwargs,
    ):
        """Initialize totalVI model with weights from pretrained :class:`~scvi.model.TOTALVI` model

        Parameters
        ----------
        totalvi_model
            Pretrained totalvi model
        labels_key
            key in `adata.obs` for label information. Label categories can not be different if
            labels_key was used to setup the TOTALVI model. If None, uses the `labels_key` used to
            setup the TOTALVI model. If that was None, and error is raised.
        unlabeled_category
            Value used for unlabeled cells in `labels_key` used to setup AnnData with scvi.
        adata
            AnnData object that has been registered via
            :meth:`~scvi.model.external.TOTALANVI.setup_anndata`.
        totalanvi_kwargs
            kwargs for totalANVI model
        """
        totalvi_model._check_if_trained(message="Passed in scvi model hasn't been trained yet.")

        totalanvi_kwargs = dict(totalanvi_kwargs)
        init_params = totalvi_model.init_params_
        non_kwargs = init_params["non_kwargs"]
        kwargs = init_params["kwargs"]
        kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        for k, v in {**non_kwargs, **kwargs}.items():
            if k in totalanvi_kwargs.keys():
                warnings.warn(
                    f"Ignoring param '{k}' as it was already passed in to pretrained "
                    f"SCVI model with value {v}.",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )
                del totalanvi_kwargs[k]

        if adata is None:
            adata = totalvi_model.adata
        else:
            totalvi_model._validate_anndata(adata)

        totalvi_setup_args = deepcopy(totalvi_model.adata_manager.registry[_SETUP_ARGS_KEY])
        if labels_key is None:
            raise ValueError(
                "A `labels_key` is necessary as the TOTALVI model was initialized without one."
            )
        totalvi_setup_args.update({"labels_key": labels_key})
        cls.setup_anndata(
            adata,
            unlabeled_category=unlabeled_category,
            **totalvi_setup_args,
        )
        totalanvi_model = cls(adata, **non_kwargs, **kwargs, **totalanvi_kwargs)
        totalvi_state_dict = totalvi_model.module.state_dict()
        totalanvi_model.module.load_state_dict(totalvi_state_dict, strict=False)
        logging.info(
            "Sample level parameters are not optimized during supervised model. "
            "This helps with integration"
        )
        totalanvi_model.module.background_pro_alpha.requires_grad = False
        totalanvi_model.module.background_pro_log_beta.requires_grad = False
        totalanvi_model.module.log_per_batch_efficiency.requires_grad = False
        totalanvi_model.was_pretrained = True

        return totalanvi_model

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        protein_expression_obsm_key: str,
        protein_names_uns_key: str | None = None,
        batch_key: str | None = None,
        panel_key: str | None = None,
        layer: str | None = None,
        size_factor_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        labels_key: str | None = None,
        unlabeled_category: str | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        protein_expression_obsm_key
            key in `adata.obsm` for protein expression data.
        protein_names_uns_key
            key in `adata.uns` for protein names. If None, will use the column names of
            `adata.obsm[protein_expression_obsm_key]` if it is a DataFrame, else will assign
            sequential names to proteins.
        %(param_batch_key)s
        panel_key
            key in 'adata.obs' for the various panels used to measure proteins.
        %(param_layer)s
        %(param_size_factor_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        %(param_labels_key)s
        %(param_unlabeled_category)s

        Returns
        -------
        %(returns)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        if panel_key is not None:
            batch_field = fields.CategoricalObsField("panel", panel_key)
        else:
            batch_field = fields.CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key)
        anndata_fields = [
            fields.LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            fields.LabelsWithUnlabeledObsField(
                REGISTRY_KEYS.LABELS_KEY, labels_key, unlabeled_category
            ),
            fields.CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            fields.NumericalObsField(
                REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False
            ),
            fields.CategoricalJointObsField(
                REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys
            ),
            fields.NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
            fields.ProteinObsmField(
                REGISTRY_KEYS.PROTEIN_EXP_KEY,
                protein_expression_obsm_key,
                use_batch_mask=True,
                batch_field=batch_field,
                colnames_uns_key=protein_names_uns_key,
                is_count_data=True,
            ),
        ]
        if panel_key is not None:
            anndata_fields.insert(0, fields.CategoricalObsField("panel", panel_key))

        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_mudata(
        cls,
        mdata: MuData,
        rna_layer: str | None = None,
        protein_layer: str | None = None,
        batch_key: str | None = None,
        panel_key: str | None = None,
        size_factor_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        modalities: dict[str, str] | None = None,
        **kwargs,
    ):
        """%(summary_mdata)s.

        Parameters
        ----------
        %(param_mdata)s
        rna_layer
            RNA layer key. If `None`, will use `.X` of specified modality key.
        protein_layer
            Protein layer key. If `None`, will use `.X` of specified modality key.
        %(param_batch_key)s
        panel_key
            key in 'adata.obs' for the various panels used to measure proteins.
        %(param_size_factor_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        %(param_modalities)s

        Examples
        --------
        >>> mdata = muon.read_10x_h5("pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5")
        >>> scvi.external.model.TOTALANVI.setup_mudata(
                mdata, modalities={"rna_layer": "rna": "protein_layer": "prot"}
            )
        >>> vae = scvi.external.model.TOTALANVI(mdata)
        """
        setup_method_args = cls._get_setup_method_args(**locals())

        if modalities is None:
            raise ValueError("Modalities cannot be None.")
        modalities = cls._create_modalities_attr_dict(modalities, setup_method_args)

        if panel_key is not None:
            batch_field = fields.MuDataCategoricalObsField(
                "panel",
                panel_key,
                mod_key=modalities.batch_key,
            )
        else:
            batch_field = fields.MuDataCategoricalObsField(
                REGISTRY_KEYS.BATCH_KEY,
                batch_key,
                mod_key=modalities.batch_key,
            )
        mudata_fields = [
            fields.MuDataLayerField(
                REGISTRY_KEYS.X_KEY,
                rna_layer,
                mod_key=modalities.rna_layer,
                is_count_data=True,
                mod_required=True,
            ),
            fields.MuDataCategoricalObsField(
                REGISTRY_KEYS.LABELS_KEY,
                None,
                mod_key=None,
            ),  # Default labels field for compatibility with TOTALVAE
            fields.MuDataCategoricalObsField(
                REGISTRY_KEYS.BATCH_KEY,
                batch_key,
                mod_key=modalities.batch_key,
            ),
            fields.MuDataNumericalObsField(
                REGISTRY_KEYS.SIZE_FACTOR_KEY,
                size_factor_key,
                mod_key=modalities.size_factor_key,
                required=False,
            ),
            fields.MuDataCategoricalJointObsField(
                REGISTRY_KEYS.CAT_COVS_KEY,
                categorical_covariate_keys,
                mod_key=modalities.categorical_covariate_keys,
            ),
            fields.MuDataNumericalJointObsField(
                REGISTRY_KEYS.CONT_COVS_KEY,
                continuous_covariate_keys,
                mod_key=modalities.continuous_covariate_keys,
            ),
            fields.MuDataProteinLayerField(
                REGISTRY_KEYS.PROTEIN_EXP_KEY,
                protein_layer,
                mod_key=modalities.protein_layer,
                use_batch_mask=True,
                batch_field=batch_field,
                is_count_data=True,
                mod_required=True,
            ),
        ]

        if panel_key:
            mudata_fields.insert(
                0,
                fields.MuDataCategoricalObsField(
                    "panel",
                    panel_key,
                    mod_key=modalities.batch_key,
                ),
            )

        adata_manager = AnnDataManager(fields=mudata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(mdata, **kwargs)
        cls.register_manager(adata_manager)
