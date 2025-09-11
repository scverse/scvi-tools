from __future__ import annotations

import logging
import warnings
from functools import partial
from typing import TYPE_CHECKING

import anndata
import joblib
import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from rich import print

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data._constants import _DATA_REGISTRY_KEY, _FIELD_REGISTRIES_KEY, _STATE_REGISTRY_KEY
from scvi.data._utils import _get_adata_minify_type
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
    ObsmField,
)
from scvi.model._utils import _init_library_size, scrna_raw_counts_properties
from scvi.model.base import (
    ArchesMixin,
    BaseMinifiedModeModelClass,
    EmbeddingMixin,
    RNASeqMixin,
    UnsupervisedTrainingMixin,
    VAEMixin,
)
from scvi.model.base._archesmixin import _get_loaded_data
from scvi.model.base._de_core import _de_core
from scvi.utils import de_dsp, setup_anndata_dsp, unsupported_if_adata_minified

from ._constants import SCVIVA_REGISTRY_KEYS
from ._module import nicheVAE
from .differential_expression import _niche_de_core

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence
    from typing import Literal

    from torch import Tensor

    from scvi.model.base import (
        BaseModelClass,
    )

    from .differential_expression import DifferentialExpressionResults

from scipy.sparse import csr_matrix

_SCVI_LATENT_QZM = "_scvi_latent_qzm"
_SCVI_LATENT_QZV = "_scvi_latent_qzv"
_SCVI_OBSERVED_LIB_SIZE = "_scvi_observed_lib_size"

logger = logging.getLogger(__name__)


class SCVIVA(
    EmbeddingMixin,
    RNASeqMixin,
    VAEMixin,
    ArchesMixin,
    UnsupervisedTrainingMixin,
    BaseMinifiedModeModelClass,
):
    """scVIVA: variational auto-encoder with niche decoders for ST:cite:p:`Levy25`.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.external.SCVIVA.setup_anndata`. If
        ``None``, then the underlying module will not be initialized until training, and a
        :class:`~lightning.pytorch.core.LightningDataModule` must be passed in during training.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    dropout_rate
        Dropout rate for neural networks.
    dispersion
        One of the following:

        * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - dispersion can differ between different batches
        * ``'gene-label'`` - dispersion can differ between different labels
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
    gene_likelihood
        One of:

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
        * ``'poisson'`` - Poisson distribution
    latent_distribution
        One of:

        * ``'normal'`` - Normal distribution
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)
    **kwargs
        Additional keyword arguments for :class:`~scvi.module.VAE`.

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.external.SCVIVA.preprocessing_anndata(
        adata,
        k_nn = 20,
        sample_key = 'slide_ID',
        labels_key = "cell_type",
        cell_coordinates_key = "spatial",
        expression_embedding_key = "X_scVI",
        **kwargs
    )
    >>> scvi.external.SCVIVA.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.external.SCVIVA(adata)
    >>> vae.train()
    >>> adata.obsm["X_scVIVA"] = vae.get_latent_representation()
    >>> adata.obsm["X_normalized_scVIVA"] = vae.get_normalized_expression()

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/quick_start/api_overview`
    2. :doc:`/tutorials/notebooks/spatial/scVIVA_tutorial`


    See Also
    --------
    :class:`~scvi.external.scviva._module.nicheVAE`
    """

    _module_cls = nicheVAE

    def __init__(
        self,
        adata: AnnData | None = None,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "poisson",
        latent_distribution: Literal["normal", "ln"] = "normal",
        **kwargs,
    ):
        super().__init__(adata)

        self.n_labels = self.summary_stats.n_labels

        self._module_kwargs = {
            "n_hidden": n_hidden,
            "n_latent": n_latent,
            "n_layers": n_layers,
            "dropout_rate": dropout_rate,
            "dispersion": dispersion,
            "gene_likelihood": gene_likelihood,
            "latent_distribution": latent_distribution,
            **kwargs,
        }
        self._model_summary_string = (
            "scVIVA model with the following parameters: \n"
            f"n_hidden: {n_hidden}, n_latent: {n_latent}, n_layers: {n_layers}, "
            f"dropout_rate: {dropout_rate}, dispersion: {dispersion}, "
            f"gene_likelihood: {gene_likelihood}, latent_distribution: {latent_distribution}."
        )

        if self._module_init_on_train:
            self.module = None
            warnings.warn(
                "Model was initialized without `adata`. The module will be initialized when "
                "calling `train`. This behavior is experimental and may change in the future.",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )
        else:
            n_cats_per_cov = (
                self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
                if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
                else None
            )
            n_batch = self.summary_stats.n_batch
            use_size_factor_key = REGISTRY_KEYS.SIZE_FACTOR_KEY in self.adata_manager.data_registry
            library_log_means, library_log_vars = None, None
            if not use_size_factor_key and self.minified_data_type is None:
                library_log_means, library_log_vars = _init_library_size(
                    self.adata_manager, n_batch
                )
            self.module = self._module_cls(
                n_input=self.summary_stats.n_vars,
                n_output_niche=self.summary_stats.n_latent_mean,
                n_batch=n_batch,
                n_labels=self.summary_stats.n_labels,
                n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
                n_cats_per_cov=n_cats_per_cov,
                n_hidden=n_hidden,
                n_latent=n_latent,
                n_layers=n_layers,
                dropout_rate=dropout_rate,
                dispersion=dispersion,
                gene_likelihood=gene_likelihood,
                latent_distribution=latent_distribution,
                use_size_factor_key=use_size_factor_key,
                library_log_means=library_log_means,
                library_log_vars=library_log_vars,
                **kwargs,
            )
            self.module.minified_data_type = self.minified_data_type

        self.init_params_ = self._get_init_params(locals())

    @staticmethod
    def preprocessing_anndata(
        adata: AnnData,
        k_nn: int = 20,
        sample_key: str | None = None,
        labels_key: str = "cell_type",
        cell_coordinates_key: str = "spatial",
        expression_embedding_key: str = "X_scVI",
        expression_embedding_niche_key: str = "niche_activation",
        niche_composition_key: str = "niche_composition",
        niche_indexes_key: str = "niche_indexes",
        niche_distances_key: str = "niche_distances",
        log1p: bool = False,
    ) -> None:
        """Preprocess an AnnData object for scVIVA analysis.

         This function prepares the input AnnData object by computing niche indexes,
         neighborhood composition, and average latent space embeddings per cell type.

        Parameters
        ----------
        adata : AnnData
             The annotated data matrix of shape `n_obs` x `n_vars`. Rows correspond to cells
             and columns to genes.
        k_nn : int, optional
             Number of nearest neighbors for niche computation. Default is 20.
        sample_key : str or None, optional
             Key in `adata.obs` for sample identifiers. Default is None.
        labels_key : str, optional
             Key in `adata.obs` for cell type labels. Default is "cell_type".
        cell_coordinates_key : str, optional
             Key in `adata.obsm` for spatial coordinates. Default is "spatial".
        expression_embedding_key : str, optional
             Key in `adata.obsm` for latent space embeddings. Default is "X_scVI".
        expression_embedding_niche_key : str, optional
             Key in `adata.obsm` where average latent embeddings per cell type are stored.
             Default is "niche_activation".
        niche_composition_key : str, optional
             Key in `adata.obsm` where neighborhood composition is stored.
             Default is "niche_composition".
        niche_indexes_key : str, optional
             Key in `adata.obsm` where niche indexes are stored. Default is "niche_indexes".
        niche_distances_key : str, optional
             Key in `adata.obsm` where neighbor distances are stored. Default is "niche_distances".
        log1p : bool, optional
             Whether to apply log1p to latent space embeddings. Default is False.

        Returns
        -------
        None
             The function modifies the input AnnData object in place.
        """
        get_niche_indexes(
            adata=adata,
            sample_key=sample_key,
            niche_indexes_key=niche_indexes_key,
            niche_distances_key=niche_distances_key,
            cell_coordinates_key=cell_coordinates_key,
            k_nn=k_nn,
        )

        get_neighborhood_composition(
            adata=adata,
            cell_type_column=labels_key,
            indices_key=niche_indexes_key,
            niche_composition_key=niche_composition_key,
        )

        get_average_latent_per_celltype(
            adata=adata,
            labels_key=labels_key,
            niche_indexes_key=niche_indexes_key,
            latent_mean_key=expression_embedding_key,
            latent_mean_ct_key=expression_embedding_niche_key,
            log1p=log1p,
        )

        return None

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        batch_key: str | None = None,
        size_factor_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        sample_key: str | None = None,
        labels_key: str = "cell_type",
        cell_coordinates_key: str = "spatial",
        expression_embedding_key: str = "X_scVI",
        expression_embedding_niche_key: str = "niche_activation",
        niche_composition_key: str = "niche_composition",
        niche_indexes_key: str = "niche_indexes",
        niche_distances_key: str = "niche_distances",
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
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
            NumericalObsField(REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
            CategoricalObsField(SCVIVA_REGISTRY_KEYS.SAMPLE_KEY, sample_key),
            ObsmField(SCVIVA_REGISTRY_KEYS.NICHE_COMPOSITION_KEY, niche_composition_key),
            ObsmField(SCVIVA_REGISTRY_KEYS.CELL_COORDINATES_KEY, cell_coordinates_key),
            ObsmField(SCVIVA_REGISTRY_KEYS.NICHE_INDEXES_KEY, niche_indexes_key),
            ObsmField(SCVIVA_REGISTRY_KEYS.NICHE_DISTANCES_KEY, niche_distances_key),
            ObsmField(SCVIVA_REGISTRY_KEYS.Z1_MEAN_KEY, expression_embedding_key),
            ObsmField(SCVIVA_REGISTRY_KEYS.Z1_MEAN_CT_KEY, expression_embedding_niche_key),
        ]
        # register new fields if the adata is minified
        adata_minify_type = _get_adata_minify_type(adata)
        if adata_minify_type is not None:
            anndata_fields += cls._get_fields_for_adata_minification(adata_minify_type)
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    @torch.inference_mode()
    def predict_neighborhood(
        self,
        adata: AnnData | None = None,
        indices: np.ndarray | None = None,
        batch_size: int | None = 1024,
    ) -> np.ndarray:
        """
        Predict the cell type composition of each cell niche in the dataset.

        Parameters
        ----------
        adata
            AnnData object. If ``None``, the model's ``adata`` will be used.
        indices
            Indices of cells to use. If ``None``, all cells will be used.
        batch_size
            Minibatch size to use during inference.

        Returns
        -------
        ct_prediction
            Predicted cell type composition of each cell niche in the dataset.
            It is computed as the expectation of the Dirichlet distribution.
        """
        self._check_if_trained(warn=False)

        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        ct_prediction = []
        for tensors in scdl:
            inference_inputs = self.module._get_inference_input(tensors)
            outputs = self.module.inference(**inference_inputs)

            batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
            decoder_input = outputs["qz"].loc

            # put batch_index in the same device as decoder_input
            batch_index = batch_index.to(decoder_input.device)

            predicted_ct_prob = self.module.composition_decoder(
                decoder_input,
                batch_index,
            )  # no batch correction here

            ct_prediction.append(
                (
                    predicted_ct_prob.concentration
                    / predicted_ct_prob.concentration.sum(dim=1).unsqueeze(1)
                )
                .detach()
                .cpu()
            )

        return torch.cat(ct_prediction).numpy()

    @torch.inference_mode()
    def predict_niche_activation(
        self,
        adata: AnnData | None = None,
        indices: np.ndarray | None = None,
        batch_size: int | None = 1024,
    ) -> np.ndarray:
        """
        Predict the activation of each cell niche in the dataset.

        Parameters
        ----------
        adata
            AnnData object. If ``None``, the model's ``adata`` will be used.
        indices
            Indices of cells to use. If ``None``, all cells will be used.
        batch_size
            Minibatch size to use during inference.

        Returns
        -------
        niche_activation
            Predicted activation of each cell niche in the dataset.
        """
        self._check_if_trained(warn=False)

        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        ct_prediction = []
        for tensors in scdl:
            inference_inputs = self.module._get_inference_input(tensors)
            outputs = self.module.inference(**inference_inputs)

            batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
            decoder_input = outputs["qz"].loc

            # put batch_index in the same device as decoder_input
            batch_index = batch_index.to(decoder_input.device)

            p_m, p_v = self.module.niche_decoder(
                decoder_input,
                batch_index,
            )  # no batch correction here

            ct_prediction.append(p_m.detach().cpu())

        return torch.cat(ct_prediction).numpy()

    @de_dsp.dedent
    def differential_expression(
        self,
        adata: AnnData | None = None,
        groupby: str | None = None,
        group1: list[str] | None = None,
        group2: str | None = None,
        idx1: list[int] | list[bool] | str | None = None,
        idx2: list[int] | list[bool] | str | None = None,
        mode: Literal["vanilla", "change"] = "change",
        delta: float | list[float] = 0.15,
        batch_size: int | None = None,
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: list[str] | None = None,
        batchid2: list[str] | None = None,
        fdr_target: float | list[float] = 0.05,
        silent: bool = False,
        weights: Literal["uniform", "importance"] | None = "uniform",
        filter_outlier_cells: bool = False,
        importance_weighting_kwargs: dict | None = None,
        ###### scVIVA specific ######
        niche_mode: bool = True,
        radius: int | None = None,
        k_nn: int | None = None,
        n_restarts_optimizer_gpc: int = 10,
        path_to_save: str | None = None,
        **kwargs,
    ) -> DifferentialExpressionResults:
        r"""A unified method for differential expression analysis.

        Implements ``'vanilla'`` DE :cite:p:`Lopez18` and ``'change'`` mode DE :cite:p:`Boyeau19`.
        Adds a neighborhood component to the DE analysis :cite:p:`Levy25`.

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
            :meth:`~scvi.model.base.RNASeqMixin.get_importance_weights`.
        niche_mode
            Whether to use scVIVA DE or SCVI DE.
        radius
            Radius for scVIVA DE.
        k_nn
            Number of nearest neighbors for scVIVA DE.
        n_restarts_optimizer_gpc
            Number of restarts for the Gaussian Process Classifier optimization.
        path_to_save
            Path to save the results to, as a pickle file.
        **kwargs
            Keyword args for :meth:`scvi.model.base.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential expression dataclass.
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

        if niche_mode:
            result = _niche_de_core(
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
                radius=radius,
                k_nn=k_nn,
                n_restarts_optimizer_gpc=n_restarts_optimizer_gpc,
                **kwargs,
            )

        else:
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

        if path_to_save is not None:
            joblib.dump(result, path_to_save)

        return result

    @torch.inference_mode()
    @unsupported_if_adata_minified
    def get_composition_error(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
        dataloader: Iterator[dict[str, Tensor | None]] = None,
        return_mean: bool = True,
        **kwargs,
    ) -> dict[str, float]:
        r"""Compute the composition prediction error on the data.

        The error is the negative log likelihood of the data (alpha) given the latent
        variables. This is typically written as
        :math:`p(alpha \mid z)`, the likelihood term given one posterior sample.

        Parameters
        ----------
        adata
            :class:`~anndata.AnnData` object with :attr:`~anndata.AnnData.var_names` in the same
            order as the ones used to train the model. If ``None`` and ``dataloader`` is also
            ``None``, it defaults to the object used to initialize the model.
        indices
            Indices of observations in ``adata`` to use. If ``None``, defaults to all observations.
            Ignored if ``dataloader`` is not ``None``
        batch_size
            Minibatch size for the forward pass. If ``None``, defaults to
            ``scvi.settings.batch_size``. Ignored if ``dataloader`` is not ``None``
        dataloader
            An iterator over minibatches of data on which to compute the metric. The minibatches
            should be formatted as a dictionary of :class:`~torch.Tensor` with keys as expected by
            the model. If ``None``, a dataloader is created from ``adata``.
        return_mean
            Whether to return the mean reconstruction loss or the reconstruction loss
            for each observation.
        **kwargs
            Additional keyword arguments to pass into the forward method of the module.

        Returns
        -------
        The composition prediction error on the data.

        Notes
        -----
        This is not the negative reconstruction error, so higher is better.
        """
        from ._log_likelihood import compute_composition_error

        if adata is not None and dataloader is not None:
            raise ValueError("Only one of `adata` or `dataloader` can be provided.")

        if dataloader is None:
            adata = self._validate_anndata(adata)
            dataloader = self._make_data_loader(
                adata=adata, indices=indices, batch_size=batch_size
            )

        return compute_composition_error(
            self.module, dataloader, return_mean=return_mean, **kwargs
        )

    @torch.inference_mode()
    @unsupported_if_adata_minified
    def get_niche_error(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
        dataloader: Iterator[dict[str, Tensor | None]] = None,
        return_mean: bool = True,
        **kwargs,
    ) -> dict[str, float]:
        r"""Compute the niche state prediction error on the data.

        The  error is the negative log likelihood of the data (eta) given the latent
        variables. This is typically written as
        :math:`p(eta \mid z)`, the likelihood term given one posterior sample.

        Parameters
        ----------
        adata
            :class:`~anndata.AnnData` object with :attr:`~anndata.AnnData.var_names` in the same
            order as the ones used to train the model. If ``None`` and ``dataloader`` is also
            ``None``, it defaults to the object used to initialize the model.
        indices
            Indices of observations in ``adata`` to use. If ``None``, defaults to all observations.
            Ignored if ``dataloader`` is not ``None``
        batch_size
            Minibatch size for the forward pass. If ``None``, defaults to
            ``scvi.settings.batch_size``. Ignored if ``dataloader`` is not ``None``
        dataloader
            An iterator over minibatches of data on which to compute the metric. The minibatches
            should be formatted as a dictionary of :class:`~torch.Tensor` with keys as expected by
            the model. If ``None``, a dataloader is created from ``adata``.
        return_mean
            Whether to return the mean reconstruction loss or the reconstruction loss
            for each observation.
        **kwargs
            Additional keyword arguments to pass into the forward method of the module.

        Returns
        -------
        The niche state prediction error of the data.

        Notes
        -----
        This is not the negative reconstruction error, so higher is better.
        """
        from ._log_likelihood import compute_niche_error

        if adata is not None and dataloader is not None:
            raise ValueError("Only one of `adata` or `dataloader` can be provided.")

        if dataloader is None:
            adata = self._validate_anndata(adata)
            dataloader = self._make_data_loader(
                adata=adata, indices=indices, batch_size=batch_size
            )

        return compute_niche_error(self.module, dataloader, return_mean=return_mean, **kwargs)

    def preprocessing_query_anndata(
        self,
        adata: AnnData,
        reference_model: str | BaseModelClass,
        return_reference_var_names: bool = False,
        inplace: bool = True,
        k_nn: int = 20,
        sample_key: str | None = None,
        labels_key: str = "cell_type",
        cell_coordinates_key: str = "spatial",
        expression_embedding_key: str = "X_scVI",
        expression_embedding_niche_key: str = "niche_activation",
        niche_composition_key: str = "niche_composition",
        niche_indexes_key: str = "niche_indexes",
        niche_distances_key: str = "niche_distances",
        log1p: bool = False,
    ) -> AnnData | pd.Index | None:
        """Prepare data for query integration.

        Merges SCVIVA.preprocessing_anndata and ArchesMixin.prepare_query_anndata.

        This function will return a new AnnData object with padded zeros
        for missing features (genes, alpha and eta), as well as correctly sorted features.

        Parameters
        ----------
        adata
            AnnData organized in the same way as data used to train model.
            It is not necessary to run setup_anndata,
            as AnnData is validated against the ``registry``.
        reference_model
            Either an already instantiated model of the same class, or a path to
            saved outputs for reference model.
        return_reference_var_names
            Only load and return reference var names if True.
        inplace
            Whether to subset and rearrange query vars inplace or return new AnnData.
        k_nn : int, optional
             Number of nearest neighbors for niche computation. Default is 20.
        sample_key : str or None, optional
             Key in `adata.obs` for sample identifiers. Default is None.
        labels_key : str, optional
             Key in `adata.obs` for cell type labels. Default is "cell_type".
        cell_coordinates_key : str, optional
             Key in `adata.obsm` for spatial coordinates. Default is "spatial".
        expression_embedding_key : str, optional
             Key in `adata.obsm` for latent space embeddings. Default is "X_scVI".
        expression_embedding_niche_key : str, optional
             Key in `adata.obsm` where average latent embeddings per cell type are stored.
             Default is "niche_activation".
        niche_composition_key : str, optional
             Key in `adata.obsm` where neighborhood composition is stored.
             Default is "niche_composition".
        niche_indexes_key : str, optional
             Key in `adata.obsm` where niche indexes are stored. Default is "niche_indexes".
        niche_distances_key : str, optional
             Key in `adata.obsm` where neighbor distances are stored. Default is "niche_distances".
        log1p : bool, optional
             Whether to apply log1p to latent space embeddings. Default is False.


        Returns
        -------
        Query adata ready to use in `load_query_data` unless `return_reference_var_names`
        in which case a pd.Index of reference var names is returned.
        """
        SCVIVA.preprocessing_anndata(
            adata=adata,
            k_nn=k_nn,
            sample_key=sample_key,
            labels_key=labels_key,
            cell_coordinates_key=cell_coordinates_key,
            expression_embedding_key=expression_embedding_key,
            expression_embedding_niche_key=expression_embedding_niche_key,
            niche_composition_key=niche_composition_key,
            niche_indexes_key=niche_indexes_key,
            niche_distances_key=niche_distances_key,
            log1p=log1p,
        )

        _, var_names, _, _ = _get_loaded_data(reference_model, device="cpu")
        var_names = pd.Index(var_names)

        if return_reference_var_names:
            return var_names

        reference_niche_composition_key = reference_model.registry[_FIELD_REGISTRIES_KEY][
            SCVIVA_REGISTRY_KEYS.NICHE_COMPOSITION_KEY
        ][_DATA_REGISTRY_KEY]["attr_key"]
        assert reference_niche_composition_key == niche_composition_key, (
            f"niche_composition_key in query ({niche_composition_key}) must match that "
            f"of reference ({reference_niche_composition_key})"
        )
        reference_expression_embedding_niche_key = reference_model.registry[_FIELD_REGISTRIES_KEY][
            SCVIVA_REGISTRY_KEYS.Z1_MEAN_CT_KEY
        ][_DATA_REGISTRY_KEY]["attr_key"]
        assert reference_expression_embedding_niche_key == expression_embedding_niche_key, (
            f"expression_embedding_niche_key in query ({expression_embedding_niche_key}) must "
            f"match that of reference ({reference_expression_embedding_niche_key})"
        )

        reference_label_names = reference_model.registry[_FIELD_REGISTRIES_KEY][
            SCVIVA_REGISTRY_KEYS.NICHE_COMPOSITION_KEY
        ][_STATE_REGISTRY_KEY]["column_names"]
        reference_label_names = pd.Index(reference_label_names)

        query_label_names = pd.Index(adata.obsm[niche_composition_key].columns)

        return _pad_and_sort_query_anndata(
            adata,
            var_names,
            reference_label_names,
            query_label_names,
            niche_composition_key,
            expression_embedding_niche_key,
            inplace,
        )


def get_niche_indexes(
    adata: AnnData,
    sample_key: str,
    niche_indexes_key: str,
    cell_coordinates_key: str,
    k_nn: int,
    niche_distances_key: str,
) -> None:
    """Get the k nearest neighbors of each cell in the dataset, grouped per sample.

    The indexes of the neighbors are stored in adata.obsm[niche_indexes_key] and the distances to
    the neighbors are stored in adata.obsm[niche_distances_key].

    Parameters
    ----------
    adata
        Anndata object
    sample_key
        Key in adata.obs that contains the sample of each cell (i.e. the donor slice)
    niche_indexes_key
        Key in adata.obsm where the indexes of the neighbors will be stored
    cell_coordinates_key
        Key in adata.obsm that contains the spatial coordinates of each cell
    k_nn
        Number of nearest neighbors to compute
    niche_distances_key
        Key in adata.obsm where the distances to the neighbors will be stored

    Returns
    -------
    None

    """
    from sklearn.neighbors import NearestNeighbors

    adata.obsm[niche_indexes_key] = np.zeros(
        (adata.n_obs, k_nn)
    )  # for each cell, store the indexes of its k_nn neighbors
    adata.obsm[niche_distances_key] = np.zeros(
        (adata.n_obs, k_nn)
    )  # for each cell, store the distances to its k_nn neighbors

    adata.obs["index"] = np.arange(adata.shape[0])
    # build a dictionnary giving the index of each 'donor_slice' observation:
    donor_slice_index = {}
    for sample in adata.obs[sample_key].unique():
        donor_slice_index[sample] = adata.obs[adata.obs[sample_key] == sample]["index"].values

    for sample in adata.obs[sample_key].unique():
        sample_coord = adata.obsm[cell_coordinates_key][adata.obs[sample_key] == sample]

        # Create a NearestNeighbors object
        knn = NearestNeighbors(n_neighbors=k_nn + 1)

        # Fit the kNN model to the points
        knn.fit(sample_coord)

        # Find the indices of the kNN for each point
        distances, indices = knn.kneighbors(sample_coord)

        # Store the indices in the adata object
        sample_global_index = donor_slice_index[sample][indices].astype(int)

        adata.obsm[niche_indexes_key][adata.obs[sample_key] == sample] = sample_global_index[:, 1:]

        adata.obsm[niche_indexes_key] = adata.obsm[niche_indexes_key].astype(int)

        adata.obsm[niche_distances_key][adata.obs[sample_key] == sample] = distances[:, 1:]

    print(
        f"[bold cyan]Saved {niche_indexes_key} and {niche_distances_key} in adata.obsm[/bold cyan]"
    )

    return None


def get_neighborhood_composition(
    adata: AnnData,
    cell_type_column: str,
    indices_key: str = "niche_indexes",
    niche_composition_key: str = "niche_composition",
) -> None:
    """Get the composition of each neighborhood in the dataset (alpha)."""
    n_cell_types = len(adata.obs[cell_type_column].unique())  # number of cell types
    adata.obsm[niche_composition_key] = np.zeros(
        (adata.n_obs, n_cell_types)
    )  # for each cell, store the composition of its neighborhood as a convex vector of cell type
    # proportions

    indices = adata.obsm[indices_key].astype(int)

    cell_types = adata.obs[cell_type_column].unique().tolist()
    cell_type_to_int = {cell_types[i]: i for i in range(len(cell_types))}
    adata.uns["cell_type_to_int"] = cell_type_to_int

    # Transform the query vector into an integer-valued vector
    integer_vector = np.vectorize(cell_type_to_int.get)(adata.obs[cell_type_column])

    n_cells = adata.n_obs
    # For each cell, get the cell types of its neighbors
    cell_types_in_the_neighborhood = [integer_vector[indices[cell, :]] for cell in range(n_cells)]

    # Compute the composition of each neighborhood
    composition = np.array(
        [
            np.bincount(
                cell_types_in_the_neighborhood[cell],
                minlength=len(cell_type_to_int),
            )
            for cell in range(n_cells)
        ]
    )

    # Normalize the composition of each neighborhood
    composition = composition / indices.shape[1]
    composition = np.array(composition)

    neighborhood_composition_df = pd.DataFrame(
        data=composition,
        columns=cell_types,
        index=adata.obs_names,
    )

    adata.obsm[niche_composition_key] = neighborhood_composition_df

    print(f"[bold green]Saved {niche_composition_key} in adata.obsm[/bold green]")

    return None


def get_average_latent_per_celltype(
    adata: AnnData,
    labels_key: str,
    niche_indexes_key: str,
    latent_mean_key: str,
    latent_mean_ct_key: str = "niche_activation",
    log1p: bool = False,
) -> None:
    """Get the average embedding per cell type in the dataset.

    For this one needs to provide the cell type of each cell in the
    dataset and an embedding for each cell, computed for instance with PCA, or scVI.

    Parameters
    ----------
    adata
        Anndata object
    labels_key
        Key in adata.obs that contains the cell type of each cell
    niche_indexes_key
        Key in adata.obsm that contains the indexes of the neighbors of each cell
    latent_mean_key
        Key in adata.obsm that contains the expression embedding of each cell
    latent_mean_ct_key
        Key in adata.obsm where the average embedding per cell type will be stored
    log1p
        Whether the latent space is log-transformed

    Returns
    -------
    None

    """
    n_cells = adata.n_obs
    n_cell_types = len(adata.obs[labels_key].unique())
    n_latent_z1 = adata.obsm[latent_mean_key].shape[1]
    niche_indexes = adata.obsm[niche_indexes_key]

    if log1p:
        z1_mean_niches = np.log1p(adata.obsm[latent_mean_key][niche_indexes])
    else:
        z1_mean_niches = adata.obsm[latent_mean_key][niche_indexes]

    cell_type_to_int = adata.uns["cell_type_to_int"]
    integer_vector = np.vectorize(cell_type_to_int.get)(adata.obs[labels_key])

    # For each cell, get the cell types of its neighbors (as integers)
    cell_types_in_the_neighborhood = np.vstack(
        [integer_vector[niche_indexes[cell, :]] for cell in range(n_cells)]
    )

    dict_of_cell_type_indices = {}

    for cell_type, cell_type_idx in cell_type_to_int.items():
        ct_row_indices, ct_col_indices = np.where(cell_types_in_the_neighborhood == cell_type_idx)

        # dict of cells:local index of the cells of cell_type in the neighborhood.
        result_dict = {}
        for row_idx, col_idx in zip(ct_row_indices, ct_col_indices, strict=False):
            result_dict.setdefault(row_idx, []).append(col_idx)

        dict_of_cell_type_indices[cell_type] = result_dict

    latent_mean_ct_prior = np.zeros((n_cell_types, n_latent_z1))

    z1_mean_niches_ct = np.stack(
        [latent_mean_ct_prior] * n_cells, axis=0
    )  # batch times n_cell_types times n_latent. Initialize your prior with a non-spatial average.

    # outer loop over cell types
    for cell_type, cell_type_idx in cell_type_to_int.items():
        ct_dict = dict_of_cell_type_indices[cell_type]
        # inner loop over every cell that has this cell type in its neighborhood.
        for cell_idx, neighbor_idxs in ct_dict.items():
            z1_mean_niches_ct[cell_idx, cell_type_idx, :] = np.mean(
                z1_mean_niches[cell_idx, neighbor_idxs, :], axis=0
            )

    adata.obsm[latent_mean_ct_key] = z1_mean_niches_ct

    print(f"[bold green]Saved {latent_mean_ct_key} in adata.obsm[/bold green]")

    return None


def _pad_and_sort_query_anndata(
    adata: AnnData,
    reference_var_names: pd.Index,
    reference_label_names: pd.Index,
    query_label_names: pd.Index,
    niche_composition_key: str,
    expression_embedding_niche_key: str,
    inplace: bool,
    min_var_name_ratio: float = 0.8,
) -> AnnData | None:
    r"""
    Pad and sort anndata to match reference var names.

    Also covers \alpha and \eta niche features.
    """
    intersection_genes = adata.var_names.intersection(reference_var_names)
    inter_len_genes = len(intersection_genes)
    if inter_len_genes == 0:
        raise ValueError(
            "No reference var names found in query data. "
            "Please rerun with return_reference_var_names=True "
            "to see reference var names."
        )

    ratio = inter_len_genes / len(reference_var_names)
    logger.info(f"Found {ratio * 100}% reference vars in query data.")
    if ratio < min_var_name_ratio:
        warnings.warn(
            f"Query data contains less than {min_var_name_ratio:.0%} of reference "
            "var names. This may result in poor performance.",
            UserWarning,
            stacklevel=settings.warnings_stacklevel,
        )

    missing_in_query = reference_label_names.difference(query_label_names)
    extra_in_query = query_label_names.difference(reference_label_names)

    if len(missing_in_query) > 0:
        logger.info(f"Labels not observed in query: {sorted(missing_in_query.tolist())}")
    if len(extra_in_query) > 0:
        raise ValueError(
            f"Label(s) observed in query but not in reference: {sorted(extra_in_query.tolist())}"
        )

    # pad alpha composition if needed
    if len(missing_in_query) > 0:
        # add missing columns
        for ct in missing_in_query:
            adata.obsm[niche_composition_key][ct] = 0.0
    # sort columns anyway to match reference
    adata.obsm[niche_composition_key] = adata.obsm[niche_composition_key][reference_label_names]
    pad_and_sorted_query_label_names = pd.Index(adata.obsm[niche_composition_key].columns)
    assert pad_and_sorted_query_label_names.equals(reference_label_names), (
        "Error when sorting query label names to match reference."
    )

    # pad eta niche activation if needed
    cell_type_to_int = adata.uns["cell_type_to_int"]
    if len(missing_in_query) > 0:
        # Update with missing labels: add zero embeddings
        for ct in missing_in_query:
            cell_type_to_int[ct] = len(cell_type_to_int)
            # Add zeros for this new cell type
            zeros = np.zeros(
                (adata.n_obs, 1, adata.obsm[expression_embedding_niche_key].shape[-1])
            )
            z_niche = np.concatenate([adata.obsm[expression_embedding_niche_key], zeros], axis=1)

    else:
        z_niche = adata.obsm[expression_embedding_niche_key]
    # reorder axis 1 to match reference_label_names
    idxs = np.array([cell_type_to_int[ct] for ct in reference_label_names])
    adata.obsm[expression_embedding_niche_key] = z_niche[:, idxs, :]

    genes_to_add = reference_var_names.difference(adata.var_names)
    needs_padding = len(genes_to_add) > 0
    if needs_padding:
        padding_mtx = csr_matrix(np.zeros((adata.n_obs, len(genes_to_add))))
        adata_padding = AnnData(
            X=padding_mtx.copy(),
            layers={layer: padding_mtx.copy() for layer in adata.layers},
        )
        adata_padding.var_names = genes_to_add
        adata_padding.obs_names = adata.obs_names
        # Concatenate object
        adata_out = anndata.concat(
            [adata, adata_padding],
            axis=1,
            join="outer",
            index_unique=None,
            merge="unique",
        )
    else:
        adata_out = adata

    # also covers the case when new adata has more var names than old
    if not reference_var_names.equals(adata_out.var_names):
        adata_out._inplace_subset_var(reference_var_names)

    if inplace:
        if adata_out is not adata:
            adata._init_as_actual(adata_out)
    else:
        return adata_out
