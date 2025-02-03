from __future__ import annotations

import logging
import warnings
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import torch
from rich import print

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data._constants import _ADATA_MINIFY_TYPE_UNS_KEY, ADATA_MINIFY_TYPE
from scvi.data._utils import _get_adata_minify_type
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
    ObsmField,
    StringUnsField,
)
from scvi.model._utils import _init_library_size
from scvi.model.base import (
    ArchesMixin,
    BaseMinifiedModeModelClass,
    EmbeddingMixin,
    UnsupervisedTrainingMixin,
)

# RNASeqMixin,
# VAEMixin,
from scvi.model.utils import get_minified_adata_scrna
from scvi.utils import setup_anndata_dsp

from ._constants import NICHEVI_REGISTRY_KEYS
from ._module import nicheVAE
from ._rnamixin import NicheRNASeqMixin

# from ._training_mixin import UnsupervisedTrainingMixin
from ._vaemixin import NicheVAEMixin

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData

    from scvi._types import MinifiedDataType
    from scvi.data.fields import (
        BaseAnnDataField,
    )

_SCVI_LATENT_QZM = "_scvi_latent_qzm"
_SCVI_LATENT_QZV = "_scvi_latent_qzv"
_SCVI_OBSERVED_LIB_SIZE = "_scvi_observed_lib_size"

logger = logging.getLogger(__name__)


class nicheSCVI(
    EmbeddingMixin,
    NicheRNASeqMixin,
    NicheVAEMixin,
    ArchesMixin,
    UnsupervisedTrainingMixin,
    BaseMinifiedModeModelClass,
):
    """single-cell Variational Inference :cite:p:`Lopez18`.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.SCVI.setup_anndata`. If
        ``None``, then the underlying module will not be initialized until training, and a
        :class:`~lightning.pytorch.core.LightningDataModule` must be passed in during training
        (``EXPERIMENTAL``).
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
    >>> scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.model.SCVI(adata)
    >>> vae.train()
    >>> adata.obsm["X_scVI"] = vae.get_latent_representation()
    >>> adata.obsm["X_normalized_scVI"] = vae.get_normalized_expression()

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/quick_start/api_overview`
    2. :doc:`/tutorials/notebooks/scrna/harmonization`
    3. :doc:`/tutorials/notebooks/scrna/scarches_scvi_tools`
    4. :doc:`/tutorials/notebooks/scrna/scvi_in_R`

    See Also
    --------
    :class:`~scvi.module.VAE`
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
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb",
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
            "nicheVI model with the following parameters: \n"
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

    def preprocessing_anndata(
        adata: AnnData,
        k_nn: int = 20,
        sample_key: str | None = None,
        labels_key: str = "cell_type",
        cell_coordinates_key: str = "spatial",
        expression_embedding_key: str = "X_scVI",
        expression_embedding_niche_key: str = "X_scVI_niche",
        niche_composition_key: str = "niche_composition",
        niche_indexes_key: str = "niche_indexes",
        niche_distances_key: str = "niche_distances",
        log1p: bool = False,
    ):
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
        ############################
        sample_key: str | None = None,
        labels_key: str = "cell_type",
        cell_coordinates_key: str = "spatial",
        expression_embedding_key: str = "X_scVI",
        expression_embedding_niche_key: str = "X_scVI_niche",
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
            CategoricalObsField(NICHEVI_REGISTRY_KEYS.SAMPLE_KEY, sample_key),
            ObsmField(NICHEVI_REGISTRY_KEYS.NICHE_COMPOSITION_KEY, niche_composition_key),
            ObsmField(NICHEVI_REGISTRY_KEYS.CELL_COORDINATES_KEY, cell_coordinates_key),
            ObsmField(NICHEVI_REGISTRY_KEYS.NICHE_INDEXES_KEY, niche_indexes_key),
            ObsmField(NICHEVI_REGISTRY_KEYS.NICHE_DISTANCES_KEY, niche_distances_key),
            ObsmField(NICHEVI_REGISTRY_KEYS.Z1_MEAN_KEY, expression_embedding_key),
            ObsmField(NICHEVI_REGISTRY_KEYS.Z1_MEAN_CT_KEY, expression_embedding_niche_key),
        ]
        # register new fields if the adata is minified
        adata_minify_type = _get_adata_minify_type(adata)
        if adata_minify_type is not None:
            anndata_fields += cls._get_fields_for_adata_minification(adata_minify_type)
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    @staticmethod
    def _get_fields_for_adata_minification(
        minified_data_type: MinifiedDataType,
    ) -> list[BaseAnnDataField]:
        """Return the fields required for adata minification of the given minified_data_type."""
        if minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR:
            fields = [
                ObsmField(
                    REGISTRY_KEYS.LATENT_QZM_KEY,
                    _SCVI_LATENT_QZM,
                ),
                ObsmField(
                    REGISTRY_KEYS.LATENT_QZV_KEY,
                    _SCVI_LATENT_QZV,
                ),
                NumericalObsField(
                    REGISTRY_KEYS.OBSERVED_LIB_SIZE,
                    _SCVI_OBSERVED_LIB_SIZE,
                ),
            ]
        else:
            raise NotImplementedError(f"Unknown MinifiedDataType: {minified_data_type}")
        fields.append(
            StringUnsField(
                REGISTRY_KEYS.MINIFY_TYPE_KEY,
                _ADATA_MINIFY_TYPE_UNS_KEY,
            ),
        )
        return fields

    def minify_adata(
        self,
        minified_data_type: MinifiedDataType = ADATA_MINIFY_TYPE.LATENT_POSTERIOR,
        use_latent_qzm_key: str = "X_latent_qzm",
        use_latent_qzv_key: str = "X_latent_qzv",
    ) -> None:
        """Minifies the model's adata.

        Minifies the adata, and registers new anndata fields: latent qzm, latent qzv, adata uns
        containing minified-adata type, and library size.
        This also sets the appropriate property on the module to indicate that the adata is
        minified.

        Parameters
        ----------
        minified_data_type
            How to minify the data. Currently only supports `latent_posterior_parameters`.
            If minified_data_type == `latent_posterior_parameters`:

            * the original count data is removed (`adata.X`, adata.raw, and any layers)
            * the parameters of the latent representation of the original data is stored
            * everything else is left untouched
        use_latent_qzm_key
            Key to use in `adata.obsm` where the latent qzm params are stored
        use_latent_qzv_key
            Key to use in `adata.obsm` where the latent qzv params are stored

        Notes
        -----
        The modification is not done inplace -- instead the model is assigned a new (minified)
        version of the adata.
        """
        # TODO(adamgayoso): Add support for a scenario where we want to cache the latent posterior
        # without removing the original counts.
        if minified_data_type != ADATA_MINIFY_TYPE.LATENT_POSTERIOR:
            raise NotImplementedError(f"Unknown MinifiedDataType: {minified_data_type}")

        if self.module.use_observed_lib_size is False:
            raise ValueError("Cannot minify the data if `use_observed_lib_size` is False")

        minified_adata = get_minified_adata_scrna(self.adata, minified_data_type)
        minified_adata.obsm[_SCVI_LATENT_QZM] = self.adata.obsm[use_latent_qzm_key]
        minified_adata.obsm[_SCVI_LATENT_QZV] = self.adata.obsm[use_latent_qzv_key]
        counts = self.adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)
        minified_adata.obs[_SCVI_OBSERVED_LIB_SIZE] = np.squeeze(np.asarray(counts.sum(axis=1)))
        self._update_adata_and_manager_post_minification(minified_adata, minified_data_type)
        self.module.minified_data_type = minified_data_type

    @torch.inference_mode()
    def predict_neighborhood(
        self,
        adata: AnnData | None = None,
        indices: np.ndarray | None = None,
        batch_size: int | None = 1024,
    ):
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

            ct_prediction.append(predicted_ct_prob.concentration.detach().cpu())

        return torch.cat(ct_prediction).numpy()

    @torch.inference_mode()
    def predict_niche_activation(
        self,
        adata: AnnData | None = None,
        indices: np.ndarray | None = None,
        batch_size: int | None = 1024,
    ):
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

    @torch.inference_mode()
    def get_niche_attention(
        self,
        adata: AnnData | None = None,
        indices: np.ndarray | None = None,
        batch_size: int = 1024,
    ) -> np.ndarray:
        """Description

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
        niche_attention
            Attention weights for each cell in the dataset.
        """
        self._check_if_trained(warn=False)

        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        if self.module.n_heads is None:
            raise ValueError(
                "The model was not trained with the attention_decoder parameter set to True. "
                "Please retrain the model with the attention_decoder parameter set to True."
            )

        attention_weights = []
        for tensors in scdl:
            inference_inputs = self.module._get_inference_input(tensors)
            outputs = self.module.inference(**inference_inputs)

            batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
            decoder_input = outputs["qz"].loc

            # put batch_index in the same device as decoder_input
            batch_index = batch_index.to(decoder_input.device)

            niche_mean, niche_variance, niche_attention = self.module.niche_decoder(
                decoder_input,
                batch_index,
            )

            attention_weights.append(niche_attention.detach().cpu())

        return torch.cat(attention_weights).numpy()

    def get_cell_type_attention(
        self,
        adata: AnnData | None = None,
        attention_key: str = "attention_weights",
        cell_type_key: str = "cell_type",
        compute_attention: bool = True,
    ):
        if self.module.n_heads is None:
            raise ValueError(
                "The model was not trained with the attention_decoder parameter set to True. "
                "Please retrain the model with the attention_decoder parameter set to True."
            )

        self._check_if_trained(warn=False)

        adata = self._validate_anndata(adata)

        cell_types = adata.obs[cell_type_key].unique().tolist()
        cell_type_to_int = {cell_types[i]: i for i in range(len(cell_types))}

        if compute_attention:
            attention_weights = self.get_niche_attention(adata=adata)

        else:
            attention_weights = adata.obsm[attention_key]

        attention_weights = attention_weights[:, 1:, 1:]

        token_attention_weights = {
            token_name: attention_weights[:, token_idx, :]
            for token_name, token_idx in cell_type_to_int.items()
        }

        return token_attention_weights


def get_niche_indexes(
    adata: AnnData,
    sample_key: str,
    niche_indexes_key: str,
    cell_coordinates_key: str,
    k_nn: int,
    niche_distances_key: str | None = None,
):
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

    print("[bold cyan]Saved niche_indexes and niche_distances in adata.obsm[/bold cyan]")

    return None


def get_neighborhood_composition(
    adata: AnnData,
    cell_type_column: str,
    indices_key: str = "niche_indexes",
    niche_composition_key: str = "niche_composition",
):
    n_cell_types = len(adata.obs[cell_type_column].unique())  # number of cell types
    adata.obsm[niche_composition_key] = np.zeros(
        (adata.n_obs, n_cell_types)
    )  # for each cell, store the composition of its neighborhood as a convex vector of cell type
    # proportions

    indices = adata.obsm[indices_key].astype(int)

    cell_types = adata.obs[cell_type_column].unique().tolist()
    cell_type_to_int = {cell_types[i]: i for i in range(len(cell_types))}

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

    print("[bold green]Saved niche_composition in adata.obsm[/bold green]")

    return None


def get_average_latent_per_celltype(
    adata: AnnData,
    labels_key: str,
    niche_indexes_key: str,
    latent_mean_key: str | None = None,
    latent_mean_ct_key: str = "qz1_m_niche_ct",
    log1p: bool = False,
):
    # for each cell, take the average of the latent space for each label, namely the label-averaged
    # latent_mean obsm

    if latent_mean_key is None:
        adata.obsm["qz1_m_niche_ct"] = np.empty(
            (adata.n_obs, adata.obsm[latent_mean_key].shape[1])
        )

        return None

    n_cells = adata.n_obs
    n_cell_types = len(adata.obs[labels_key].unique())
    n_latent_z1 = adata.obsm[latent_mean_key].shape[1]
    niche_indexes = adata.obsm[niche_indexes_key]

    if log1p:
        z1_mean_niches = np.log1p(adata.obsm[latent_mean_key])[niche_indexes]

    else:
        z1_mean_niches = adata.obsm[latent_mean_key][niche_indexes]

    cell_types = adata.obs[labels_key].unique().tolist()

    cell_type_to_int = {cell_types[i]: i for i in range(len(cell_types))}
    integer_vector = np.vectorize(cell_type_to_int.get)(adata.obs[labels_key])

    # For each cell, get the cell types of its neighbors (as integers)
    cell_types_in_the_neighborhood = np.vstack(
        [integer_vector[niche_indexes[cell, :]] for cell in range(n_cells)]
    )

    dict_of_cell_type_indices = {}

    for cell_type, cell_type_idx in cell_type_to_int.items():
        ct_row_indices, ct_col_indices = np.where(
            cell_types_in_the_neighborhood == cell_type_idx
        )  # [1]

        # dict of cells:local index of the cells of cell_type in the neighborhood.
        result_dict = {}
        for row_idx, col_idx in zip(ct_row_indices, ct_col_indices, strict=False):
            result_dict.setdefault(row_idx, []).append(col_idx)

        dict_of_cell_type_indices[cell_type] = result_dict

    # print(dict_of_cell_type_indices)

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

    print("[bold green]Saved qz1_m_niche_ct in adata.obsm[/bold green]")

    return None
