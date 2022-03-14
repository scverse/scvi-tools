from copy import deepcopy
from typing import Optional

import numpy as np
from anndata import AnnData

from scvi import REGISTRY_KEYS

from . import _constants
from ._manager import AnnDataManager
from .fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LabelsWithUnlabeledObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
    ProteinObsmField,
)

LEGACY_REGISTRY_KEY_MAP = {
    "X": REGISTRY_KEYS.X_KEY,
    "batch_indices": REGISTRY_KEYS.BATCH_KEY,
    "labels": REGISTRY_KEYS.LABELS_KEY,
    "cat_covs": REGISTRY_KEYS.CAT_COVS_KEY,
    "cont_covs": REGISTRY_KEYS.CONT_COVS_KEY,
    "protein_expression": REGISTRY_KEYS.PROTEIN_EXP_KEY,
    "ind_x": REGISTRY_KEYS.INDICES_KEY,
}


def registry_from_setup_dict(
    setup_dict: dict, unlabeled_category: Optional[str] = None
) -> dict:
    """
    Converts old setup dict format to new registry dict format.

    Only to be used for backwards compatibility when loading setup dictionaries for models.
    Takes old hard-coded setup dictionary structure and fills in the analogous registry structure.

    Parameters
    ----------
    setup_dict
        Setup dictionary created after registering an AnnData with former ``setup_anndata`` implementation.
    """
    registry = {
        _constants._SCVI_VERSION_KEY: setup_dict[_constants._SCVI_VERSION_KEY],
        _constants._FIELD_REGISTRIES_KEY: {},
    }
    data_registry = setup_dict[_constants._DATA_REGISTRY_KEY]
    categorical_mappings = setup_dict["categorical_mappings"]
    summary_stats = setup_dict[_constants._SUMMARY_STATS_KEY]
    field_registries = registry[_constants._FIELD_REGISTRIES_KEY]
    for (
        registry_key,
        adata_mapping,
    ) in data_registry.items():  # Note: this does not work for empty fields.
        if registry_key not in LEGACY_REGISTRY_KEY_MAP:
            continue
        new_registry_key = LEGACY_REGISTRY_KEY_MAP[registry_key]

        attr_name = adata_mapping[_constants._DR_ATTR_NAME]
        attr_key = adata_mapping[_constants._DR_ATTR_KEY]

        field_registries[new_registry_key] = {
            _constants._DATA_REGISTRY_KEY: adata_mapping,
            _constants._STATE_REGISTRY_KEY: dict(),
            _constants._SUMMARY_STATS_KEY: dict(),
        }
        field_registry = field_registries[new_registry_key]
        field_state_registry = field_registry[_constants._STATE_REGISTRY_KEY]
        field_summary_stats = field_registry[_constants._SUMMARY_STATS_KEY]

        if attr_name in (_constants._ADATA_ATTRS.X, _constants._ADATA_ATTRS.LAYERS):
            field_state_registry[LayerField.N_CELLS_KEY] = summary_stats["n_cells"]
            field_state_registry[LayerField.N_VARS_KEY] = summary_stats["n_vars"]
            field_summary_stats.update(field_state_registry)
        elif attr_name == _constants._ADATA_ATTRS.OBS:
            categorical_mapping = categorical_mappings[attr_key]
            field_state_registry[
                CategoricalObsField.CATEGORICAL_MAPPING_KEY
            ] = categorical_mapping["mapping"]
            if new_registry_key == REGISTRY_KEYS.BATCH_KEY:
                field_summary_stats[f"n_{new_registry_key}"] = summary_stats["n_batch"]
            elif new_registry_key == REGISTRY_KEYS.LABELS_KEY:
                field_summary_stats[f"n_{new_registry_key}"] = summary_stats["n_labels"]
                if unlabeled_category is not None:
                    field_state_registry[
                        LabelsWithUnlabeledObsField.UNLABELED_CATEGORY
                    ] = unlabeled_category
        elif attr_name == _constants._ADATA_ATTRS.OBSM:
            if new_registry_key == REGISTRY_KEYS.CONT_COVS_KEY:
                columns = setup_dict["extra_continuous_keys"].copy()
                field_state_registry[NumericalJointObsField.COLUMNS_KEY] = columns
                field_summary_stats[f"n_{new_registry_key}"] = columns.shape[0]
            elif new_registry_key == REGISTRY_KEYS.CAT_COVS_KEY:
                extra_categoricals_mapping = deepcopy(setup_dict["extra_categoricals"])
                field_state_registry.update(deepcopy(setup_dict["extra_categoricals"]))
                field_summary_stats[f"n_{new_registry_key}"] = len(
                    extra_categoricals_mapping["keys"]
                )
            elif new_registry_key == REGISTRY_KEYS.PROTEIN_EXP_KEY:
                field_state_registry[ProteinObsmField.COLUMN_NAMES_KEY] = setup_dict[
                    "protein_names"
                ].copy()
                if "totalvi_batch_mask" in setup_dict:
                    field_state_registry[
                        ProteinObsmField.PROTEIN_BATCH_MASK
                    ] = setup_dict["totalvi_batch_mask"].copy()
                field_summary_stats[f"n_{new_registry_key}"] = len(
                    setup_dict["protein_names"]
                )
    return registry


def manager_from_setup_dict(
    cls,
    adata: AnnData,
    setup_dict: dict,
    unlabeled_category: Optional[str] = None,
    **transfer_kwargs,
) -> AnnDataManager:
    """
    Creates an :class:`~scvi.data.AnnDataManager` given only a scvi-tools setup dictionary.

    Only to be used for backwards compatibility when loading setup dictionaries for models.
    Infers the AnnDataField instances used to define the :class:`~scvi.data.AnnDataManager` instance,
    then uses the :meth:`~scvi.data.AnnDataManager.transfer_fields` method to register the new AnnData object.

    Parameters
    ----------
    adata
        AnnData object to be registered.
    setup_dict
        Setup dictionary created after registering an AnnData with former ``setup_anndata`` implementation.
    **kwargs
        Keyword arguments to modify transfer behavior.
    """
    fields = []
    setup_args = dict()
    data_registry = setup_dict[_constants._DATA_REGISTRY_KEY]
    categorical_mappings = setup_dict["categorical_mappings"]
    for registry_key, adata_mapping in data_registry.items():
        if registry_key not in LEGACY_REGISTRY_KEY_MAP:
            continue
        new_registry_key = LEGACY_REGISTRY_KEY_MAP[registry_key]

        field = None
        attr_name = adata_mapping[_constants._DR_ATTR_NAME]
        attr_key = adata_mapping[_constants._DR_ATTR_KEY]
        if attr_name == _constants._ADATA_ATTRS.X:
            field = LayerField(REGISTRY_KEYS.X_KEY, None)
            setup_args["layer"] = None
        elif attr_name == _constants._ADATA_ATTRS.LAYERS:
            field = LayerField(REGISTRY_KEYS.X_KEY, attr_key)
            setup_args["layer"] = attr_key
        elif attr_name == _constants._ADATA_ATTRS.OBS:
            if new_registry_key in {REGISTRY_KEYS.BATCH_KEY, REGISTRY_KEYS.LABELS_KEY}:
                original_key = categorical_mappings[attr_key]["original_key"]
                if (
                    unlabeled_category is not None
                    and new_registry_key == REGISTRY_KEYS.LABELS_KEY
                ):
                    field = LabelsWithUnlabeledObsField(
                        new_registry_key, original_key, unlabeled_category
                    )
                else:
                    field = CategoricalObsField(new_registry_key, original_key)
                setup_args[f"{new_registry_key}_key"] = original_key
            elif new_registry_key == REGISTRY_KEYS.INDICES_KEY:
                adata.obs[attr_key] = np.arange(adata.n_obs).astype("int64")
                field = NumericalObsField(new_registry_key, attr_key)
        elif attr_name == _constants._ADATA_ATTRS.OBSM:
            if new_registry_key == REGISTRY_KEYS.CONT_COVS_KEY:
                obs_keys = setup_dict["extra_continuous_keys"]
                field = NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, obs_keys)
                setup_args["continuous_covariate_keys"] = obs_keys
            elif new_registry_key == REGISTRY_KEYS.CAT_COVS_KEY:
                obs_keys = setup_dict["extra_categoricals"]["keys"]
                field = CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, obs_keys)
                setup_args["categorical_covariate_keys"] = obs_keys
            elif new_registry_key == REGISTRY_KEYS.PROTEIN_EXP_KEY:
                protein_names = setup_dict["protein_names"]
                adata.uns["_protein_names"] = protein_names
                field = ProteinObsmField(
                    REGISTRY_KEYS.PROTEIN_EXP_KEY,
                    attr_key,
                    use_batch_mask=True,
                    batch_key="_scvi_batch",
                    colnames_uns_key="_protein_names",
                )
                setup_args["protein_expression_obsm_key"] = attr_key
                setup_args["protein_names_uns_key"] = "_protein_names"
            else:
                raise NotImplementedError(
                    f"Unrecognized .obsm attribute {attr_key} registered as {new_registry_key}. Backwards compatibility unavailable."
                )
        else:
            raise NotImplementedError(
                f"Backwards compatibility for attribute {attr_name} is not implemented."
            )
        fields.append(field)

    setup_method_args = {
        _constants._MODEL_NAME_KEY: cls.__name__,
        _constants._SETUP_ARGS_KEY: setup_args,
    }
    adata_manager = AnnDataManager(fields=fields, setup_method_args=setup_method_args)

    source_registry = registry_from_setup_dict(
        setup_dict, unlabeled_category=unlabeled_category
    )
    adata_manager.register_fields(
        adata, source_registry=source_registry, **transfer_kwargs
    )
    return adata_manager
