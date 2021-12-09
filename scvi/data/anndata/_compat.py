from copy import deepcopy

from anndata import AnnData

from . import _constants
from .fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LayerField,
    NumericalJointObsField,
)
from .manager import AnnDataManager


def registry_from_setup_dict(setup_dict: dict) -> dict:
    """
    Converts old setup dict format to new registry dict format.

    Only to be used for backwards compatibility when loading setup dictionaries for models.
    Takes old hard-coded setup dictionary structure and fills in the analogous registry structure.

    Parameters
    ----------
    setup_dict
        Setup dictionary created after registering an AnnData with former `setup_anndata(...)` implementation.
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
        attr_name = adata_mapping[_constants._DR_ATTR_NAME]
        attr_key = adata_mapping[_constants._DR_ATTR_KEY]

        field_registries[registry_key] = {
            _constants._DATA_REGISTRY_KEY: adata_mapping,
            _constants._STATE_REGISTRY_KEY: dict(),
            _constants._SUMMARY_STATS_KEY: dict(),
        }
        field_registry = field_registries[registry_key]
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
            if attr_key == "_scvi_batch":
                field_summary_stats[f"n_{registry_key}"] = summary_stats["n_batch"]
            elif attr_key == "_scvi_labels":
                field_summary_stats[f"n_{registry_key}"] = summary_stats["n_labels"]
        elif attr_name == _constants._ADATA_ATTRS.OBSM:
            if attr_key == "_scvi_extra_continuous":
                columns = setup_dict["extra_continuous_keys"].copy()
                field_state_registry[NumericalJointObsField.COLUMNS_KEY] = columns
                field_summary_stats[f"n_{registry_key}"] = columns.shape[0]
            elif attr_key == "_scvi_extra_categoricals":
                extra_categoricals_mapping = deepcopy(setup_dict["extra_categoricals"])
                field_state_registry.update(deepcopy(setup_dict["extra_categoricals"]))
                field_summary_stats[f"n_{registry_key}"] = len(
                    extra_categoricals_mapping["keys"]
                )
    return registry


def manager_from_setup_dict(
    adata: AnnData, setup_dict: dict, **transfer_kwargs
) -> AnnDataManager:
    """
    Creates an AnnDataManager given only a scvi-tools setup dictionary.

    Only to be used for backwards compatibility when loading setup dictionaries for models.
    Infers the AnnDataField instances used to define the AnnDataManager instance,
    then uses the `AnnDataManager.transfer_setup(...)` method to register the new AnnData object.

    Parameters
    ----------
    adata
        AnnData object to be registered.
    setup_dict
        Setup dictionary created after registering an AnnData with former `setup_anndata(...)` implementation.
    **kwargs
        Keyword arguments to modify transfer behavior.
    """
    source_adata_manager = AnnDataManager()
    data_registry = setup_dict[_constants._DATA_REGISTRY_KEY]
    categorical_mappings = setup_dict["categorical_mappings"]
    for registry_key, adata_mapping in data_registry.items():
        field = None
        attr_name = adata_mapping[_constants._DR_ATTR_NAME]
        attr_key = adata_mapping[_constants._DR_ATTR_KEY]
        if attr_name == _constants._ADATA_ATTRS.X:
            field = LayerField(registry_key, None)
        elif attr_name == _constants._ADATA_ATTRS.LAYERS:
            field = LayerField(registry_key, attr_key)
        elif attr_name == _constants._ADATA_ATTRS.OBS:
            original_key = categorical_mappings[attr_key]["original_key"]
            field = CategoricalObsField(registry_key, original_key)
        elif attr_name == _constants._ADATA_ATTRS.OBSM:
            if attr_key == "_scvi_extra_continuous":
                obs_keys = setup_dict["extra_continuous_keys"]
                field = NumericalJointObsField(registry_key, obs_keys)
            elif attr_key == "_scvi_extra_categoricals":
                obs_keys = setup_dict["extra_categoricals"]["keys"]
                field = CategoricalJointObsField(registry_key, obs_keys)
            else:
                raise NotImplementedError(
                    f"Unrecognized .obsm attribute {attr_key}. Backwards compatibility unavailable."
                )
        else:
            raise NotImplementedError(
                f"Backwards compatibility for attribute {attr_name} is not implemented yet."
            )
        source_adata_manager.add_field(field)
    source_registry = registry_from_setup_dict(setup_dict)
    return source_adata_manager.transfer_setup(
        adata, source_registry=source_registry, **transfer_kwargs
    )
