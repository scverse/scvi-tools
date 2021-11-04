from anndata import AnnData

from . import _constants
from ._fields import CategoricalObsField, LayerField
from ._manager import AnnDataManager


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
        Setup dictionary created after registering an AnnData using an AnnDataManager object.
    **kwargs
        Keyword arguments to modify transfer behavior.
    """
    source_adata_manager = AnnDataManager()
    data_registry = setup_dict[_constants._DATA_REGISTRY_KEY]
    categorical_mappings = setup_dict[_constants._CATEGORICAL_MAPPINGS_KEY]
    for registry_key, adata_mapping in data_registry.items():
        field = None
        attr_name = adata_mapping[_constants._DR_ATTR_NAME]
        attr_key = adata_mapping[_constants._DR_ATTR_KEY]
        if attr_name == _constants._ADATA_ATTRS.X:
            field = LayerField(registry_key, None)
        elif attr_name == _constants._ADATA_ATTRS.LAYERS:
            field = LayerField(registry_key, attr_key)
        elif attr_name == _constants._ADATA_ATTRS.OBS:
            original_key = categorical_mappings[attr_key][_constants._CM_ORIGINAL_KEY]
            field = CategoricalObsField(registry_key, original_key)
        else:
            raise NotImplementedError(
                f"Backwards compatibility for attribute {attr_name} is not implemented yet."
            )
        source_adata_manager.add_field(field)
    return source_adata_manager.transfer_setup(
        adata, source_setup_dict=setup_dict, **transfer_kwargs
    )
