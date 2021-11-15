from abc import ABC, abstractmethod
from typing import Optional

import numpy as np
from anndata import AnnData

from scvi.data.anndata import _constants
from scvi.data.anndata._utils import _get_field


class BaseAnnDataField(ABC):
    """
    Abstract class for a single AnnData field.

    An AnnDataField class defines how scvi-tools will map a data field used by a model
    to an attribute in an AnnData object.
    """

    def __init__(self) -> None:
        super().__init__()

    @property
    @abstractmethod
    def registry_key(self):
        """The key that is referenced by models via a data loader."""
        pass

    @property
    @abstractmethod
    def attr_name(self) -> str:
        """The name of the AnnData attribute where the data is stored (e.g. obs)."""
        pass

    @property
    @abstractmethod
    def attr_key(self) -> Optional[str]:
        """The key of the data field within the relevant AnnData attribute."""
        pass

    @property
    @abstractmethod
    def is_empty(self) -> bool:
        """
        Returns True if the field is empty as a function of its kwargs.

        A field can be empty if it is composed of a collection of variables, and for a given
        instance of a model, the collection is empty. If empty, the field will be omitted from
        the registry, but included in the summary stats dictionary.
        """
        pass

    @abstractmethod
    def validate_field(self, adata: AnnData) -> None:
        """Validates whether an AnnData object is compatible with this field definition."""
        pass

    @abstractmethod
    def register_field(self, adata: AnnData) -> None:
        """Sets up the AnnData object and creates a mapping for scvi-tools models to use."""
        self.validate_field(adata)

    @abstractmethod
    def transfer_field(self, setup_dict: dict, adata_target: AnnData, **kwargs) -> None:
        """
        Takes an existing scvi-tools setup dictionary and transfers the same setup to the target AnnData.

        Used when one is running a pretrained model on a new AnnData object, which
        requires the mapping from the original data to be applied to the new AnnData object.

        Parameters
        ----------
        setup_dict
            Setup dictionary created after registering an AnnData using an AnnDataManager object.
        adata_target
            AnnData object that is being registered.
        **kwargs
            Keyword arguments to modify transfer behavior.
        """
        pass

    def data_registry_mapping(self) -> dict:
        """
        Returns a nested dictionary which describes the mapping to the AnnData data field.

        The dictionary is of the form {registry_key: {"attr_name": attr_name, "attr_key": attr_key}}.
        This mapping is then combined with the mappings of other fields to make up the data registry.
        """
        if self.is_empty:
            return dict()

        return {
            self.registry_key: {
                _constants._DR_ATTR_NAME: self.attr_name,
                _constants._DR_ATTR_KEY: self.attr_key,
            }
        }

    def get_field(self, adata: AnnData) -> np.ndarray:
        """Returns the data field as a NumPy array for a given AnnData object."""
        assert not self.is_empty
        return _get_field(adata, self.attr_name, self.attr_key)

    def compute_summary_stats(self, adata: AnnData) -> dict:
        """Returns a dictionary comprising of summary statistics relevant to the field."""
        return dict()
