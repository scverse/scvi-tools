from abc import ABC, abstractmethod
from typing import Optional, Type, Union

import numpy as np
import pandas as pd
import rich

from scvi._types import AnnOrMuData
from scvi.data import _constants
from scvi.data._utils import get_anndata_attribute


class BaseAnnDataField(ABC):
    """
    Abstract class for a single AnnData/MuData field.

    A Field class defines how scvi-tools will map a data field used by a model
    to an attribute in an AnnData/MuData object.
    """

    @property
    @abstractmethod
    def registry_key(self) -> str:
        """The key that is referenced by models via a data loader."""

    @property
    @abstractmethod
    def attr_name(self) -> str:
        """The name of the AnnData attribute where the data is stored."""

    @property
    @abstractmethod
    def attr_key(self) -> Optional[str]:
        """The key of the data field within the relevant AnnData attribute."""

    @property
    def mod_key(self) -> Optional[str]:
        """The modality key of the data field within the MuData (if applicable)."""
        return None

    @property
    @abstractmethod
    def is_empty(self) -> bool:
        """
        Returns True if the field is empty as a function of its kwargs.

        A field can be empty if it is composed of a collection of variables, and for a given
        instance of a model, the collection is empty. If empty, the field will be omitted from
        the registry, but included in the summary stats dictionary.
        """

    @abstractmethod
    def validate_field(self, adata: AnnOrMuData) -> None:
        """Validates whether an AnnData/MuData object is compatible with this field definition."""

    @abstractmethod
    def register_field(self, adata: AnnOrMuData) -> dict:
        """
        Sets up the AnnData/MuData object and creates a mapping for scvi-tools models to use.

        Returns
        -------
        dict
            A dictionary containing any additional state required for scvi-tools models not
            stored directly on the AnnData/MuData object.
        """
        self.validate_field(adata)
        return dict()

    @abstractmethod
    def transfer_field(
        self, state_registry: dict, adata_target: AnnOrMuData, **kwargs
    ) -> dict:
        """
        Takes an existing scvi-tools setup dictionary and transfers the same setup to the target AnnData.

        Used when one is running a pretrained model on a new AnnData object, which
        requires the mapping from the original data to be applied to the new AnnData object.

        Parameters
        ----------
        state_registry
            state_registry dictionary created after registering an AnnData using an :class:`~scvi.data.AnnDataManager` object.
        adata_target
            AnnData/MuData object that is being registered.
        **kwargs
            Keyword arguments to modify transfer behavior.

        Returns
        -------
        dict
            A dictionary containing any additional state required for scvi-tools models not
            stored directly on the AnnData object.
        """
        return dict()

    @abstractmethod
    def get_summary_stats(self, state_registry: dict) -> dict:
        """
        Returns a dictionary comprising of summary statistics relevant to the field.

        Parameters
        ----------
        state_registry
            Dictionary returned by :meth:`~scvi.data.fields.BaseAnnDataField.register_field`.
            Summary stats should always be a function of information stored in this dictionary.

        Returns
        -------
        summary_stats_dict
            The dictionary is of the form {summary_stat_name: summary_stat_value}.
            This mapping is then combined with the mappings of other fields to make up
            the summary stats mapping.
        """

    @abstractmethod
    def view_state_registry(self, state_registry: dict) -> Optional[rich.table.Table]:
        """
        Returns a :class:`rich.table.Table` summarizing a state registry produced by this field.

        Parameters
        ----------
        state_registry
            Dictionary returned by :meth:`~scvi.data.fields.BaseAnnDataField.register_field`.
            Printed summary should always be a function of information stored in this dictionary.

        Returns
        -------
        state_registry_summary
            Optional :class:`rich.table.Table` summarizing the ``state_registry``.
        """

    def get_field_data(self, adata: AnnOrMuData) -> Union[np.ndarray, pd.DataFrame]:
        """Returns the requested data as determined by the field for a given AnnData/MuData object."""
        if self.is_empty:
            raise AssertionError(f"The {self.registry_key} field is empty.")
        return get_anndata_attribute(
            adata, self.attr_name, self.attr_key, mod_key=self.mod_key
        )

    def get_data_registry(self) -> dict:
        """
        Returns a nested dictionary which describes the mapping to the data field.

        The dictionary is of the form {"mod_key": mod_key, "attr_name": attr_name, "attr_key": attr_key}.
        This mapping is then combined with the mappings of other fields to make up the data registry.
        """
        if self.is_empty:
            return dict()

        data_registry = {
            _constants._DR_ATTR_NAME: self.attr_name,
            _constants._DR_ATTR_KEY: self.attr_key,
        }
        if self.mod_key is not None:
            data_registry[_constants._DR_MOD_KEY] = self.mod_key
        return data_registry


# Convenience type
AnnDataField = Type[BaseAnnDataField]
