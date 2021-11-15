import logging
from typing import List, Optional

import pandas as pd
from anndata import AnnData

from scvi.data.anndata import _constants

from ._base_field import BaseAnnDataField

logger = logging.getLogger(__name__)


class BaseObsmField(BaseAnnDataField):
    """An abstract AnnDataField for .obsm attributes in the AnnData data structure."""

    _attr_name = _constants._ADATA_ATTRS.OBSM

    def __init__(
        self,
        registry_key: str,
    ) -> None:
        super().__init__()
        self._registry_key = registry_key

    @property
    def registry_key(self):
        return self._registry_key

    @property
    def attr_name(self):
        return self._attr_name


class JointObsField(BaseObsmField):
    """
    An abstract AnnDataField for a collection of .obs fields in the AnnData data structure.

    Creates an .obsm field containing each .obs field to be referenced as a whole a model.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    obs_keys
        Sequence of keys to combine to form the obsm field.
    """

    def __init__(self, registry_key: str, obs_keys: Optional[List[str]]) -> None:
        super().__init__(registry_key)
        self._attr_key = f"_scvi_{registry_key}"
        self._columns_key = f"{self.registry_key}_keys"
        self._obs_keys = obs_keys if obs_keys is not None else []
        self._is_empty = len(self.obs_keys) == 0

    def _combine_obs_fields(self, adata: AnnData) -> None:
        adata.obsm[self.attr_key] = pd.concat(
            (adata.obs[key] for key in self.obs_keys), axis=1
        )

    @property
    def obs_keys(self):
        return self._obs_keys

    @property
    def attr_key(self):
        return self._attr_key

    @property
    def is_empty(self) -> bool:
        return self._is_empty


class NumericalJointObsField(JointObsField):
    """
    An AnnDataField for a collection of numerical .obs fields in the AnnData data structure.

    Creates an .obsm field containing each .obs field to be referenced as a whole a model.

    Parameters
    ----------
    registry_key
        Key to register field under in data registry.
    obs_keys
        Sequence of keys to combine to form the obsm field.
    """

    def validate_field(self, adata: AnnData) -> None:
        super().validate_field(adata)
        for obs_key in self._obs_keys:
            assert obs_key in adata.obs, f"{obs_key} not found in adata.obs."

    def register_field(self, adata: AnnData) -> None:
        super().register_field(adata)
        self._combine_obs_fields(adata)
        adata.uns[_constants._SETUP_DICT_KEY][self._columns_key] = self.get_field(
            adata
        ).columns.to_numpy()

    def transfer_field(
        self,
        setup_dict: dict,
        adata_target: AnnData,
        **kwargs,
    ) -> None:
        super().transfer_field(setup_dict, adata_target, **kwargs)
        self.validate_field(adata_target)
        self.register_field(adata_target)

    def compute_summary_stats(self, adata: AnnData) -> dict:
        super().compute_summary_stats(adata)
        n_continuous_covariates = len(self.obs_keys)
        stat_name = f"n_{self.registry_key}"
        return {stat_name: n_continuous_covariates}
