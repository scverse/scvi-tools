from typing import Callable, Optional

import rich
from anndata import AnnData
from mudata import MuData

from scvi._types import AnnOrMuData

from ._base_field import AnnDataField, BaseAnnDataField


class BaseMuDataWrapperClass(BaseAnnDataField):
    """
    A wrapper class that adds MuData support for an AnnDataField.

    The wrapper class makes calls to the ``self.adata_field`` instance using the appropriate
    modality and allows for manipulation at the MuData object level.

    Parameters
    ----------
    mod_key
        Modality key where data is stored. If ``None``, uses the top level MuData object attributes.
    mod_required
        If ``True``, raises ``ValueError`` when ``mod_key`` is ``None``.
    """

    def __init__(
        self, mod_key: Optional[str] = None, mod_required: bool = False
    ) -> None:
        super().__init__()
        if mod_required and mod_key is None:
            raise ValueError(
                f"Modality required for {self.__class__.__name__} but not provided."
            )
        self._mod_key = mod_key
        self._preregister = lambda _self, _mdata: None

    @property
    def adata_field(self) -> AnnDataField:
        """:class"`~scvi.data.fields.AnnDataField` instance that this class instance wraps."""
        return self._adata_field

    @property
    def registry_key(self) -> str:
        """The key that is referenced by models via a data loader."""
        return self.adata_field.registry_key

    @property
    def mod_key(self) -> Optional[str]:
        """The modality key of the data field within the MuData (if applicable)."""
        return self._mod_key

    @property
    def attr_name(self) -> str:
        """The name of the AnnData/MuData attribute where the data is stored."""
        return self.adata_field.attr_name

    @property
    def attr_key(self) -> Optional[str]:
        """The key of the data field within the relevant AnnData/MuData attribute."""
        return self.adata_field.attr_key

    @property
    def is_empty(self) -> bool:  # noqa: D102
        return self.adata_field.is_empty

    def get_modality(self, mdata: MuData) -> AnnOrMuData:
        """Fetches the appropriate modality from the MuData object using ``self.mod_key``."""
        if isinstance(mdata, AnnData):
            raise AssertionError("`get_modality` can only be called on MuData objects.")
        bdata = mdata
        if self.mod_key is not None:
            if self.mod_key not in mdata.mod:
                raise ValueError(f"Modality {self.mod_key} not found in mdata.mod.")
            bdata = mdata.mod[self.mod_key]
        return bdata

    def validate_field(self, mdata: MuData) -> None:
        """Validate the field."""
        if isinstance(mdata, AnnData):
            raise ValueError("`get_modality` can only be called on MuData objects.")
        bdata = self.get_modality(mdata)
        return self.adata_field.validate_field(bdata)

    def preregister(self, mdata: MuData) -> None:
        """
        Function that is called prior to registering fields.

        Function that is be called at the beginning of :meth:`~scvi.data.fields.BaseMuDataWrapperClass.register_field`
        and :meth:`~scvi.data.fields.BaseMuDataWrapperClass.transfer_field`.
        Used when data manipulation is necessary across modalities.
        """
        return self._preregister(self, mdata)

    def register_field(self, mdata: MuData) -> dict:
        """Register the field."""
        self.preregister(mdata)
        bdata = self.get_modality(mdata)
        return self.adata_field.register_field(bdata)

    def transfer_field(
        self, state_registry: dict, mdata_target: MuData, **kwargs
    ) -> dict:
        """Transfer the field."""
        self.preregister(mdata_target)
        bdata_target = self.get_modality(mdata_target)
        return self.adata_field.transfer_field(state_registry, bdata_target, **kwargs)

    def get_summary_stats(self, state_registry: dict) -> dict:
        """Get summary stats."""
        return self.adata_field.get_summary_stats(state_registry)

    def view_state_registry(self, state_registry: dict) -> Optional[rich.table.Table]:
        """View the state registry."""
        return self.adata_field.view_state_registry(state_registry)


def MuDataWrapper(
    adata_field_cls: AnnDataField, preregister_fn: Optional[Callable] = None
) -> AnnDataField:
    """
    Wraps an AnnDataField with :class:`~scvi.data.fields.BaseMuDataWrapperClass`.

    Parameters
    ----------
    adata_field_cls
        AnnDataField class to wrap.
    preregister_fn
        Function that will be called at the beginning of :meth:`~scvi.data.fields.BaseMuDataWrapperClass.register_field`
        and :meth:`~scvi.data.fields.BaseMuDataWrapperClass.transfer_field`.
    """
    if not isinstance(adata_field_cls, type):
        raise ValueError("`adata_field_cls` must be a class, not an instance.")

    def mudata_field_init(
        self, *args, mod_key: Optional[str] = None, mod_required: bool = False, **kwargs
    ):
        BaseMuDataWrapperClass.__init__(
            self, mod_key=mod_key, mod_required=mod_required
        )
        self._adata_field = adata_field_cls(*args, **kwargs)
        if preregister_fn is not None:
            self._preregister = preregister_fn

    return type(
        f"MuData{adata_field_cls.__name__}",
        (BaseMuDataWrapperClass,),
        {
            "__init__": mudata_field_init,
        },
    )
