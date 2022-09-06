from functools import wraps
from typing import Callable

import anndata

from scvi.data._utils import _get_latent_adata_type


def unsupported_in_latent_mode(fn: Callable) -> Callable:
    @wraps(fn)
    def wrapper(self, *args, **kwargs):
        adata = None
        if len(args) > 0:
            assert isinstance(args[0], anndata.AnnData)
            adata = args[0]
        if adata is None and len(kwargs) > 0:
            adata = kwargs.get("adata", None)
        adata = self._validate_anndata(adata)
        if _get_latent_adata_type(adata) is not None:
            raise ValueError(
                f"Latent mode currently not supported for the {fn.__qualname__} function."
            )
        return fn(self, *args, **kwargs)

    return wrapper
