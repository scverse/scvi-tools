from collections.abc import Callable
from functools import wraps

from scvi.data._constants import ADATA_MINIFY_TYPE


def unsupported_if_adata_minified(fn: Callable) -> Callable:
    """Decorator to raise an error if the model's `adata` is minified."""

    @wraps(fn)
    def wrapper(self, *args, **kwargs):
        if getattr(self, "minified_data_type", None) == ADATA_MINIFY_TYPE.LATENT_POSTERIOR:
            raise ValueError(
                f"The {fn.__qualname__} function currently does not support minified data."
            )
        return fn(self, *args, **kwargs)

    return wrapper
