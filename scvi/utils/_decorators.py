from functools import wraps
from typing import Callable


def unsupported_in_latent_mode(fn: Callable) -> Callable:
    """Decorator to raise an error if the model is in latent mode."""

    @wraps(fn)
    def wrapper(self, *args, **kwargs):
        if getattr(self, "latent_data_type", None) is not None:
            raise ValueError(
                f"Latent mode currently not supported for the {fn.__qualname__} function."
            )
        return fn(self, *args, **kwargs)

    return wrapper
