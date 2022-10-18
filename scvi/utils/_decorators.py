from functools import wraps
from typing import Callable

import flax.linen as nn


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


def flax_configure(cls: nn.Module) -> Callable:
    """Decorator to raise an error if the model is in latent mode."""
    original_init = cls.__init__

    @wraps(original_init)
    def init(self, *args, **kwargs):
        self.configure()
        original_init(self, *args, **kwargs)

    cls.__init__ = init
    return cls
