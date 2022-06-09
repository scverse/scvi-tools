from typing import Callable

import jax
from jax import random


def device_selecting_PRNGKey(use_gpu: bool) -> Callable:
    # if key is generated on CPU, model params will be on CPU
    # we have to pay the price of a JIT compilation though
    if use_gpu is False:
        key = jax.jit(lambda i: random.PRNGKey(i), backend="cpu")
    else:
        # dummy function
        def key(i: int):
            return random.PRNGKey(i)

    return key
