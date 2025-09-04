from collections.abc import Callable

from scvi.utils import dependencies


@dependencies("jax")
def device_selecting_PRNGKey(use_cpu: bool = True) -> Callable:
    """Returns a PRNGKey that is either on CPU or GPU."""
    # if key is generated on CPU, model params will be on CPU
    import jax
    from jax import random

    if use_cpu is True:

        def key(i: int):
            return jax.device_put(random.PRNGKey(i), jax.devices("cpu")[0])
    else:
        # dummy function
        def key(i: int):
            return random.PRNGKey(i)

    return key
