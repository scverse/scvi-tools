from typing import Optional

import jax
import jax.numpy as jnp


def _one_hot_jax(indexes: jnp.ndarray, n_classes: Optional[int] = None) -> jnp.ndarray:
    """
    One-hot encode a tensor of category indexes.

    Parameters
    ----------
    indexes
        A :class:`~jax.numpy.ndarray` of shape `(n_samples,)` or `(n_samples, 1)`
        containing the category index for each sample.
    n_classes
        The number of categories. If `None`, the number of categories is inferred from
        the maximum value in `indexes`.

    Returns
    -------
    one_hot
        A :class:`~jax.numpy.ndarray` of shape `(n_samples, n_classes)` containing the
        one-hot encoding of `indexes`.
    """
    indexes = jnp.reshape(indexes, (-1,))
    n_classes = n_classes or jnp.max(indexes) + 1
    return jax.nn.one_hot(indexes, n_classes, dtype=jnp.float32)
