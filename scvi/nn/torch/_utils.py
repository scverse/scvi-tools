from typing import Optional

import torch


def _one_hot_torch(
    indexes: torch.Tensor, n_classes: Optional[int] = None
) -> torch.Tensor:
    """
    One-hot encode a tensor of category indexes.

    Parameters
    ----------
    indexes
        A :class:`~torch.Tensor` of shape `(n_samples,)` or `(n_samples, 1)` containing
        the category index for each sample.
    n_classes
        The number of categories. If `None`, the number of categories is inferred from
        the maximum value in `indexes`.

    Returns
    -------
    one_hot
        A :class:`~torch.Tensor` of shape `(n_samples, n_classes)` containing the one-hot
        encoding of `indexes`.
    """
    indexes = indexes.view(-1).type(torch.long)
    n_classes = n_classes or -1
    one_hot = torch.nn.functional.one_hot(indexes, num_classes=n_classes)
    one_hot = one_hot.type(torch.float32).to(indexes.device)
    return one_hot
