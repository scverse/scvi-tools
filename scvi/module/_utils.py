import torch

from scvi.nn import one_hot


def iterate(obj, func):
    """Iterates over an object and applies a function to each element."""
    t = type(obj)
    if t is list or t is tuple:
        return t([iterate(o, func) for o in obj])
    else:
        return func(obj) if obj is not None else None


def broadcast_labels(y, *o, n_broadcast=-1):
    """
    Utility for the semi-supervised setting.

    If y is defined(labelled batch) then one-hot encode the labels (no broadcasting needed)
    If y is undefined (unlabelled batch) then generate all possible labels (and broadcast other arguments if not None)
    """
    if not len(o):
        raise ValueError("Broadcast must have at least one reference argument")
    if y is None:
        ys = enumerate_discrete(o[0], n_broadcast)
        new_o = iterate(
            o,
            lambda x: x.repeat(n_broadcast, 1)
            if len(x.size()) == 2
            else x.repeat(n_broadcast),
        )
    else:
        ys = one_hot(y, n_broadcast)
        new_o = o
    return (ys,) + new_o


def enumerate_discrete(x, y_dim):
    """Enumerate discrete variables."""

    def batch(batch_size, label):
        labels = torch.ones(batch_size, 1, device=x.device, dtype=torch.long) * label
        return one_hot(labels, y_dim)

    batch_size = x.size(0)
    return torch.cat([batch(batch_size, i) for i in range(y_dim)])


def masked_softmax(weights, mask, dim=-1, eps=1e-30):
    """
    Computes a softmax of ``weights`` along ``dim`` where ``mask is True``.

    Adds a small ``eps`` term in the numerator and denominator to avoid zero division.
    Taken from: https://discuss.pytorch.org/t/apply-mask-softmax/14212/15.
    Pytorch issue tracked at: https://github.com/pytorch/pytorch/issues/55056.
    """
    weight_exps = torch.exp(weights)
    masked_exps = weight_exps.masked_fill(mask == 0, eps)
    masked_sums = masked_exps.sum(dim, keepdim=True) + eps
    return masked_exps / masked_sums
