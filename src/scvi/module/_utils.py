import torch


def broadcast_labels(o, n_broadcast=-1):
    """Utility for the semi-supervised setting.

    If y is defined(labelled batch) then one-hot encode the labels (no broadcasting needed)
    If y is undefined (unlabelled batch) then generate all possible labels (and broadcast other
    arguments if not None)
    """
    ys_ = torch.nn.functional.one_hot(
        torch.arange(n_broadcast, device=o.device, dtype=torch.long), n_broadcast
    )
    ys = ys_.repeat_interleave(o.size(-2), dim=0)
    if o.ndim == 2:
        new_o = o.repeat(n_broadcast, 1)
    elif o.ndim == 3:
        new_o = o.repeat(1, n_broadcast, 1)
    return ys, new_o


def masked_softmax(weights, mask, dim=-1, eps=1e-30):
    """Computes a softmax of ``weights`` along ``dim`` where ``mask is True``.

    Adds a small ``eps`` term in the numerator and denominator to avoid zero division.
    Taken from: https://discuss.pytorch.org/t/apply-mask-softmax/14212/15.
    Pytorch issue tracked at: https://github.com/pytorch/pytorch/issues/55056.
    """
    weight_exps = torch.exp(weights)
    masked_exps = weight_exps.masked_fill(mask == 0, eps)
    masked_sums = masked_exps.sum(dim, keepdim=True) + eps
    return masked_exps / masked_sums
