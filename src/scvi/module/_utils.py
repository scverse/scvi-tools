import torch


def broadcast_labels(o, n_broadcast=-1):
    """Utility for the semi-supervised setting.

    If y is defined(labelled batch) then one-hot encode the labels (no broadcasting needed)
    If y is undefined (unlabelled batch) then generate all possible labels (and broadcast other
    arguments if not None)
    """

    def batch(batch_size, label):
        labels = torch.ones(batch_size, 1, device=o.device, dtype=torch.long) * label
        # return size (batch_size, n_broadcast)
        return torch.nn.functional.one_hot(labels.squeeze(-1), n_broadcast)

    batch_size = o.size(-2)
    if o.ndim == 2:
        ys = torch.cat([batch(batch_size, i) for i in range(n_broadcast)])
        new_o = o.repeat(n_broadcast, 1)
    elif o.ndim == 3:
        n_samples = o.size(0)
        ys = torch.cat([batch(batch_size, i) for i in range(n_broadcast)])
        ys = ys.unsqueeze(0).repeat(n_samples, 1, 1)
        new_o = o.repeat(n_broadcast, 1, 1)

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
