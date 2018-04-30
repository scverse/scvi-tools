import torch
from torch.distributions import register_kl, Multinomial


def one_hot(index, n_cat):
    onehot = torch.zeros(index.size(0), n_cat, device=index.device)
    onehot.scatter_(1, index, 1)
    return onehot.type(torch.float32)


def enumerate_discrete(x, y_dim):
    def batch(batch_size, label):
        labels = (torch.ones(batch_size, 1, device=x.device, dtype=torch.long) * label)
        return one_hot(labels, y_dim)

    batch_size = x.size(0)
    return torch.cat([batch(batch_size, i) for i in range(y_dim)])


@register_kl(Multinomial, Multinomial)
def kl_multinomial_multinomial(p, q):
    return torch.sum(torch.mul(p.probs, torch.log(q.probs) - torch.log(p.probs + 1e-8)), dim=-1)
