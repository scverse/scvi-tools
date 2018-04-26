import torch


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


def to_cuda(tensor_list, async=True):
    return [t.cuda(async=async) for t in tensor_list]
