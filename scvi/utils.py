import torch
from torch.autograd import Variable


def one_hot(index, n_cat, dtype):
    onehot = (index.new(index.size(0), n_cat).fill_(0))
    onehot.scatter_(1, index, 1)
    return Variable(onehot.type(dtype))


def enumerate_discrete(x, y_dim):
    def batch(batch_size, label):
        labels = (torch.ones(batch_size, 1) * label).type(torch.LongTensor)
        return one_hot(labels, y_dim, x.data.type())

    batch_size = x.size(0)
    return torch.cat([batch(batch_size, i) for i in range(y_dim)])
