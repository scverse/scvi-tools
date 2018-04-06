import numpy as np

from . import GeneExpressionDataset


class SyntheticDataset(GeneExpressionDataset):
    def __init__(self, batch_size=20, nb_genes=100):
        # Generating samples according to a ZINB process
        data = np.random.negative_binomial(5, 0.3, size=(batch_size, nb_genes))
        newdata = np.ones((batch_size, nb_genes))
        mask = np.random.binomial(n=1, p=0.7, size=(batch_size, nb_genes))
        for i in range(batch_size):
            newdata[i, :] = data[i, :] / np.sum(data[i, :])
            newdata[i, :] = newdata[i, :] * mask[i, :]
        super(SyntheticDataset, self).__init__([newdata])
