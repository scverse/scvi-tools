import numpy as np
import scipy.sparse as sp_sparse

from . import GeneExpressionDataset


class SyntheticDataset(GeneExpressionDataset):
    @classmethod
    def get_attributes(self, batch_size=20, nb_genes=100, n_batches=2):
        # Generating samples according to a ZINB process
        data = np.random.negative_binomial(5, 0.3, size=(batch_size, nb_genes, n_batches))
        mask = np.random.binomial(n=1, p=0.7, size=(batch_size, nb_genes, n_batches))
        newdata = (data * mask).swapaxes(0, 2)  # We put the batch index first
        return [sp_sparse.csr_matrix(data) for data in newdata]

    @classmethod
    def get_dataset(cls, type='train'):
        return cls.from_list_batches(cls.get_attributes())
