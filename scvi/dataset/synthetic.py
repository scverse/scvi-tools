import numpy as np
import scipy.sparse as sp_sparse

from . import GeneExpressionDataset


class SyntheticDataset(GeneExpressionDataset):
    def __init__(self, batch_size=20, nb_genes=100, n_batches=2):
        # Generating samples according to a ZINB process
        data = np.random.negative_binomial(5, 0.3, size=(batch_size, nb_genes, n_batches))
        mask = np.random.binomial(n=1, p=0.7, size=(batch_size, nb_genes, n_batches))
        newdata = (data * mask).swapaxes(0, 2)  # We put the batch index first
        super(SyntheticDataset, self).__init__([sp_sparse.csr_matrix(data) for data in newdata])
