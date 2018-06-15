from .dataset import GeneExpressionDataset
import loompy
import numpy as np
import scipy.sparse as sp_sparse


class LoomDataset(GeneExpressionDataset):

    def __init__(self, filename, save_path, gene_row=None, cluster_col=None, batch_col=None, url=None):
        self.download_name = filename
        self.save_path = save_path
        self.gene_row = gene_row
        self.batch_col = batch_col
        self.cluster_col = cluster_col
        self.url = url

        data, gene_names, labels, cell_batches = self.download_and_preprocess()
        print("Finished preprocessing dataset")

        if cluster_col is not None and batch_col is not None:  # File has cluster and batch info
            data_with_info = np.hstack((data, labels, cell_batches))
            first_batch = data_with_info[data_with_info[:, -1] == 0.0]
            second_batch = data_with_info[data_with_info[:, -1] == 1.0]
            first_batch = first_batch[:, :-1]
            second_batch = second_batch[:, :-1]
            super(LoomDataset, self).__init__(*GeneExpressionDataset.get_attributes_from_list(
                [sp_sparse.csr_matrix(first_batch[:, :-1]), sp_sparse.csr_matrix(second_batch[:, :-1])],
                list_labels=[first_batch[:, -1], second_batch[:, -1]]))
        elif gene_names is not None:  # File has gene names info
            super(LoomDataset, self).__init__(
                *GeneExpressionDataset.get_attributes_from_matrix(
                    data, labels=labels), gene_names=gene_names)

    def preprocess(self):
        gene_names, labels, cell_batches = None, None, None

        ds = loompy.connect(self.save_path + self.download_name)
        select = ds[:, :].sum(axis=0) > 0  # Take out cells that doesn't express any gene

        if self.gene_row:
            gene_names = ds.ra[self.gene_row]

        if self.batch_col:
            cell_batches = ds.ca[self.batch_col]
            cell_batches = np.reshape(cell_batches, (cell_batches.shape[0], 1))[select]

        if self.cluster_col:
            labels = np.array(ds.ca[self.cluster_col])
            labels = np.reshape(labels, (labels.shape[0], 1))[select]

        data = ds[:, select].T  # change matrix to cells by genes
        ds.close()

        return data, gene_names, labels, cell_batches
