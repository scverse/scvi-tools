from .dataset import GeneExpressionDataset
import loompy
import numpy as np
import scipy.sparse as sp_sparse

gene_row = 'Gene'
batch_col = 'BatchID'
cluster_col = 'ClusterID'


class LoomDataset(GeneExpressionDataset):

    def __init__(self, filename, save_path, url=None):
        self.download_name = filename
        self.save_path = save_path
        self.url = url

        self.has_gene, self.has_batch, self.has_cluster = False, False, False

        data, gene_names, labels, cell_batches = self.download_and_preprocess()
        print("Finished preprocessing dataset")

        if self.has_cluster and self.has_batch:  # File has cluster and batch info
            data_with_info = np.hstack((data, labels, cell_batches))
            first_batch = data_with_info[data_with_info[:, -1] == 0.0]
            second_batch = data_with_info[data_with_info[:, -1] == 1.0]
            first_batch = first_batch[:, :-1]
            second_batch = second_batch[:, :-1]
            super(LoomDataset, self).__init__(*GeneExpressionDataset.get_attributes_from_list(
                [sp_sparse.csr_matrix(first_batch[:, :-1]), sp_sparse.csr_matrix(second_batch[:, :-1])],
                list_labels=[first_batch[:, -1], second_batch[:, -1]]))
        elif self.has_gene:  # File has gene names info
            super(LoomDataset, self).__init__(
                *GeneExpressionDataset.get_attributes_from_matrix(
                    data, labels=labels), gene_names=gene_names)

    def preprocess(self):
        gene_names, labels, cell_batches = None, None, None
        ds = loompy.connect(self.save_path + self.download_name)
        select = ds[:, :].sum(axis=0) > 0  # Take out cells that doesn't express any gene

        self.has_gene = gene_row in ds.ra.keys()
        self.has_batch = batch_col in ds.ca.keys()
        self.has_cluster = cluster_col in ds.ca.keys()

        if self.has_gene:
            gene_names = ds.ra[gene_row]

        if self.has_batch:
            cell_batches = ds.ca[batch_col]
            cell_batches = np.reshape(cell_batches, (cell_batches.shape[0], 1))[select]

        if self.has_cluster:
            labels = np.array(ds.ca[cluster_col])
            labels = np.reshape(labels, (labels.shape[0], 1))[select]

        data = ds[:, select].T  # change matrix to cells by genes
        ds.close()

        return data, gene_names, labels, cell_batches
