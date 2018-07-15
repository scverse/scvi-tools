import loompy
import numpy as np

from .dataset import GeneExpressionDataset


class LoomDataset(GeneExpressionDataset):
    def __init__(self, filename, save_path='data/', url=None):
        self.download_name = filename
        self.save_path = save_path
        self.url = url

        self.has_gene, self.has_batch, self.has_cluster = False, False, False

        data, batch_indices, labels, gene_names, cell_types = self.download_and_preprocess()

        X, local_means, local_vars, batch_indices_, labels = \
            GeneExpressionDataset.get_attributes_from_matrix(data, labels=labels)
        batch_indices = batch_indices if batch_indices is not None else batch_indices_
        super(LoomDataset, self).__init__(X, local_means, local_vars, batch_indices, labels,
                                          gene_names=gene_names, cell_types=cell_types)

    def preprocess(self):
        print("Preprocessing dataset")
        gene_names, labels, batch_indices, cell_types = None, None, None, None
        ds = loompy.connect(self.save_path + self.download_name)
        select = ds[:, :].sum(axis=0) > 0  # Take out cells that doesn't express any gene

        if 'Gene' in ds.ra:
            gene_names = ds.ra['Gene']

        if 'BatchID' in ds.ca:
            batch_indices = ds.ca['BatchID']
            batch_indices = np.reshape(batch_indices, (batch_indices.shape[0], 1))[select]

        if 'ClusterID' in ds.ca:
            labels = np.array(ds.ca['ClusterID'])
            labels = np.reshape(labels, (labels.shape[0], 1))[select]

        if 'CellTypes' in ds.attrs:
            cell_types = np.array(ds.attrs['CellTypes'])

        data = ds[:, select].T  # change matrix to cells by genes
        ds.close()

        print("Finished preprocessing dataset")
        return data, batch_indices, labels, gene_names, cell_types
