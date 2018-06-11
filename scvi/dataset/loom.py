from .dataset import GeneExpressionDataset
import loompy
import os
import numpy as np
import scipy.sparse as sp_sparse


class LoomDataset(GeneExpressionDataset):

    def __init__(self, filepath, unit_test=False):
        path, file = os.path.split(filepath)
        self.save_path = path + '/' if not unit_test else 'tests/data/'
        self.download_name = file if not unit_test else 'retina.loom'

        data, labels, cell_batches = self.preprocess()

        # Currently only supports retina.loom dataset
        data_with_info = np.hstack((data, labels, cell_batches))
        first_batch = data_with_info[data_with_info[:, -1] == 0.0]
        second_batch = data_with_info[data_with_info[:, -1] == 1.0]
        first_batch = first_batch[:, :-1]
        second_batch = second_batch[:, :-1]
        print("Finished preprocessing dataset")
        super(LoomDataset, self).__init__(*GeneExpressionDataset.get_attributes_from_list(
            [sp_sparse.csr_matrix(first_batch[:, :-1]), sp_sparse.csr_matrix(second_batch[:, :-1])],
            list_labels=[first_batch[:, -1], second_batch[:, -1]]))

    def preprocess(self):
        ds = loompy.connect(self.save_path + self.download_name)
        cell_batches = ds.ca['Batch_id']
        cell_batches = np.reshape(cell_batches, (cell_batches.shape[0], 1))

        labels = np.array(ds.ca['Clusters'])
        labels = np.reshape(labels, (labels.shape[0], 1))

        data = np.array(ds[:, :]).T  # change matrix to cells by genes
        ds.close()

        return data, labels, cell_batches
