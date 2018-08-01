import numpy as np
import loompy
from .dataset import GeneExpressionDataset


class SmfishDataset(GeneExpressionDataset):
    def __init__(self, save_path='data/'):

        self.download_name = 'osmFISH_SScortex_mouse_all_cell.loom'
        self.save_path = save_path
        self.url = 'http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom'

        data, labels, cell_types, x_coord, y_coord = self.download_and_preprocess()
        super(SmfishDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                data,
                labels=labels),
            x_coord=x_coord, y_coord=y_coord)

    def preprocess(self):
        print("Preprocessing smFISH_utils dataset")
        ds = loompy.connect(self.save_path + self.download_name)
        select = ds[:, :].sum(axis=0) > 0  # Take out cells that doesn't express any gene

        labels, cell_types = np.array(ds.ca['ClusterID']), np.array(ds.ca['ClusterName'])
        labels = np.reshape(labels, (labels.shape[0], 1))[select]
        cell_types = np.reshape(cell_types, (cell_types.shape[0], 1))[select]

        x_coord, y_coord = np.array(ds.ca['X']), np.array(ds.ca['Y'])
        x_coord = np.reshape(x_coord, (x_coord.shape[0], 1))[select]
        y_coord = np.reshape(y_coord, (y_coord.shape[0], 1))[select]

        data = ds[:, select].T  # change matrix to cells by genes
        ds.close()

        print("Finished preprocessing smFISH_utils dataset")
        return data, labels, cell_types, x_coord, y_coord
