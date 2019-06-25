from . import GeneExpressionDataset
import numpy as np
import loompy


class StarmapDataset(GeneExpressionDataset):
    r""" Loads the starmap dataset


    """

    def __init__(self, save_path='data/', file='starmap.loom', without_positions=False, scrna_seq_genes=None):

        self.save_path = save_path
        self.url = ['https://github.com/YosefLab/scVI-data/raw/master/starmap.loom']
        self.download_name = file
        self.seq_genes = scrna_seq_genes

        data, batch_indices, labels, gene_names, cell_types, x_coord, y_coord = self.preprocess()
        if without_positions:
            x_coord = None
            y_coord = None

        X, local_means, local_vars, batch_indices, labels = \
            GeneExpressionDataset.get_attributes_from_matrix(data, batch_indices=batch_indices, labels=labels)
        super(StarmapDataset, self).__init__(X, local_means, local_vars, batch_indices, labels,
                                             gene_names=gene_names, cell_types=cell_types, x_coord=x_coord,
                                             y_coord=y_coord)

    def preprocess(self):
        print("Preprocessing dataset")
        gene_names, labels, batch_indices, cell_types = None, None, 0, None
        try:
            ds = loompy.connect(self.save_path + self.download_name)
            select = ds[:, :].sum(axis=0) > 0  # Take out cells that doesn't express any gene
            x_coord = None
            y_coord = None

            if 'Gene' in ds.ra:
                gene_names = ds.ra['Gene']

            if 'BatchID' in ds.ca:
                batch_indices = ds.ca['BatchID']
                batch_indices = np.reshape(batch_indices, (batch_indices.shape[0], 1))[select]

            if 'Clusters' in ds.ca:
                labels = np.array(ds.ca['Clusters'])
                labels = np.reshape(labels, (labels.shape[0], 1))[select]

            if 'Spatial_coordinates' in ds.ca:
                positions = np.array(ds.ca['Spatial_coordinates'])
                x_coord = positions[select, 0]
                y_coord = positions[select, 1]
                x_coord = np.reshape(x_coord, (x_coord.shape[0], 1))
                y_coord = np.reshape(y_coord, (y_coord.shape[0], 1))

            if 'CellTypes' in ds.attrs:
                cell_types = np.array(ds.attrs['CellTypes'])

            indices = np.arange(len(gene_names))
            if self.seq_genes is not None:
                gene_names, indices, _ = np.intersect1d(gene_names, self.seq_genes, return_indices=True)

            data = np.array(ds[:, select].T)  # change matrix to cells by genes
            data = data[:, indices]
            ds.close()
        except OSError:
            print("Error: the file " + self.download_name + " should be in the " + self.save_path + " directory.")
            raise

        print("Finished preprocessing dataset")
        return data, batch_indices, labels, gene_names, cell_types, x_coord, y_coord
