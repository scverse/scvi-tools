from . import GeneExpressionDataset
import numpy as np
import loompy


class StarmapDataset(GeneExpressionDataset):
    r""" Loads the dropseq dataset

    This is only 70 000 annotated cells located in the cortex of adult mouses among the 700 000 cells studied by
    Saunders et al using the Drop-seq method. From the 29000 genes sequenced in the study, we only kept those whose
    variance in the expression matrix were above a certain threshold and we also kept all the genes who were studied in
    the Starmap experiment. We work with a 70000*6000 gene expression matrix.
    """

    def __init__(self, save_path='data/', genes_starmap=[], without_positions=False):
        # If we want to harmonize the dataset with the StarMAP dataset, we need to keep
        # StarMAP genes and order the genes from Dropseq accordingly
        self.genes_starmap = genes_starmap

        self.save_path = save_path
        self.url = ['https://github.com/YosefLab/scVI-data/raw/master/starmap.loom']
        self.download_name = 'dropseq.loom'

        data, batch_indices, labels, gene_names, cell_types, x_coord, y_coord = self.download_and_preprocess()
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

        if 'CellTypes' in ds.attrs:
            cell_types = np.array(ds.attrs['CellTypes'])

        data = ds[:, select].T  # change matrix to cells by genes
        ds.close()

        print("Finished preprocessing dataset")
        return data, batch_indices, labels, gene_names, cell_types, x_coord, y_coord
