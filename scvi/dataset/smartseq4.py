from . import GeneExpressionDataset
import numpy as np
import loompy
import os


class Smartseq4Dataset(GeneExpressionDataset):
    r""" Loads the smartseq V4 dataset
    """

    def __init__(self, save_path='data/', gene_normalized=True):

        self.save_path = save_path
        self.url = ['https://github.com/YosefLab/scVI-data/raw/master/smartseq4.loom']
        self.download_name = "smartseq4 .loom"

        data, labels, gene_names, cell_types, complete_clusters, complete_classes, complete_labels = \
            self.preprocess()
        self.complete_clusters = complete_clusters
        self.complete_classes = complete_classes
        self.complete_labels = complete_labels
        if gene_normalized:
            None

        X, local_means, local_vars, batch_indices, labels = \
            GeneExpressionDataset.get_attributes_from_matrix(data, labels=labels)
        super(Smartseq4Dataset, self).__init__(X, local_means, local_vars, batch_indices, labels,
                                               gene_names=gene_names, cell_types=cell_types)

    def preprocess(self):
        print("Preprocessing dataset")
        try:
            ds = loompy.connect(os.path.join(self.save_path, self.download_name))
            select = ds[:, :].sum(axis=0) > 0  # Take out cells that doesn't express any gene
            gene_names = list(ds.ra['Gene'])
            labels = np.array(ds.ca['Clusters'])
            labels = np.reshape(labels, (labels.shape[0], 1))[select]
            cell_types = np.array(ds.attrs['CellTypes'])

            # additional metadata attributes that we might use later
            complete_clusters = ds.ca["precise_clusters"]
            complete_labels = ds.ca["precise_labels"]
            complete_classes = ds.ca["precise_classes"]

            data = np.array(ds[:, select].T)  # change matrix to cells by genes
            fish_genes = np.load('data/Fish_genes.npy')
            fish_genes = np.array([gene_names.index(gene) for gene in fish_genes if gene in gene_names])
            selected_genes = np.std(data, axis=0).argsort()[-7000:]
            selected_genes = np.unique(np.concatenate((selected_genes, fish_genes)))
            gene_names = np.array(gene_names)[selected_genes]
            data = data[:, selected_genes]
            ds.close()
        except OSError:
            print("Error: the file " + self.download_name + " should be in the " + self.save_path + " directory.")
            raise

        print("Finished preprocessing dataset")
        return data, labels, gene_names, cell_types, complete_clusters, complete_classes, complete_labels
