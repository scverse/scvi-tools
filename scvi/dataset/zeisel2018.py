from . import GeneExpressionDataset
import numpy as np
import loompy
import os


class Zeisel2018Dataset(GeneExpressionDataset):
    r""" Loads the mouse cortex dataset published by Zeisel in 2018


    """

    def __init__(self, save_path='data/'):

        self.save_path = save_path
        self.url = ['http://mousebrain.org/tissues.html/l1_cortex1.loom']
        self.download_name = "l1_cortex1.loom"

        data, labels, gene_names, batch_indices, cell_types, complete_clusters, global_classes = \
            self.preprocess()
        self.complete_clusters = complete_clusters
        self.global_classes = global_classes

        X, local_means, local_vars, batch_indices, labels = \
            GeneExpressionDataset.get_attributes_from_matrix(data, labels=labels, batch_indices=batch_indices)
        super(Zeisel2018Dataset, self).__init__(X, local_means, local_vars, batch_indices, labels,
                                                gene_names=gene_names, cell_types=cell_types)

    def preprocess(self):
        print("Preprocessing dataset")
        try:
            ds = loompy.connect(os.path.join(self.save_path, self.download_name))
            select = ds[:, :].sum(axis=0) > 0  # Take out cells that doesn't express any gene
            gene_names = list(ds.ra['Gene'])
            gene_names = [gene.lower() for gene in gene_names]
            labels = np.array(ds.ca['Subclass'])
            labels = np.reshape(labels, (labels.shape[0], 1))[select]
            cell_types = np.unique(labels)
            batch_indices = np.array(ds.ca["SampleID"])
            batch_indices = np.reshape(batch_indices, (batch_indices.shape[0], 1))[select]

            # additional metadata attributes that we might use later
            complete_clusters = np.array(ds.ca["Clusters"])[select]
            global_classes = np.array(ds.ca["Class"])[select]

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
        return data, labels, gene_names, batch_indices, cell_types, complete_clusters, global_classes
