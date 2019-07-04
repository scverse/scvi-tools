from . import GeneExpressionDataset
import numpy as np
import loompy


class DropseqDataset(GeneExpressionDataset):
    r""" Loads the dropseq dataset

    This is only the 71639 annotated cells located in the frontal cortex of adult mouses among the 690,000 cells
    studied by Saunders et al using the Drop-seq method. From the 29000 genes sequenced in the study, we only kept
    those whose variance in the expression matrix were above a certain threshold and we also kept all the genes who
    were studied in the Starmap experiment. We have a 71639*7611 gene expression matrix from which we subsample to
    keep only `subsample` cells
    """

    def __init__(self, save_path='data/', genes_to_keep=[], subsample=15000):
        self.genes_to_keep = genes_to_keep
        self.subsample = subsample

        self.url = ['https://github.com/YosefLab/scVI-data/raw/master/dropseq.loom']
        self.download_name = 'dropseq.loom'
        self.save_path = save_path

        data, batch_indices, labels, gene_names, cell_types, subclusters = self.preprocess()

        X, local_means, local_vars, batch_indices, labels, subclusters = \
            GeneExpressionDataset.get_attributes_from_matrix(data, batch_indices=batch_indices, labels=labels,
                                                             subclusters=subclusters)
        self.subclusters = subclusters
        super(DropseqDataset, self).__init__(X, local_means, local_vars, batch_indices, labels,
                                             gene_names=gene_names, cell_types=cell_types)

    def preprocess(self):
        print("Preprocessing dataset")
        gene_names, labels, batch_indices, cell_types = None, None, 0, None
        try:
            ds = loompy.connect(self.save_path + self.download_name)
            select = ds[:, :].sum(axis=0) > 0  # Take out cells that doesn't express any gene

            if 'Gene' in ds.ra:
                gene_names = [str(s).upper() for s in ds.ra['Gene']]

            if 'BatchID' in ds.ca:
                batch_indices = ds.ca['BatchID']
                batch_indices = np.reshape(batch_indices, (batch_indices.shape[0], 1))[select]

            if 'SubClusters' in ds.ca:
                subclusters = np.array(ds.ca['SubClusters'])
                subclusters = np.reshape(subclusters, (subclusters.shape[0], 1))[select]
            else:
                subclusters = None

            if 'Clusters' in ds.ca:
                labels = np.array(ds.ca['Clusters'])
                labels = np.reshape(labels, (labels.shape[0], 1))[select]

            if 'CellTypes' in ds.attrs:
                cell_types = np.array(ds.attrs['CellTypes'])

            data = ds[:, select].T  # change matrix to cells by genes
            ds.close()

            random_state = np.random.RandomState(0)
            cells = random_state.choice(np.arange(data.shape[0]), self.subsample, replace=False)
            data = data[cells]
            labels = labels[cells]
            if subclusters is not None:
                subclusters = subclusters[cells]
            if batch_indices != 0:
                batch_indices = batch_indices[cells]

            # reorder labels so that layers of the cortex appear in a right order when we plot cell types
            order_labels = [5, 6, 3, 2, 4, 0, 1, 8, 7, 9, 10, 11, 12, 13]

            labels = np.vectorize(lambda x: order_labels.index(x))(labels)
            cell_types = cell_types[order_labels]

            if len(self.genes_to_keep) > 0:
                gene_indices = []
                for g in self.genes_to_keep:
                    if g in gene_names:
                        gene_indices.append(gene_names.index(g))
                gene_indices = np.array(gene_indices)
                gene_names = np.array(gene_names)[gene_indices]
                data = data[:, gene_indices]

        except OSError:
            print("Error: the file " + self.download_name + " should be in the " + self.save_path + " directory.")
            raise

        print("Finished preprocessing dataset")
        return data, batch_indices, labels, gene_names, cell_types, subclusters
