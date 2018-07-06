# On the 10X website:
# The main categories (Cell Ranger 1.1.0 / Cell Ranger 2.1.0 / ...)
# have same access suffix for each of their dataset.
# For dataset name (eg. 'pbmc8k', 'pbmc4k', ect...) their are two available specifications,
# either filtered or raw data
import os
import tarfile
import pandas as pd
import numpy as np
from scipy import io


from scvi.dataset import GeneExpressionDataset

available_datasets = {"1.1.0":
                      ["frozen_pbmc_donor_a",
                       "frozen_pbmc_donor_b",
                       "frozen_pbmc_donor_c"],
                      "2.1.0":
                      ["pbmc8k",
                       "pbmc4k",
                       "t_3k",
                       "t_4k",
                       "neuron_9k"]}

to_groups = dict([(dataset_name, group)
                  for group, list_datasets in available_datasets.items()
                  for dataset_name in list_datasets])
available_specification = ['filtered', 'raw']


class Dataset10X(GeneExpressionDataset):
    r""" Loads a file from `10x`_ website.

    Args:
        :filename: Name of the dataset file.
        :save_path: Save path of the dataset. Default: ``'data/'``.
        :type: Either `filtered` data or `raw` data. Default: ``'filtered'``.
        :new_n_genes: Number of subsampled genes. Default: ``3000``.
        :subset_genes: List of genes for subsampling. Default: ``None``.

    Examples:
        >>> tenX_dataset = Dataset10X("neuron_9k")

    .. _10x:
        http://cf.10xgenomics.com/

    """
    def __init__(self, filename, save_path='data/', type='filtered', new_n_genes=3000, subset_genes=None):
        group = to_groups[filename]
        self.url = ("http://cf.10xgenomics.com/samples/cell-exp/%s/%s/%s_%s_gene_bc_matrices.tar.gz" %
                    (group, filename, filename, type))
        self.save_path = save_path + '10X/%s/' % filename
        self.save_name = '%s_gene_bc_matrices' % type

        self.download_name = self.save_name + '.tar.gz'
        expression_data, gene_names = self.download_and_preprocess()
        super(Dataset10X, self).__init__(*GeneExpressionDataset.get_attributes_from_matrix(
            expression_data), gene_names=gene_names)

        if new_n_genes is not None:
            self.subsample_genes(new_n_genes, subset_genes)

    def preprocess(self):
        print("Preprocessing dataset")
        if len(os.listdir(self.save_path)) == 1:  # nothing extracted yet
            print("Extracting tar file")
            tar = tarfile.open(self.save_path + self.download_name, "r:gz")
            tar.extractall(path=self.save_path)
            tar.close()

        path = self.save_path + self.save_name + '/'
        path += os.listdir(path)[0] + '/'
        genes_info = pd.read_csv(path + 'genes.tsv', sep='\t', header=None)
        gene_names = genes_info.values[:, 0].astype(np.str).ravel()
        expression_data = io.mmread(path + 'matrix.mtx').T.A

        print("Finished preprocessing dataset")
        return expression_data, gene_names


class BrainSmallDataset(Dataset10X):
    r"""
    This dataset consists in 9,128 mouse brain cells profiled using `10x`_ is used as a complement of PBMC for our
    study of zero abundance and quality control metrics correlation with our generative posterior parameters. We
    derived quality control metrics using the cellrangerRkit R package (v.1.1.0). Quality metrics were extracted from
    CellRanger throughout the molecule specific information file. We kept the top 3000 genes by variance. We used the
    clusters provided by cellRanger for the correlation analysis of zero probabilities.

    Args:
        :save_path: Save path of raw data file. Default: ``'data/'``.

    Examples:
        >>> gene_dataset = BrainSmallDataset()

    .. _10x:
        https://support.10xgenomics.com/single-cell-gene-expression/datasets

    """
    def __init__(self, save_path='data/'):
        super(BrainSmallDataset, self).__init__(filename="neuron_9k",
                                                save_path=save_path)
