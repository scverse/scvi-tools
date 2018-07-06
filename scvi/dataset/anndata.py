from .dataset import GeneExpressionDataset
import anndata
import numpy as np


class AnnDataset(GeneExpressionDataset):

    def __init__(self, download_name, save_path='data/', url=None, new_n_genes=False, subset_genes=None):
        self.download_name = download_name
        self.save_path = save_path
        self.url = url

        data, gene_names = self.download_and_preprocess()

        super(AnnDataset, self).__init__(*GeneExpressionDataset.get_attributes_from_matrix(data),
                                         gene_names=gene_names)

        self.subsample_genes(new_n_genes=new_n_genes, subset_genes=subset_genes)

    def preprocess(self):
        print("Preprocessing gene_dataset")

        ad = anndata.read_h5ad(self.save_path + self.download_name)  # obs = cells, var = genes
        gene_names = np.array(ad.var.index.values, dtype=str)
        data = ad.X
        select = data.sum(axis=1) > 0  # Take out cells that doesn't express any gene
        data = data[select, :]

        print("Finished preprocessing gene_dataset")
        return data, gene_names
