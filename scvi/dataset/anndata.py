from .dataset import GeneExpressionDataset
import anndata
import numpy as np
import os


class AnnDataset(GeneExpressionDataset):
    r"""Loads a `.h5ad` file .

    ``AnnDataset`` class supports loading `Anndata`_ object.

    Args:
        :filename: Name of the `.h5ad` file.
        :save_path: Save path of the dataset. Default: ``'data/'``.
        :url: Url of the remote dataset. Default: ``None``.
        :new_n_genes: Number of subsampled genes. Default: ``False``.
        :subset_genes: List of genes for subsampling. Default: ``None``.


    Examples:
        >>> # Loading a local dataset
        >>> local_ann_dataset = AnnDataset("TM_droplet_mat.h5ad", save_path = 'data/')

    .. _Anndata:
        http://anndata.readthedocs.io/en/latest/

    """

    def __init__(
        self,
        filename_or_anndata,
        save_path="data/",
        url=None,
        new_n_genes=False,
        subset_genes=None,
    ):
        """


        """

        if type(filename_or_anndata) == str:
            self.download_name = filename_or_anndata
            self.save_path = save_path
            self.url = url

            data, gene_names = self.download_and_preprocess()

        elif isinstance(filename_or_anndata, anndata.AnnData):
            ad = filename_or_anndata
            data, gene_names = self.extract_data_from_anndata(ad)

        else:
            raise Exception(
                "Please provide a filename of an AnnData file or an already loaded AnnData object"
            )

        super().__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(data),
            gene_names=gene_names
        )

        self.subsample_genes(new_n_genes=new_n_genes, subset_genes=subset_genes)

    def preprocess(self):
        print("Preprocessing dataset")

        ad = anndata.read_h5ad(
            os.path.join(self.save_path, self.download_name)
        )  # obs = cells, var = genes
        data, gene_names = self.extract_data_from_anndata(ad)

        print("Finished preprocessing dataset")
        return data, gene_names

    def extract_data_from_anndata(self, ad: anndata.AnnData):
        self.obs = (
            ad.obs
        )  # provide access to observation annotations from the underlying AnnData object.
        gene_names = np.array(ad.var.index.values, dtype=str)
        if isinstance(ad.X, np.ndarray):
            data = ad.X.copy()  # Dense
        else:
            data = ad.X.toarray()  # Sparse
        select = data.sum(axis=1) > 0  # Take out cells that doesn't express any gene
        data = data[select, :]

        return data, gene_names
