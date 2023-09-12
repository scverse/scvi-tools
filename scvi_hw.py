# from scvi.module._vae import _compute_mmd, _compute_mmd_loss, _compute_fast_mmd
from scvi.model import SCVI
import scanpy as sc
import gdown
import torch


def download_dataset():
    url = 'https://drive.google.com/uc?id=1ehxgfHTsMZXy6YzlFKGJOsBKQ5rrvMnd'
    output = 'pancreas.h5ad'
    gdown.download(url, output, quiet=False)


def read_dataset():
    adata_all = sc.read("pancreas.h5ad")
    adata = adata_all.raw.to_adata()
    return adata


def preprocess_ann_data(adata):
    print("# cells, # genes before filtering:", adata.shape)

    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.filter_cells(adata, min_counts=3)

    print("# cells, # genes after filtering:", adata.shape)

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.raw = adata
    return adata


def main():
    # download_dataset()
    adata = read_dataset()
    adata = preprocess_ann_data(adata)
    SCVI.setup_anndata(adata, layer="counts", batch_key="study")
    model = SCVI(adata)
    model.view_anndata_setup(adata)




if __name__ == '__main__':
    main()
    # ts = torch.tensor([1,2,3,4])
    # ts.to('cuda')
    # print(ts)


#   TBD (to be done):
#   1.  Figure out what's wrong with the _compute_fast_mmd function
#   2.  Figure out how to implement the SCVI loss function?
#   3.  How to run the model?
#   4.  How to get the results?
#   5.  How to plot the results? the links are broken
#   6.  What is the smoke test? How to run it?
#   7.  How to use the colab notebook?

