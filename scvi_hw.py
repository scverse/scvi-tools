# from scvi.module._vae import _compute_mmd, _compute_mmd_loss, _compute_fast_mmd
import numpy as np
import scvi.module._vae
from scvi.model import SCVI
import scanpy as sc
import gdown
import matplotlib.pyplot as plt
import torch
import os
import tempfile
import warnings


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


def error_function_printout(epoch_num,error):
    epoch = np.arange(1,epoch_num)
    loss_epoch = error
    plt.plot(epoch, loss_epoch)
    plt.xlabel('Epoch')
    plt.ylabel('Error')
    plt.title('Error as function of the Epoch')
    plt.grid(True)
    plt.show()



def UMAP():
    pass
def main():
    # download_dataset()
    adata = read_dataset()
    adata = preprocess_ann_data(adata)
    SCVI.setup_anndata(adata, layer="counts", batch_key="study")
    model = SCVI(adata)
    model.view_anndata_setup(adata)
    print(model.history_)
    model.train(max_epochs=1)
    print(model.history_)
    # sc.set_figure_params(figsize=(4, 4), frameon=False)
    # torch.set_float32_matmul_precision("high")
    # save_dir = tempfile.TemporaryDirectory()
    # warnings.simplefilter(action="ignore", category=FutureWarning)
    # pancreas_adata_path = os.path.join(save_dir.name, "pancreas.h5ad")
    # pancreas_adata = sc.read(
    #     pancreas_adata_path,
    #     backup_url="https://figshare.com/ndownloader/files/24539828",
    # )
    # pancreas_adata
    # query_mask = np.array(
    #     [s in ["smartseq2", "celseq2"] for s in pancreas_adata.obs["tech"]]
    # )
    # # pancreas_ref = pancreas_adata[~query_mask].copy()

    # pancreas_query = pancreas_adata[query_mask].copy()
    # sc.pp.highly_variable_genes(
    #     pancreas_ref, n_top_genes=2000, batch_key="tech", subset=True
    # )
    # pancreas_query = pancreas_query[:, pancreas_ref.var_names].copy()
    # SCVI_LATENT_KEY = "X_scVI"
    # pancreas_ref.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()
    # sc.pp.neighbors(pancreas_ref, use_rep=SCVI_LATENT_KEY)
    # sc.tl.leiden(pancreas_ref)
    # sc.tl.umap(pancreas_ref)
    # sc.pl.umap(
    #     pancreas_ref,
    #     color=["tech", "celltype"],
    #     frameon=False,
    #     ncols=1,
    # )


if __name__ == '__main__':

    main()

    # z1 = torch.tensor([1,2,3])




    # for i in range(2):
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

