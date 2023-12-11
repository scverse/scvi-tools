# from scvi.module._vae import _compute_mmd, _compute_mmd_loss, _compute_fast_mmd
import numpy as np
import scvi.module._vae
from scvi.model import SCVI
import scanpy as sc
import gdown
import matplotlib.pyplot as plt
import pandas as pd
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



def umap(pancreas_ref, scvi_ref):

    SCVI_LATENT_KEY = "X_scVI"
    pancreas_ref.obsm[SCVI_LATENT_KEY] = scvi_ref.get_latent_representation()
    sc.pp.neighbors(pancreas_ref, use_rep=SCVI_LATENT_KEY)
    sc.tl.leiden(pancreas_ref)
    sc.tl.umap(pancreas_ref)
    sc.pl.umap(
        pancreas_ref,
        color=["study", "cell_type"],
        frameon=False,
        ncols=1,
    )

def main():
    # download_dataset()
    adata = read_dataset()
    adata = preprocess_ann_data(adata)
    SCVI.setup_anndata(adata, layer="counts", batch_key="study")
    model = SCVI(adata)
    model.view_anndata_setup(adata)
    model.train(max_epochs=10)
    umap(adata, model)

if __name__ == '__main__':
    main()



