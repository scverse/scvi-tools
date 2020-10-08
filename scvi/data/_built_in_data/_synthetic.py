import logging

import numpy as np
import pandas as pd
from anndata import AnnData

from scvi.data import setup_anndata

logger = logging.getLogger(__name__)


def _generate_synthetic(
    batch_size: int = 128,
    n_genes: int = 100,
    n_proteins: int = 100,
    n_batches: int = 2,
    n_labels: int = 3,
    run_setup_anndata: bool = True,
) -> AnnData:

    data = np.random.negative_binomial(5, 0.3, size=(batch_size * n_batches, n_genes))
    mask = np.random.binomial(n=1, p=0.7, size=(batch_size * n_batches, n_genes))
    data = data * mask  # We put the batch index first
    labels = np.random.randint(0, n_labels, size=(batch_size * n_batches,))
    labels = np.array(["label_%d" % i for i in labels])

    batch = []
    for i in range(n_batches):
        batch += ["batch_{}".format(i)] * batch_size

    adata = AnnData(data)
    adata.obs["batch"] = pd.Categorical(batch)
    adata.obs["labels"] = pd.Categorical(labels)

    # Protein measurements
    p_data = np.random.negative_binomial(5, 0.3, size=(adata.shape[0], n_proteins))
    adata.obsm["protein_expression"] = p_data
    adata.uns["protein_names"] = np.arange(n_proteins).astype(str)
    if run_setup_anndata:
        setup_anndata(
            adata,
            batch_key="batch",
            labels_key="labels",
            protein_expression_obsm_key="protein_expression",
            protein_names_uns_key="protein_names",
        )

    return adata
