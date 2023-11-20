import logging

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix

logger = logging.getLogger(__name__)


def _generate_synthetic(
    batch_size: int = 128,
    n_genes: int = 100,
    n_proteins: int = 100,
    n_batches: int = 2,
    n_labels: int = 3,
    sparse: bool = False,
    dropout_ratio: float = 0.7,
) -> AnnData:
    n_obs = batch_size * n_batches

    rna = np.random.negative_binomial(5, 0.3, size=(n_obs, n_genes))
    mask = np.random.binomial(n=1, p=dropout_ratio, size=(n_obs, n_genes))
    rna = rna * mask
    rna = csr_matrix(rna) if sparse else rna

    labels = np.random.randint(0, n_labels, size=(n_obs,))
    labels = np.array(["label_%d" % i for i in labels])

    batch = []
    for i in range(n_batches):
        batch += [f"batch_{i}"] * batch_size

    protein = np.random.negative_binomial(5, 0.3, size=(n_obs, n_proteins))
    protein_names = np.arange(n_proteins).astype(str)

    adata = AnnData(rna)
    adata.obs["batch"] = pd.Categorical(batch)
    adata.obs["labels"] = pd.Categorical(labels)
    adata.obsm["protein_expression"] = protein
    adata.uns["protein_names"] = protein_names

    return adata
