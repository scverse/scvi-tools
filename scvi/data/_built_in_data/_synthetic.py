import logging

import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData
from scipy.sparse import csr_matrix

from scvi._types import AnnOrMuData

logger = logging.getLogger(__name__)


def _generate_synthetic(
    *,
    batch_size,
    n_genes,
    n_proteins,
    n_regions,
    n_batches,
    n_labels,
    sparse,
    dropout_ratio,
    return_mudata,
    batch_key: str = "batch",
    labels_key: str = "labels",
    rna_key: str = "rna",
    protein_expression_key: str = "protein_expression",
    protein_names_key: str = "protein_names",
    accessibility_key: str = "accessibility",
) -> AnnOrMuData:
    n_obs = batch_size * n_batches

    rna = np.random.negative_binomial(5, 0.3, size=(n_obs, n_genes))
    mask = np.random.binomial(n=1, p=dropout_ratio, size=(n_obs, n_genes))
    rna = rna * mask
    rna = csr_matrix(rna) if sparse else rna

    if n_proteins > 0:
        protein = np.random.negative_binomial(5, 0.3, size=(n_obs, n_proteins))
        protein_names = np.arange(n_proteins).astype(str)

    if n_regions > 0:
        accessibility = np.random.negative_binomial(5, 0.3, size=(n_obs, n_regions))
        mask = np.random.binomial(n=1, p=dropout_ratio, size=(n_obs, n_regions))
        accessibility = accessibility * mask
        accessibility = csr_matrix(accessibility) if sparse else accessibility

    if n_labels > 0:
        labels = np.random.choice(
            [f"label_{i}" for i in range(n_labels)], size=(n_obs,)
        )

    batch = np.random.choice([f"batch_{i}" for i in range(n_batches)], size=(n_obs,))

    if return_mudata:
        mod_dict = {rna_key: AnnData(rna, dtype=np.float32)}

        if n_proteins > 0:
            protein_adata = AnnData(protein, dtype=np.float32)
            protein_adata.var_names = protein_names
            mod_dict[protein_expression_key] = protein_adata
        if n_regions > 0:
            mod_dict[accessibility_key] = AnnData(accessibility, dtype=np.float32)

        adata = MuData(mod_dict)
    else:
        adata = AnnData(rna, dtype=np.float32)
        if n_proteins > 0:
            adata.obsm[protein_expression_key] = protein
            adata.uns[protein_names_key] = protein_names
        if n_regions > 0:
            adata.obsm[accessibility_key] = accessibility

    adata.obs[batch_key] = pd.Categorical(batch)
    adata.obs[labels_key] = pd.Categorical(labels)

    return adata
