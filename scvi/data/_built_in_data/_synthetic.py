import logging
from typing import Optional

import numpy as np
import pandas as pd
import scipy
from anndata import AnnData
from mudata import MuData

from scvi._types import AnnOrMuData

logger = logging.getLogger(__name__)


def _generate_synthetic(
    *,
    batch_size: int,
    n_genes: int,
    n_proteins: int,
    n_regions: int,
    n_batches: int,
    n_labels: int,
    dropout_ratio: float,
    sparse_format: Optional[str],
    return_mudata: bool,
    batch_key: str = "batch",
    labels_key: str = "labels",
    rna_key: str = "rna",
    protein_expression_key: str = "protein_expression",
    protein_names_key: str = "protein_names",
    accessibility_key: str = "accessibility",
) -> AnnOrMuData:
    n_obs = batch_size * n_batches

    def sparsify_data(data: np.ndarray):
        if sparse_format is not None:
            data = getattr(scipy.sparse, sparse_format)(data)
        return data

    rna = np.random.negative_binomial(5, 0.3, size=(n_obs, n_genes))
    mask = np.random.binomial(n=1, p=dropout_ratio, size=(n_obs, n_genes))
    rna = rna * mask
    rna = sparsify_data(rna)

    if n_proteins > 0:
        protein = np.random.negative_binomial(5, 0.3, size=(n_obs, n_proteins))
        protein_names = np.arange(n_proteins).astype(str)
        protein = sparsify_data(protein)

    if n_regions > 0:
        accessibility = np.random.negative_binomial(5, 0.3, size=(n_obs, n_regions))
        mask = np.random.binomial(n=1, p=dropout_ratio, size=(n_obs, n_regions))
        accessibility = accessibility * mask
        accessibility = sparsify_data(accessibility)

    batch = []
    for i in range(n_batches):
        batch += [f"batch_{i}"] * batch_size

    if n_labels > 0:
        labels = np.random.randint(0, n_labels, size=(n_obs,))
        labels = np.array([f"label_{i}" for i in labels])

    adata = AnnData(rna)
    if return_mudata:
        mod_dict = {rna_key: adata}

        if n_proteins > 0:
            protein_adata = AnnData(protein)
            protein_adata.var_names = protein_names
            mod_dict[protein_expression_key] = protein_adata
        if n_regions > 0:
            mod_dict[accessibility_key] = AnnData(accessibility)

        adata = MuData(mod_dict)
    else:
        if n_proteins > 0:
            adata.obsm[protein_expression_key] = protein
            adata.uns[protein_names_key] = protein_names
        if n_regions > 0:
            adata.obsm[accessibility_key] = accessibility

    adata.obs[batch_key] = pd.Categorical(batch)
    if n_labels > 0:
        adata.obs[labels_key] = pd.Categorical(labels)

    return adata
