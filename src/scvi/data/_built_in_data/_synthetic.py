from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import scipy
from anndata import AnnData
from mudata import MuData

if TYPE_CHECKING:
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
    sparse_format: str | None,
    generate_coordinates: bool,
    return_mudata: bool,
    batch_key: str = "batch",
    labels_key: str = "labels",
    rna_key: str = "rna",
    gene_names_prefix: str = "gene",
    protein_expression_key: str = "protein_expression",
    protein_names_key: str = "protein_names",
    protein_names_prefix: str = "protein",
    accessibility_key: str = "accessibility",
    region_names_prefix: str = "region",
    coordinates_key: str = "coordinates",
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
    gene_names = np.array([f"{gene_names_prefix}_{i}" for i in range(n_genes)])

    if n_proteins > 0:
        protein = np.random.negative_binomial(5, 0.3, size=(n_obs, n_proteins))
        protein = sparsify_data(protein)
        protein_names = np.array([f"{protein_names_prefix}_{i}" for i in range(n_proteins)])

    if n_regions > 0:
        accessibility = np.random.negative_binomial(5, 0.3, size=(n_obs, n_regions))
        mask = np.random.binomial(n=1, p=dropout_ratio, size=(n_obs, n_regions))
        accessibility = accessibility * mask
        accessibility = sparsify_data(accessibility)
        region_names = np.array([f"{region_names_prefix}_{i}" for i in range(n_regions)])

    batch = []
    for i in range(n_batches):
        batch += [f"batch_{i}"] * batch_size

    if n_labels > 0:
        labels = np.random.randint(0, n_labels, size=(n_obs,))
        labels = np.array([f"label_{i}" for i in labels])

    if generate_coordinates:
        coords = np.random.normal(size=(n_obs, 2))

    adata = AnnData(rna)
    adata.var_names = gene_names
    if return_mudata:
        mod_dict = {rna_key: adata}

        if n_proteins > 0:
            protein_adata = AnnData(protein)
            protein_adata.var_names = protein_names
            mod_dict[protein_expression_key] = protein_adata
        if n_regions > 0:
            accessibility_adata = AnnData(accessibility)
            accessibility_adata.var_names = region_names
            mod_dict[accessibility_key] = accessibility_adata

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
    if generate_coordinates:
        adata.obsm[coordinates_key] = coords

    return adata
