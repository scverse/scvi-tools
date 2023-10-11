import logging

import numpy as np
import pandas as pd
import scipy
from anndata import AnnData
from mudata import MuData

from scvi._types import AnnOrMuData

logger = logging.getLogger(__name__)


def _sparsify_array(array: np.ndarray, sparse_format: str | None) -> np.ndarray:
    if sparse_format is not None:
        array = getattr(scipy.sparse, sparse_format)(array)
    return array


def _generate_zinb(
    n_obs: int,
    n_vars: int,
    n_nb: int = 5,
    p_nb: float = 0.3,
    p_zero: float = 0.7,
    sparse_format: str | None = None,
) -> np.ndarray:
    generator = np.random.default_rng()

    counts = generator.negative_binomial(n_nb, p_nb, size=(n_obs, n_vars))
    mask = generator.binomial(n=1, p=1 - p_zero, size=(n_obs, n_vars))
    counts = counts * mask

    return _sparsify_array(counts, sparse_format)


def _generate_categorical(
    n_obs: int,
    low: int,
    high: int,
    key: str | None,
    as_categorical: bool = True,
) -> np.ndarray:
    categorical = np.random.default_rng().integers(low, high, size=(n_obs,))
    if key is not None:
        categorical = np.array([f"{key}_{i}" for i in categorical])
    if as_categorical:
        categorical = pd.Categorical(categorical)
    return categorical


def _generate_continuous(n_obs: int) -> np.ndarray:
    return np.random.default_rng().normal(size=(n_obs,))


def _generate_synthetic(
    *,
    batch_size: int,
    n_genes: int,
    n_proteins: int,
    n_regions: int,
    n_batches: int,
    n_labels: int,
    n_categorical_covariates: list[int],
    n_continuous_covariates: int,
    sparse_format: str | None,
    return_mudata: bool,
    batch_key: str = "batch",
    labels_key: str = "labels",
    categorical_covariates_key: str = "categorical_covariate",
    continuous_covariates_key: str = "continuous_covariate",
    rna_key: str = "rna",
    protein_key: str = "protein_expression",
    protein_names_key: str = "protein_names",
    accessibility_key: str = "accessibility",
    rna_kwargs: dict | None = None,
    protein_kwargs: dict | None = None,
    accessibility_kwargs: dict | None = None,
) -> AnnOrMuData:
    rna_kwargs = rna_kwargs or {}
    protein_kwargs = protein_kwargs or {}
    accessibility_kwargs = accessibility_kwargs or {}

    n_obs = batch_size * n_batches
    gene_names = [f"gene_{i}" for i in range(n_genes)]
    protein_names = [f"protein_{i}" for i in range(n_proteins)]

    rna = _generate_zinb(n_obs, n_genes, sparse_format=sparse_format, **rna_kwargs)
    if n_proteins > 0:
        protein = _generate_zinb(
            n_obs, n_proteins, sparse_format=sparse_format, **protein_kwargs
        )
    if n_regions > 0:
        accessibility = _generate_zinb(
            n_obs, n_regions, sparse_format=sparse_format, **accessibility_kwargs
        )

    obs_columns = {}

    batch = []
    for i in range(n_batches):
        batch += [f"batch_{i}"] * batch_size
    obs_columns[batch_key] = pd.Categorical(batch)

    if n_labels > 0:
        labels = _generate_categorical(
            n_obs, 0, n_labels, "label", as_categorical=True
        )
        obs_columns[labels_key] = labels

    for i, n_categories in enumerate(n_categorical_covariates):
        key = f"{categorical_covariates_key}_{i}"
        obs_columns[key] = _generate_categorical(
            n_obs, 0, n_categories, key, as_categorical=False
        )

    for i in range(n_continuous_covariates):
        key = f"{continuous_covariates_key}_{i}"
        obs_columns[key] = _generate_continuous(n_obs)

    adata = AnnData(rna)
    adata.var_names = gene_names
    if return_mudata:
        mod_dict = {rna_key: adata}

        if n_proteins > 0:
            protein_adata = AnnData(protein)
            protein_adata.var_names = protein_names
            mod_dict[protein_key] = protein_adata
        if n_regions > 0:
            mod_dict[accessibility_key] = AnnData(accessibility)

        adata = MuData(mod_dict)
    else:
        if n_proteins > 0:
            adata.obsm[protein_key] = protein
            adata.uns[protein_names_key] = protein_names
        if n_regions > 0:
            adata.obsm[accessibility_key] = accessibility

    for key, value in obs_columns.items():
        adata.obs[key] = value

    return adata
