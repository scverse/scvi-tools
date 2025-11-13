from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from anndata import AnnData
from sklearn.preprocessing import MinMaxScaler, OneHotEncoder, StandardScaler

from scvi import settings
from scvi.utils import dependencies


def validate_marker(adata: AnnData, marker: str | list[str]):
    if isinstance(marker, str):
        marker = [marker]
    for m in marker:
        if m not in adata.var_names:
            raise ValueError(f"Marker {m} not found in adata.var_names.")


def validate_obs_keys(adata: AnnData, obs_key: str | list[str]):
    if obs_key is not None:
        if isinstance(obs_key, str):
            obs_key = [obs_key]
        for key in obs_key:
            if key is not None:
                if key not in adata.obs:
                    raise ValueError(f"Key {key} not found in adata.obs.")


def validate_obsm_keys(adata, obsm_keys):
    if obsm_keys is not None:
        if isinstance(obsm_keys, str):
            obsm_keys = [obsm_keys]
        for key in obsm_keys:
            if key not in adata.obsm:
                raise KeyError(f"Key '{key}' not found in adata.obs or adata.obsm.")


def validate_layer_key(adata: AnnData, layer_key: str):
    if layer_key is not None:
        if layer_key not in adata.layers:
            raise ValueError(f"Layer key {layer_key} not found in adata.layers.")


def apply_scaling(data, method, feature_range):
    if method == "minmax":
        scaler = MinMaxScaler(feature_range=feature_range)
    elif method == "standard":
        scaler = StandardScaler()
    return scaler.fit_transform(data), scaler


def get_n_latent_heuristic(n_vars: int, latent_max: int = 20, latent_min: int = 10):
    n_latent = round(n_vars / 2)

    if n_latent > latent_max:
        n_latent = latent_max

    if n_latent < latent_min:
        n_latent = latent_min

    return n_latent


def clip_lfc_factory(min_lfc: float, max_lfc: float, pseudocount=0):
    def clip_lfc(x, y, pseudocount=pseudocount):
        x = np.clip(x, min_lfc, max_lfc)
        y = np.clip(y, min_lfc, max_lfc)
        return np.log2(x) - np.log2(y)

    return clip_lfc


def validate_expression_range(data, min_exp, max_exp):
    return np.all((data >= min_exp) & (data <= max_exp))


def encode_categories(adata, cat_key):
    """One-Hot encode the categories for the given key."""
    ohe = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
    return ohe.fit_transform(adata.obs[cat_key].values.reshape(-1, 1)), ohe


@dependencies("pynndescent")
def impute_cats_with_neighbors(
    rep_query, rep_ref, cat_encoded_ref, n_neighbors=5, compute_uncertainty=False
):
    """Use pynndescent to find nearest neighbors and impute missing categories."""
    import pynndescent

    nn_index = pynndescent.NNDescent(rep_ref, n_neighbors=n_neighbors, metric="euclidean")

    indices, distances = nn_index.query(rep_query, k=n_neighbors)
    neighbor_categories = cat_encoded_ref[indices]  # Shape: (n_query, n_neighbors, n_categories)

    category_sums = np.sum(neighbor_categories, axis=1)  # Shape: (n_query, n_categories)

    if compute_uncertainty:
        category_prop = category_sums / np.sum(category_sums, axis=1, keepdims=True)
        uncertainty = 1 - np.max(category_prop, axis=1)
    else:
        uncertainty = None

    imputed_cat_indices = np.argmax(category_sums, axis=1)  # Shape: (n_query,)

    return imputed_cat_indices, uncertainty


@dependencies("pynndescent")
def impute_expr_with_neighbors(
    rep_query, rep_ref, expr_data_ref, n_neighbors=5, compute_uncertainty=False
):
    """Use pynndescent to find nearest neighbors and impute missing expression."""
    import pynndescent

    nn_index = pynndescent.NNDescent(rep_ref, n_neighbors=n_neighbors, metric="euclidean")

    indices, distances = nn_index.query(rep_query, k=n_neighbors)
    neighbor_expr = expr_data_ref[indices]  # Shape: (n_query, n_neighbors, expr_dim)

    imputed_expr = np.mean(neighbor_expr, axis=1)  # Shape: (n_query, expr_dim)

    if compute_uncertainty:  # note: not implemented yet
        raise NotImplementedError("Uncertainty not implemented yet.")
        uncertainty = None
    else:
        uncertainty = None

    return imputed_expr, uncertainty


def log_median(x, axis=1):
    return np.log(np.median(np.exp(x), axis=axis))


def get_balanced_sample_indices(
    adata: AnnData,
    sample_key: str,
    random_state: int | None = None,
    return_boolean_mask: bool = True,
) -> np.ndarray:
    """
    Returns indices for a balanced subsample across sample_key groups..

    A warning is issued if any sample has fewer cells than the 10th percentile of sample sizes.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix.
    sample_key : str
        Key in `adata.obs` indicating sample identity.
    random_state : int or None, optional
        Seed for reproducible subsampling. Default is None.
    return_boolean_mask : bool, optional
        If True, returns a boolean mask indicating the selected indices.
        If False, returns the actual indices. Default is True.

    Returns
    -------
    np.ndarray
        Array of indices (from the original `adata`) corresponding to the balanced subset.
    """
    if sample_key not in adata.obs.columns:
        raise ValueError(f"'{sample_key}' not found in adata.obs")

    rng = np.random.default_rng(seed=random_state)
    sample_groups = adata.obs.groupby(sample_key).groups

    sample_sizes = np.array([len(v) for v in sample_groups.values()])
    threshold = np.percentile(sample_sizes, 10)

    small_samples = {k: len(v) for k, v in sample_groups.items() if len(v) < threshold}
    min_size = min(sample_sizes)

    if small_samples:
        warning_msg = (
            f"The following samples have fewer cells than the 10th percentile "
            f"({int(threshold)} cells):\n"
            + "\n".join([f"- {k}: {v} cells" for k, v in small_samples.items()])
            + (
                f"\n\nAll samples will be subsampled to the minimum number of cells observed "
                f"({min_size}) to ensure balanced representation. To change this behavior, set "
                "`balance_sample` to False for DE computation."
            )
        )
        warnings.warn(warning_msg, UserWarning, stacklevel=settings.warnings_stacklevel)

    balanced_indices = []
    for _group, indices in sample_groups.items():
        selected = rng.choice(indices, size=min_size, replace=False)
        balanced_indices.extend(selected)

    balanced_indices = np.array(balanced_indices)

    if return_boolean_mask:
        boolean_mask = adata.obs_names.isin(balanced_indices)
        return np.array(boolean_mask)

    else:
        return balanced_indices
