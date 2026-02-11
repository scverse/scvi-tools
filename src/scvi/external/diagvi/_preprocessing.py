from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import anndata as ad
import numpy as np
import scipy.sparse as sp

if TYPE_CHECKING:
    from typing import Literal
    from anndata import AnnData

from ._utils import apply_scaling, validate_layer_key, validate_marker

logger = logging.getLogger(__name__)


def transform_arcsinh(
    adata: AnnData,
    raw_layer_key: str = "raw",
    transformed_layer_key: str = "transformed",
    global_scaling_factor: float = 5,
    scaling_dict: dict[str, float] | None = None,
    inplace: bool = True,
) -> AnnData | None:
    """Apply the arcsinh transformation to the raw layer of an AnnData object.
    
    This function is adapted from CytoVI's preprocessing utilities.

    Parameters
    ----------
    adata
        AnnData object to transform.
    raw_layer_key
        Key for the raw expression layer in the AnnData object.
    transformed_layer_key
        Key for the layer where the transformed expression will be saved.
    global_scaling_factor
        Global scaling factor to apply.
    scaling_dict
        Dictionary of marker-specific scaling factors.
    inplace
        If True, apply the transformation in place; if False, return a new AnnData.

    Returns
    -------
    If ``inplace`` is False, returns the transformed ``AnnData`` object; otherwise, nothing.
    """
    validate_layer_key(adata, raw_layer_key)
    if scaling_dict is not None:
        validate_marker(adata, scaling_dict.keys())

    # Check if the transformed layer is already present
    if transformed_layer_key in adata.layers:
        logger.warning(f"Layer {transformed_layer_key} already exists. Overwriting it.")

    # Combine scaling factors into one dict
    global_dict = dict.fromkeys(adata.var_names, global_scaling_factor)

    # Overwrite scaling factors if dictionary is provided
    if scaling_dict is not None:
        global_dict.update(scaling_dict)

    # Convert layer to dense if sparse
    raw_layer = adata.layers[raw_layer_key]
    if sp.issparse(raw_layer):
        raw_layer = raw_layer.toarray()
        was_sparse = True
    else:
        raw_layer = raw_layer.copy()
        was_sparse = False

    # Apply scaling factor and arcsinh transformation
    transformed_layer = raw_layer / np.array(list(global_dict.values()))
    transformed_layer = np.arcsinh(transformed_layer)
    
    # Convert back to sparse if original was sparse
    if was_sparse:
        adata.layers[transformed_layer_key] = sp.csr_matrix(transformed_layer)
    else:
        adata.layers[transformed_layer_key] = transformed_layer

    return adata if not inplace else None


def scale(
    adata: AnnData,
    method: Literal["minmax", "standard"] = "minmax",
    feature_range: tuple[float, float] = (0.0, 1.0),
    transformed_layer_key: str = "transformed",
    scaled_layer_key: str = "scaled",
    batch_key: str = None,
    feat_eps: float = 1e-6,
    inplace: bool = True,
) -> AnnData | None:
    """Apply scaling to the transformed layer of an AnnData object.
    
    This function is adapted from CytoVI's preprocessing utilities.

    Parameters
    ----------
    adata
        The AnnData object to scale.
    method
        Scaling method: ``'minmax'`` or ``'standard'``.
    feature_range
        Desired range for min–max scaling as ``(low, high)``.
    transformed_layer_key
        Key of the source (transformed) layer in ``adata.layers``.
    scaled_layer_key
        Key to store the scaled data in ``adata.layers``.
    batch_key
        Key in ``adata.obs`` for batch labels; if provided, scaling is performed per batch.
    feat_eps
        Small epsilon added to the feature range to avoid division by zero.
    inplace
        If ``True``, apply in place; if ``False``, return a new ``AnnData``.

    Returns
    -------
    If ``inplace`` is ``False``, returns the scaled ``AnnData`` object; otherwise, nothing.
    """
    validate_layer_key(adata, transformed_layer_key)

    # Check if the scaled layer is present already
    if scaled_layer_key in adata.layers:
        logger.warning(f"Layer {scaled_layer_key} already exists. Overwriting it.")

    # Set feature range with epsilon
    feature_range = (feature_range[0] + feat_eps, feature_range[1] - feat_eps)

    # Check sparsity before branching
    was_sparse = sp.issparse(adata.layers[transformed_layer_key])


    if batch_key:
        # Scale per batch
        batches = adata.obs[batch_key].unique()
        adata_list = []
        idx = adata.obs_names

        for batch in batches:
            adata_batch = adata[adata.obs[batch_key] == batch].copy()

            # Convert layer to dense if sparse
            transformed_layer = adata_batch.layers[transformed_layer_key]
            if sp.issparse(transformed_layer):
                transformed_layer = transformed_layer.toarray()
            else:
                transformed_layer = transformed_layer.copy()

            scaled_data, _ = apply_scaling(
                transformed_layer, method, feature_range
            )
            adata_batch.layers[scaled_layer_key] = scaled_data
            adata_list.append(adata_batch)
        
        adata_temp = ad.concat(adata_list, join="outer")
        adata_temp = adata_temp[idx].copy()
        scaled_data = adata_temp.layers[scaled_layer_key].copy()
        adata.layers[scaled_layer_key] = scaled_data
    
    else:
        # Scale all data together

        # Convert layer to dense if sparse
        transformed_layer = adata.layers[transformed_layer_key]
        if sp.issparse(transformed_layer):
            transformed_layer = transformed_layer.toarray()
        else:
            transformed_layer = transformed_layer.copy()

        scaled_data, scaler = apply_scaling(
            transformed_layer, method, feature_range
        )
        adata.layers[scaled_layer_key] = scaled_data

    # Convert back to sparse if original was sparse
    if was_sparse:
        adata.layers[scaled_layer_key] = sp.csr_matrix(adata.layers[scaled_layer_key])

    return adata if not inplace else None
