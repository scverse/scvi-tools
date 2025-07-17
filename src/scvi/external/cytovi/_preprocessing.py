from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import anndata as ad
import numpy as np

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData

from scvi import settings

from ._constants import CYTOVI_SCATTER_FEATS
from ._utils import apply_scaling, validate_layer_key, validate_marker, validate_obs_keys


def arcsinh(
    adata: AnnData,
    raw_layer_key: str = "raw",
    transformed_layer_key: str = "transformed",
    global_scaling_factor: float = 5,
    scaling_dict: dict[str, float] | None = None,
    transform_scatter: bool = False,
    inplace: bool = True,
) -> AnnData | None:
    """
    Apply the arcsinh transformation to the 'raw' layer of an AnnData object.

    Parameters
    ----------
        adata (AnnData): The AnnData object to transform.
        raw_layer_key (str): Key for the raw expression layer in the AnnData object.
        transformed_layer_key (str): Key for the layer in the AnnData object, where the transformed
            expression will be saved.
        global_scaling_factor (float): The global scaling factor to apply.
        scaling_dict (Optional[Dict[str, float]]): A dictionary of specific scaling factors
            for markers.
        transform_scatter (bool): If True, scatter features are omitted from the transformation.
        inplace (bool): If True, the transformation is applied in place. If False, a new AnnData
            object is returned.

    Returns
    -------
        Optional[AnnData]: If inplace is False, returns the transformed AnnData object.
        Otherwise, returns None.
    """
    # check arguments
    validate_layer_key(adata, raw_layer_key)
    if scaling_dict is not None:
        validate_marker(adata, scaling_dict.keys())

    # check if the transformed layer is present already
    if transformed_layer_key in adata.layers:
        msg = f"Layer {transformed_layer_key} already exists. Overwriting it."
        warnings.warn(msg, UserWarning, stacklevel=settings.warnings_stacklevel)

    # combine scaling factors into one dict
    global_dict = dict.fromkeys(adata.var_names, global_scaling_factor)

    # overwrite scaling factors if dictionary is provided
    if scaling_dict is not None:
        global_dict.update(scaling_dict)

    # apply scaling and arcsinh transformation
    transformed_layer = adata.layers[raw_layer_key].copy() / np.array(list(global_dict.values()))
    adata.layers[transformed_layer_key] = np.arcsinh(transformed_layer)

    if not transform_scatter:
        is_scatter = [marker.startswith(CYTOVI_SCATTER_FEATS) for marker in adata.var_names]

        if any(is_scatter):
            scatter_str = ", ".join(adata.var_names[is_scatter])
            msg = (
                "Detected scatter features, which are omitted for transformation. \n"
                f"Scatter features: {scatter_str}"
            )
            warnings.warn(msg, UserWarning, stacklevel=settings.warnings_stacklevel)

            adata.layers[transformed_layer_key][:, is_scatter] = adata.layers[raw_layer_key][
                :, is_scatter
            ]

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
    """
    Apply scaling to the transformed layer of an AnnData object.

    Parameters
    ----------
    adata : AnnData
        The AnnData object to scale.
    method : Literal['minmax', 'standard']
        The scaling method to use. Either 'minmax' or 'standard'.
    feature_range : Tuple[float, float]
        The desired range of the scaled data for min-max scaling.
    transformed_layer_key : str
        The key of the transformed layer in `adata.layers`.
    scaled_layer_key : str
        The key to store the scaled data in `adata.layers`.
    batch_key : str, optional
        The key for batch information in `adata.obs`. If provided, scaling is done per batch.
    feat_eps : float
        Small epsilon to adjust the feature range and avoid division by zero.
    inplace : bool
        If True, the scaling is applied in place. If False, a new AnnData object is returned.

    Returns
    -------
    Optional[AnnData]
        If inplace is False, returns the scaled AnnData object. Otherwise, returns None.
    """
    validate_layer_key(adata, transformed_layer_key)

    # check if the scaled layer is present already
    if scaled_layer_key in adata.layers:
        msg = f"Layer {scaled_layer_key} already exists. Overwriting it."
        warnings.warn(msg, UserWarning, stacklevel=settings.warnings_stacklevel)

    feature_range = (feature_range[0] + feat_eps, feature_range[1] - feat_eps)

    if batch_key:
        batches = adata.obs[batch_key].unique()
        adata_list = []
        idx = adata.obs_names

        for batch in batches:
            adata_batch = adata[adata.obs[batch_key] == batch].copy()
            scaled_data, _ = apply_scaling(
                adata_batch.layers[transformed_layer_key].copy(), method, feature_range
            )
            adata_batch.layers[scaled_layer_key] = scaled_data
            adata_list.append(adata_batch)
        adata_temp = ad.concat(adata_list, join="outer")
        adata_temp = adata_temp[idx].copy()
        scaled_data = adata_temp.layers[scaled_layer_key].copy()
        adata.layers[scaled_layer_key] = scaled_data
    else:
        scaled_data, scaler = apply_scaling(
            adata.layers[transformed_layer_key].copy(), method, feature_range
        )
        adata.layers[scaled_layer_key] = scaled_data

    return adata if not inplace else None


def register_nan_layer(
    adata: AnnData,
    mask_layer_key: str = "_nan_mask",
    scaled_layer_key: str = "scaled",
    inplace: bool = True,
) -> AnnData | None:
    """
    Add a mask layer and replace NaNs by zero in the scaled layer of an AnnData object.

    Parameters
    ----------
        adata (AnnData): The AnnData object to process.
        mask_layer_key (str): The key to store the mask layer in `adata.layers`.
        scaled_layer_key (str): The key of the scaled layer.
        inplace (bool): If True, the processing is applied in place. If False, a new
            AnnData object is returned.

    Returns
    -------
        Optional[AnnData]: If inplace is False, returns the processed AnnData object.
        Otherwise, returns None.
    """
    validate_layer_key(adata, scaled_layer_key)

    # check if the mask layer is present already
    if mask_layer_key in adata.layers:
        msg = f"Masking layer {mask_layer_key} already exists. Overwriting it."
        warnings.warn(msg, UserWarning, stacklevel=settings.warnings_stacklevel)

    # add mask layer and replace nans by zero
    adata.layers[mask_layer_key] = np.ones_like(adata.layers[scaled_layer_key])
    adata.layers[mask_layer_key][np.isnan(adata.layers[scaled_layer_key])] = 0

    # replace nans by zeroes in expression layer
    adata.layers[scaled_layer_key][np.isnan(adata.layers[scaled_layer_key])] = 0

    return adata if not inplace else None


def merge_batches(
    adata_list: list[AnnData],
    mask_layer_key: str = "_nan_mask",
    batch_key="batch",
    scaled_layer_key: str = "scaled",
    nan_layer_registration: bool = True,
) -> AnnData:
    """
    Merge batches of AnnData objects and handle missing markers.

    Parameters
    ----------
    adata_list : List[AnnData]
        List of AnnData objects to be merged.
    mask_layer_key : str, optional
        The key to store the mask layer in `adata.layers`.
    batch_key : str, optional
        The key for the batch information in `adata.obs`.
    scaled_layer_key : str, optional
        The key of the scaled layer.
    nan_layer_registration : bool, optional
        Whether to register a nan layer for imputation of missing proteins.

    Returns
    -------
    AnnData
        The merged AnnData object.
    """
    # check if there are NaNs before merging and validate layers
    for batch, adata_batch in enumerate(adata_list):
        validate_layer_key(adata_batch, scaled_layer_key)
        if np.isnan(adata_batch.layers[scaled_layer_key]).any():
            error_msg = "Nan values are present in batch {}. This will interfere with"
            " downstream processing."
            raise ValueError(error_msg.format(batch))

        # check if batch key is already present in adata
        if batch_key in adata_batch.obs:
            msg = f"Batch key {batch_key} already exists in adata. Overwriting it."
            warnings.warn(msg, UserWarning, stacklevel=settings.warnings_stacklevel)

    adata = ad.concat(adata_list, join="outer", label=batch_key, fill_value=None)
    all_markers = adata.var_names

    # register missing markers after merging
    if np.isnan(adata.layers[scaled_layer_key]).any():
        backbone_markers = list(all_markers[~np.isnan(adata.layers[scaled_layer_key]).any(axis=0)])
        backbone_str = ", ".join(backbone_markers)

        msg = (
            "Not all proteins are detected across all batches. Will generate nan_layer "
            + f"for imputation of missing proteins. \nBackbone markers: {backbone_str}"
        )
        warnings.warn(msg, UserWarning, stacklevel=settings.warnings_stacklevel)

        if nan_layer_registration:
            adata = register_nan_layer(
                adata,
                mask_layer_key=mask_layer_key,
                scaled_layer_key=scaled_layer_key,
                inplace=False,
            )

        # register which markers are present in which batch
        for batch, adata_batch in enumerate(adata_list):
            batch_name = "_batch" + "_" + str(batch)
            batch_markers = adata_batch.var_names.intersection(all_markers)
            adata.var[batch_name] = all_markers.isin(batch_markers.intersection(all_markers))

    adata.obs_names_make_unique()

    return adata


def subsample(
    adata: AnnData,
    n_obs: int = 10000,
    random_state: int = 0,
    replace: bool = False,
    groupby: str = None,
    n_obs_group: int | None = None,
) -> AnnData | None:
    """
    Subsample an AnnData object.

    Parameters
    ----------
    adata : AnnData
        The AnnData object to downsample.
    n_obs : int, optional
        The number of observations to downsample to. Default is 10000.
    random_state : int, optional
        The random state to use for the downsampling. Default is 0.
    replace : bool, optional
        If True, the downsampling is applied with replacement. If False, a new AnnData
        object is returned. Default is False.
    groupby : str, optional
        The column name in `adata.obs` to group the observations by. Default is None.
    n_obs_group : int, optional
        The number of observations to downsample to within each group. If not provided,
          it is calculated as `n_obs` divided by the number of unique groups. Default is None.

    Returns
    -------
    AnnData or None
        If `replace` is False, returns the downsampled AnnData object. Otherwise, returns None.

    Raises
    ------
    ValueError
        If the observations in `adata` are not unique.
    ValueError
        If the specified `groupby` column is not found in `adata.obs`.
    UserWarning
        If a group has fewer observations than `n_obs_group` and `replace` is False.

    """
    if len(adata.obs.index) != len(set(adata.obs.index)):
        msg = (
            "Observations are not unique. Cannot subsample. Call `.obs_names_make_unique` before."
        )
        raise ValueError(msg)

    if groupby is not None:
        if groupby not in adata.obs:
            raise ValueError(f"Group {groupby} not found in adata.obs.")
        group_cats = adata.obs[groupby].drop_duplicates().values

        if n_obs_group is None:
            n_obs_group = n_obs // len(group_cats)

        if not replace:
            for group in group_cats:
                if len(adata.obs[adata.obs[groupby] == group]) < n_obs_group:
                    msg = (
                        f"Group {group} has fewer observations than {n_obs_group} observations."
                        + " Taking all group observations. Set replace to True to sample with"
                        " replacement."
                    )
                    warnings.warn(msg, UserWarning, stacklevel=settings.warnings_stacklevel)

        index = adata.obs.groupby(groupby, as_index=False).apply(
            lambda x: x.sample(n_obs_group, random_state=random_state, replace=replace)
            if len(x) > n_obs_group
            else x
        )
        index = index.reset_index()["level_1"].to_list()
        adata_subsampled = adata[index, :].copy()
    else:
        adata_subsampled = adata[
            adata.obs.sample(n_obs, random_state=random_state, replace=replace).index, :
        ].copy()

    return adata_subsampled


def mask_markers(
    adata: AnnData,
    markers: str | list[str],
    batch_key: str = "batch",
    masked_batch: str = "1",
    nan_layer_registration: bool = True,
):
    """
    Mask specific markers in an AnnData object by removing them from a specific batch.

    Parameters
    ----------
    adata : AnnData
        The AnnData object to mask markers in.
    markers : Union[str, list[str]]
        The marker(s) to be masked. Can be a single marker as a string or a list of markers.
    batch_key : str, optional
        The key in `adata.obs` that represents the batch information, by default 'batch'.
    masked_batch : str, optional
        The batch in which the marker(s) should be masked, by default '1'.
    nan_layer_registration : bool, optional
        Whether to register a nan layer for imputation of missing proteins, by default True.

    Returns
    -------
    AnnData
        The masked AnnData object.
    """
    adata_list = []
    idx = adata.obs_names

    if isinstance(markers, str):
        markers = [markers]

    validate_marker(adata, markers)
    validate_obs_keys(adata, batch_key)

    for batch in adata.obs[batch_key].cat.categories:
        adata_batch = adata[adata.obs[batch_key] == batch].copy()

        if batch == masked_batch:
            keep_columns = ~adata_batch.var_names.isin(markers)
            adata_batch = adata_batch[:, keep_columns].copy()

        adata_list.append(adata_batch)
    adata_out = merge_batches(adata_list, nan_layer_registration=nan_layer_registration)
    adata_out = adata_out[idx].copy()

    return adata_out
