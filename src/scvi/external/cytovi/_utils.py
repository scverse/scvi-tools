from __future__ import annotations
from typing import Union

import numpy as np
import pandas as pd
import warnings
import pynndescent
from anndata import AnnData
from sklearn.preprocessing import MinMaxScaler, OneHotEncoder, StandardScaler


def validate_marker(adata: AnnData, marker: Union[str, list[str]]):
    if isinstance(marker, str):
        marker = [marker]
    for m in marker:
        if m not in adata.var_names:
            raise ValueError(f"Marker {m} not found in adata.var_names.")


def validate_obs_keys(adata: AnnData, obs_key: Union[str, list[str]]):
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


def clip_lfc_factory(min_lfc: float, max_lfc: float, pseudocount = 0):
    def clip_lfc(x, y, pseudocount = pseudocount):
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


def impute_cats_with_neighbors(rep_query, rep_ref, cat_encoded_ref, n_neighbors=5, compute_uncertainty=False):
    """Use pynndescent to find nearest neighbors and impute missing categories."""
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

def impute_expr_with_neighbors(rep_query, rep_ref, expr_data_ref, n_neighbors=5, compute_uncertainty=False):
    """Use pynndescent to find nearest neighbors and impute missing expression."""
    nn_index = pynndescent.NNDescent(rep_ref, n_neighbors=n_neighbors, metric="euclidean")

    indices, distances = nn_index.query(rep_query, k=n_neighbors)
    neighbor_expr = expr_data_ref[indices]  # Shape: (n_query, n_neighbors, expr_dim)

    imputed_expr = np.mean(neighbor_expr, axis=1)  # Shape: (n_query, expr_dim)

    if compute_uncertainty: # note: not implemented yet
        raise NotImplementedError("Uncertainty not implemented yet.")
        uncertainty = None
    else:
        uncertainty = None

    return imputed_expr, uncertainty

def log_median(x, axis = 1):
    return np.log(np.median(np.exp(x), axis = axis))


def get_balanced_sample_indices(
    adata: AnnData,
    sample_key: str,
    random_state: int | None = None,
    return_boolean_mask: bool = True
) -> np.ndarray:
    """
    Returns indices of a subsampled AnnData object such that each group
    defined by `sample_key` has an equal number of observations (balanced samples).

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
            f"The following samples have fewer cells than the 10th percentile ({int(threshold)} cells):\n" +
            "\n".join([f"- {k}: {v} cells" for k, v in small_samples.items()]) +
            f"\n\nAll samples will be subsampled to the minimum number of cells observed ({min_size}) to " +
            "ensure balanced representation. To change this behavior, set `balance_sample` to False for DE computation."
        )
        warnings.warn(warning_msg, UserWarning)

    balanced_indices = []
    for group, indices in sample_groups.items():
        selected = rng.choice(indices, size=min_size, replace=False)
        balanced_indices.extend(selected)

    balanced_indices = np.array(balanced_indices)

    if return_boolean_mask:
        boolean_mask = adata.obs_names.isin(balanced_indices)
        return np.array(boolean_mask)

    else:
        return balanced_indices


#### processing utils
import warnings
from typing import Literal, Optional, Union

import anndata as ad
import numpy as np
from anndata import AnnData
from scvi import settings

# from ._utils import apply_scaling, validate_layer_key, validate_marker, validate_obs_keys
from ._constants import CYTOVI_SCATTER_FEATS

def arcsinh(
    adata: AnnData,
    raw_layer_key: str = "raw",
    transformed_layer_key: str = "transformed",
    global_scaling_factor: float = 5,
    scaling_dict: Optional[dict[str, float]] = None,
    transform_scatter: bool = False,
    inplace: bool = True,
) -> Optional[AnnData]:
    """
    Apply the arcsinh transformation to the 'raw' layer of an AnnData object with variable scaling factors and saves in a new layer.

    Parameters
    ----------
        adata (AnnData): The AnnData object to transform.
        raw_layer_key (str): Key for the raw expression layer in the AnnData object.
        transformed_layer_key (str): Key for the layer in the AnnData object, where the transformed expression will be saved.
        global_scaling_factor (float): The global scaling factor to apply.
        scaling_dict (Optional[Dict[str, float]]): A dictionary of specific scaling factors for markers.
        transform_scatter (bool): If True, scatter features are omitted from the transformation.
        inplace (bool): If True, the transformation is applied in place. If False, a new AnnData object is returned.

    Returns
    -------
        Optional[AnnData]: If inplace is False, returns the transformed AnnData object. Otherwise, returns None.
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
    global_dict = {marker: global_scaling_factor for marker in adata.var_names}

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
            msg = f"Detected scatter features, which are omitted for transformation. \nScatter features: {scatter_str}"
            warnings.warn(msg, UserWarning, stacklevel=settings.warnings_stacklevel)

            adata.layers[transformed_layer_key][:, is_scatter] = adata.layers[raw_layer_key][:, is_scatter]

    return adata if not inplace else None


def logp(
    adata: AnnData,
    raw_layer_key: str = "raw",
    transformed_layer_key: str = "transformed",
    offset: float = 1.0,
    transform_scatter: bool = False,
    inplace: bool = True,
) -> Optional[AnnData]:
    """
    Apply log transformation to the 'raw' layer of an AnnData object.

    Parameters
    ----------
    adata : AnnData
        The AnnData object to normalize.
    raw_layer_key : str, optional
        The key of the raw layer to be transformed. Default is "raw".
    transformed_layer_key : str, optional
        The key to store the transformed data in `adata.layers`. Default is "transformed".
    offset : float, optional
        The offset value to add before taking the logarithm. Default is 1.0.
    transform_scatter : bool, optional
        If True, scatter features are omitted from transformation. Default is False.
    inplace : bool, optional
        If True, the normalization is applied in place. If False, a new AnnData object is returned. Default is True.

    Returns
    -------
    Optional[AnnData]
        If inplace is False, returns the normalized AnnData object. Otherwise, returns None.
    """
    validate_layer_key(adata, raw_layer_key)

    # check if the transformed layer is present already
    if transformed_layer_key in adata.layers:
        msg = f"Layer {transformed_layer_key} already exists. Overwriting it."
        warnings.warn(msg, UserWarning, stacklevel=settings.warnings_stacklevel)

    adata.layers[transformed_layer_key] = adata.layers[raw_layer_key].copy()
    adata.layers[transformed_layer_key] = adata.layers[transformed_layer_key].astype('float')
    adata.layers[transformed_layer_key] += offset
    adata.layers[transformed_layer_key] = np.log(adata.layers[transformed_layer_key])

    if not transform_scatter:
        is_scatter = [marker.startswith(CYTOVI_SCATTER_FEATS) for marker in adata.var_names]

        if any(is_scatter):
            scatter_str = ", ".join(adata.var_names[is_scatter])
            msg = f"Detected scatter features, which are omitted for transformation. \nScatter features: {scatter_str}"
            warnings.warn(msg, UserWarning, stacklevel=settings.warnings_stacklevel)

            adata.layers[transformed_layer_key][:, is_scatter] = adata.layers[raw_layer_key][:, is_scatter]

    return adata if not inplace else None


def scale(
    adata: AnnData,
    method: Literal['minmax', 'standard'] = "minmax",
    feature_range: tuple[float, float] = (0.0, 1.0),
    transformed_layer_key: str = "transformed",
    scaled_layer_key: str = "scaled",
    batch_key: str = None,
    feat_eps: float = 1e-6,
    inplace: bool = True,
) -> Optional[AnnData]:
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
            scaled_data, _ = apply_scaling(adata_batch.layers[transformed_layer_key].copy(), method, feature_range)
            adata_batch.layers[scaled_layer_key] = scaled_data
            adata_list.append(adata_batch)
        adata_temp = ad.concat(adata_list, join="outer")
        adata_temp = adata_temp[idx].copy()
        scaled_data = adata_temp.layers[scaled_layer_key].copy()
        adata.layers[scaled_layer_key] = scaled_data
    else:
        scaled_data, scaler = apply_scaling(adata.layers[transformed_layer_key].copy(), method, feature_range)
        adata.layers[scaled_layer_key] = scaled_data

    return adata if not inplace else None


def register_nan_layer(
    adata: AnnData,
    mask_layer_key: str = "_nan_mask",
    scaled_layer_key: str = "scaled",
    inplace: bool = True,
) -> Optional[AnnData]:
    """
    Add a mask layer and replace NaNs by zero in the scaled layer of an AnnData object.

    Parameters
    ----------
        adata (AnnData): The AnnData object to process.
        mask_layer_key (str): The key to store the mask layer in `adata.layers`.
        scaled_layer_key (str): The key of the scaled layer.
        inplace (bool): If True, the processing is applied in place. If False, a new AnnData object is returned.

    Returns
    -------
        Optional[AnnData]: If inplace is False, returns the processed AnnData object. Otherwise, returns None.
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
    adata_list: list[AnnData], mask_layer_key: str = "_nan_mask", batch_key = 'batch', scaled_layer_key: str = "scaled", nan_layer_registration: bool = True,
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
            error_msg = "Nan values are present in batch {}. This will interfere with downstream processing."
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
                adata, mask_layer_key=mask_layer_key, scaled_layer_key=scaled_layer_key, inplace=False
            )

        # register which markers are present in which batch
        for batch, adata_batch in enumerate(adata_list):
            batch_name = "_batch" + "_" + str(batch)
            batch_markers = adata_batch.var_names.intersection(all_markers)
            adata.var[batch_name] = all_markers.isin(batch_markers.intersection(all_markers))

    adata.obs_names_make_unique()

    return adata


def subsample(
    adata: AnnData, n_obs: int = 10000, random_state: int = 0, replace: bool = False, groupby: str = None, n_obs_group: Optional[int] = None,
) -> Optional[AnnData]:
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
        If True, the downsampling is applied with replacement. If False, a new AnnData object is returned. Default is False.
    groupby : str, optional
        The column name in `adata.obs` to group the observations by. Default is None.
    n_obs_group : int, optional
        The number of observations to downsample to within each group. If not provided, it is calculated as `n_obs` divided by the number of unique groups. Default is None.

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
        msg = "Observations are not unique. Cannot subsample. Call `.obs_names_make_unique` before."
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
                        + " Taking all group observations. Set replace to True to sample with replacement."
                    )
                    warnings.warn(msg, UserWarning, stacklevel=settings.warnings_stacklevel)

        index = adata.obs.groupby(groupby, as_index=False).apply(
            lambda x: x.sample(n_obs_group, random_state=random_state, replace=replace) if len(x) > n_obs_group else x
        )
        index = index.reset_index()["level_1"].to_list()
        adata_subsampled = adata[index, :].copy()
    else:
        adata_subsampled = adata[adata.obs.sample(n_obs, random_state=random_state, replace=replace).index, :].copy()

    return adata_subsampled


def mask_markers(
    adata: AnnData, markers: Union[str, list[str]], batch_key: str = 'batch', masked_batch: str = '1', nan_layer_registration: bool = True):
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



#### plotting utils
import math
from typing import Union

import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from anndata import AnnData

# from cytovi._utils import validate_layer_key, validate_marker, validate_obs_keys
# from cytovi.pp.cyto_pp import subsample


def histogram(
    adata: ad.AnnData,
    marker: Union[str, list[str]] = "all",
    groupby: str = None,
    layer_key: str = "raw",
    downsample: bool = True,
    n_obs: int = 10000,
    col_wrap: int = None,
    tight_layout: bool = True,
    save: Union[bool, str] = None,
    return_plot: bool = False,
    kde_kwargs: dict = None,
    **kwargs,
):
    """
    Create a FacetGrid of histograms for specified markers in AnnData.

    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix.

    marker : Union[str, List[str]], optional
        Names of markers to plot. 'all' to plot all markers.

    groupby : str, optional
        Key for grouping or categorizing the data. E.g. key for batch.

    layer_key : str, optional
        Key for the layer in AnnData.

    downsample : bool, optional
        Whether to downsample the data if there are too many observations.

    n_obs : int, optional
        Number of observations to subsample if downsample is True.

    col_wrap : int, optional
        Number of columns to wrap the plots. If None, it is calculated based on the number of markers.

    tight_layout : bool, optional
        Whether to use tight layout for the plot.

    save : Union[bool, str], optional
        If True, the plot is saved as "marker_histogram.png". If a string is provided, the plot is saved with the given filename.

    return_plot : bool, optional
        Whether to return the FacetGrid object.

    kde_kwargs : dict, optional
        Additional keyword arguments to pass to Seaborn's kdeplot.

    **kwargs : additional keyword arguments
        Additional arguments to pass to Seaborn's FacetGrid.

    Returns
    -------
    None or sns.FacetGrid
        If return_plot is True, returns the FacetGrid object. Otherwise, returns None.

    Example:
    ----------
    # Plot density plots for specific markers
    cytovi.pl.histogram(adata, marker=['CD3', 'CD4'], group_by='Condition')

    # Plot density plots for all markers
    cytovi.pl.histogram(adata, marker='all', group_by='Batch')
    """
    if kde_kwargs is None:
        kde_kwargs = {}

    if marker == "all":
        marker = adata.var_names
    elif isinstance(marker, str):
        marker = [marker]

    validate_marker(adata, marker)
    validate_obs_keys(adata, groupby)
    validate_layer_key(adata, layer_key)

    # subsample if too many observations
    if downsample and adata.n_obs > 10000:
        adata = subsample(adata, n_obs=n_obs, groupby=groupby)

    num_plots = len(marker)

    if col_wrap is None:
        col_wrap = math.ceil(math.sqrt(num_plots))

    data_plot = adata[:, marker].to_df(layer=layer_key)

    if groupby is not None:
        data_plot[groupby] = adata.obs[groupby]

    data_plot_melt = data_plot.melt(id_vars=groupby, var_name='variable', value_name='value')

    # generate the plot
    g = sns.FacetGrid(
        data_plot_melt, col="variable", hue=groupby, col_wrap=col_wrap, sharey=False, sharex=False, **kwargs
    )
    g.map(sns.kdeplot, "value", fill=True, **kde_kwargs)
    g.set_titles("{col_name}")
    g.set(yticks=[])
    g.set_axis_labels("", "")
    g.add_legend()
    g.fig.text(0, 0.5, "Density", va="center", ha="center", rotation="vertical")

    if tight_layout:
        g.fig.tight_layout()

    if save is not None:
        if save is True:
            save = "marker_histogram.png"
        g.savefig(save)

    if return_plot:
        return g


def biaxial(
    adata: AnnData,
    marker_x: Union[str, list[str]] = None,
    marker_y: Union[str, list[str]] = None,
    color: str = None,
    n_bins: int = 10,
    layer_key: str = "raw",
    downsample: bool = True,
    n_obs: int = 10000,
    sample_color_groups: bool = False,
    save: Union[bool, str] = None,
    kde: bool = True,
    kde_kwargs: dict = None,
    scatter_kwargs: dict = None,
    **kwargs,
):
    """
    Create a PairGrid of biaxial (scatter and density) plots for specified markers in AnnData.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.

    marker_x : Union[str, List[str]], optional
        Variable name(s) to be plotted on the x-axis.

    marker_y : Union[str, List[str]], optional
        Variable name(s) to be plotted on the y-axis.

    color : str, optional
        Variable name to be used for coloring the scatter plots.

    n_bins : int, optional
        Number of levels for density contours in kdeplot.

    layer_key : str, optional
        Key for the layer in AnnData.

    downsample : bool, optional
        Whether to downsample the data if there are too many observations.

    n_obs : int, optional
        Number of observations to keep if downsampling is enabled.

    sample_color_groups : bool, optional
        Whether to sample observations within each color group if downsampling is enabled.

    save : Union[bool, str], optional
        If True, save the plot as "marker_histogram.png". If a string is provided, save the plot with the given filename.

    kde : bool, optional
        Whether to include density contours in the plot.

    kde_kwargs : dict, optional
        Additional keyword arguments to pass to Seaborn's kdeplot.

    scatter_kwargs : dict, optional
        Additional keyword arguments to pass to Seaborn's scatterplot.

    **kwargs : additional keyword arguments
        Additional arguments to pass to Seaborn's PairGrid.

    Returns
    -------
    None

    Example
    -------
    # Plot biaxial plots for specific markers
    cytovi.pl.biaxial(adata, marker_x='CD3', marker_y='CD4', color='Condition')

    # Plot biaxial plots for multiple markers
    cytovi.pl.biaxial(adata, marker_x=['CD8', 'CD20'], marker_y='CD56', color='batch')
    """
    if isinstance(marker_x, str):
        marker_x = [marker_x]
    if isinstance(marker_y, str):
        marker_y = [marker_y]

    if marker_x is None and marker_y is None:
        raise ValueError("At least one of marker_x or marker_y must be specified.")
    elif marker_x is None:
        marker_x = [*adata.var_names]
        marker_x = list(set(marker_x) - set(marker_y))
    elif marker_y is None:
        marker_y = [*adata.var_names]
        marker_y = list(set(marker_y) - set(marker_x))
    else:
        marker_x = list(set(marker_x) - set(marker_y))

    validate_marker(adata, marker_x)
    validate_marker(adata, marker_y)
    validate_obs_keys(adata, color)
    validate_layer_key(adata, layer_key)

    if kde_kwargs is None:
        kde_kwargs = {}

    if scatter_kwargs is None:
        scatter_kwargs = {}

    # subsample if too many observations
    if downsample and adata.n_obs > 10000:
        if color is not None and sample_color_groups is True:
            adata = subsample(adata, n_obs=n_obs, groupby=color)
        else:
            adata = subsample(adata, n_obs=n_obs)

    marker = marker_x + marker_y

    data_plot = adata[:, marker].to_df(layer=layer_key)

    if color is not None:
        data_plot[color] = adata.obs[color]

    g = sns.PairGrid(data_plot, x_vars=marker_x, y_vars=marker_y, hue=color, **kwargs)

    if kde is True:
        g.map(sns.kdeplot, levels=n_bins, **kde_kwargs)

    g.map(sns.scatterplot, s=5, **scatter_kwargs)
    g.add_legend()

    if save is not None:
        if save is True:
            save = "marker_histogram.png"
        g.savefig(save)

    plt.show()
