from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

import scvi
from scvi.external.drvi import DRVI


def iterate_dimensions(
    latent_dims: np.ndarray,
    latent_min: np.ndarray,
    latent_max: np.ndarray,
    n_steps: int = 10 * 2,
    n_samples: int = 100,
) -> AnnData:
    """Generate systematic traversal data for latent dimensions.

    This function creates a systematic grid of latent space traversals by
    generating combinations of dimension IDs, step IDs, and sample IDs.
    It creates a sparse matrix representation of the traversal vectors.

    Parameters
    ----------
    latent_dims
        Array of latent dimension indices to traverse.
    latent_min
        Minimum values for each latent dimension (typically negative).
        Must have same length as `latent_dims`.
    latent_max
        Maximum values for each latent dimension (typically positive).
        Must have same length as `latent_dims`.
    n_steps
        Number of steps in the traversal. Must be even (half negative, half positive).
    n_samples
        Number of samples to generate for each step.

    Returns
    -------
    AnnData
        AnnData object containing the traversal vectors in `.X` and metadata in `.obs`.
        The `.obs` contains:
        - `original_order`: Original index of each row
        - `dim_id`: Which latent dimension this row corresponds to
        - `sample_id`: Sample identifier (0 to n_samples-1)
        - `step_id`: Step identifier (0 to n_steps-1)
        - `span_value`: The actual latent value for this step

    Raises
    ------
    AssertionError
        If `n_steps` is not even.

    Notes
    -----
    The function creates a systematic traversal where:
    - For each latent dimension, it generates n_steps steps
    - Each step has n_samples samples
    - The first half of steps go from latent_min to 0
    - The second half of steps go from 0 to latent_max
    - The result is a sparse matrix of shape (n_latent * n_steps * n_samples, n_latent)

    The traversal values are linearly interpolated between the min/max bounds,
    ensuring smooth coverage of the latent space for each dimension.

    Examples
    --------
    >>> # Basic traversal
    >>> latent_dims = np.array([0, 1, 2])
    >>> latent_min = np.array([-2, -1, -3])
    >>> latent_max = np.array([2, 1, 3])
    >>> traverse_data = iterate_dimensions(latent_dims, latent_min, latent_max)
    >>> print(f"Shape: {traverse_data.X.shape}")
    >>> print(f"Unique dimensions: {traverse_data.obs['dim_id'].unique()}")
    """
    assert n_steps % 2 == 0, "n_steps must be even"
    # Sometimes it is negligibly not like below
    # assert np.all(latent_min <= 0) & np.all(latent_max >= 0)

    dim_ids = (
        latent_dims.reshape(-1, 1, 1)
        * np.ones(n_steps).astype(int).reshape(1, -1, 1)
        * np.ones(n_samples).astype(int).reshape(1, 1, -1)
    ).reshape(-1)  # n_latent * n_steps * n_samples
    sample_ids = (
        np.ones(len(latent_dims)).astype(int).reshape(-1, 1, 1)
        * np.ones(n_steps).astype(int).reshape(1, -1, 1)
        * np.arange(n_samples).astype(int).reshape(1, 1, -1)
    ).reshape(-1)  # n_latent * n_steps * n_samples
    step_ids = (
        np.ones(len(latent_dims)).astype(int).reshape(-1, 1, 1)
        * np.arange(n_steps).astype(int).reshape(1, -1, 1)
        * np.ones(n_samples).astype(int).reshape(1, 1, -1)
    ).reshape(-1)  # n_latent * n_steps * n_samples
    span_values = (
        np.concatenate(
            [
                np.linspace(latent_min, latent_min * 0.0, num=int(n_steps / 2)),
                np.linspace(latent_max * 0.0, latent_max, num=int(n_steps / 2)),
            ],
            axis=0,
        ).T.reshape(-1, 1)
        * np.ones(n_samples).reshape(1, -1)
    ).reshape(-1)  # n_latent * n_steps * n_samples

    span_vectors = sparse.coo_matrix(
        (span_values, (np.arange(len(dim_ids)), dim_ids)), dtype=np.float32
    )
    span_vectors = span_vectors.tocsr()  # n_latent * n_steps * n_samples x n_latent

    span_adata = AnnData(
        X=span_vectors,
        obs=pd.DataFrame(
            {
                "original_order": np.arange(span_vectors.shape[0]),
                "dim_id": dim_ids,
                "sample_id": sample_ids,
                "step_id": step_ids,
                "span_value": span_values,
            }
        ),
    )

    return span_adata


def make_traverse_adata(
    model: DRVI,
    embed: AnnData,
    n_steps: int = 10 * 2,
    n_samples: int = 100,
    noise_formula: callable = lambda x: x / 2,
    max_noise_std: float = 0.2,
    copy_adata_var_info: bool = True,
    **kwargs,
) -> AnnData:
    """Create traversal AnnData by decoding latent space traversals.

    This function generates systematic traversals through the latent space
    and decodes them to observe the effects on gene expression. It creates
    both control (baseline) and effect (traversal) conditions.

    Parameters
    ----------
    model
        Trained DRVI model for decoding latent representations.
    embed
        AnnData object containing latent dimension statistics in `.var`.
        Must have columns: `original_dim_id`, `min`, `max`, `std`.
    n_steps
        Number of steps in the traversal. Must be even (half negative, half positive).
    n_samples
        Number of samples to generate for each step.
    noise_formula
        Function to compute noise standard deviation from dimension std values.
        Should take a numpy array and return a numpy array.
    max_noise_std
        Maximum allowed noise standard deviation.
    copy_adata_var_info
        Whether to copy variable information from the original model's AnnData.
    **kwargs
        Additional keyword arguments passed to `model.decode_latent_samples`.

    Returns
    -------
    AnnData
        AnnData object containing the traversal results with the following structure:
        - `.X`: Difference between effect and control conditions (main result)
        - `.layers['control']`: Control condition gene expression (noise only)
        - `.layers['effect']`: Effect condition gene expression (noise + traversal)
        - `.obsm['control_latent']`: Control condition latent values (noise only)
        - `.obsm['effect_latent']`: Effect condition latent values (noise + traversal)
        - `.obsm['cat_covs']`: Categorical covariates (if model uses them)
        - `.obs['lib_size']`: Library size for each sample (set to 10,000)
        - `.obs`: Metadata including dim_id, sample_id, step_id, span_value

    Raises
    ------
    NotImplementedError
        If the model has continuous covariates (not yet supported).
    ValueError
        If required columns are missing from `embed.var`.

    Notes
    -----
    The function performs the following steps:
    1. Generates systematic traversal vectors using `iterate_dimensions`
    2. Creates baseline noise for each sample based on dimension std values
    3. Generates random categorical covariates for individual samples if the model uses them
    4. Sets library size to 10,000 for all samples (typical for single-cell data)
    5. Decodes both control (noise only) and effect (noise + traversal) conditions
    6. Returns the difference between effect and control as the main data

    The traversal goes from negative to positive values for each dimension,
    allowing observation of how each dimension affects gene expression.

    The control condition contains only noise, while the effect condition
    contains noise plus the systematic traversal. The difference shows
    the pure effect of the latent dimension traversal on gene expression.

    The noise is generated based on the standard deviation of each latent dimension,
    scaled by the `noise_formula` function and clipped to `max_noise_std`.
    """
    # Generate random delta vectors for each dimension
    span_adata = iterate_dimensions(
        latent_dims=embed.var["original_dim_id"].values,
        latent_min=embed.var["min"].values,
        latent_max=embed.var["max"].values,
        n_steps=n_steps,
        n_samples=n_samples,
    )

    # make baseline noise for each sample
    noise_std = noise_formula(embed.var["std"].values).clip(0, max_noise_std).reshape(1, -1)
    sample_wise_noises = np.random.randn(n_samples, embed.n_vars).astype(np.float32) * noise_std
    noise_vector = sample_wise_noises[span_adata.obs["sample_id"]]

    # make categorical covariates for each sample
    if model.adata_manager.get_state_registry(scvi.REGISTRY_KEYS.CAT_COVS_KEY):
        n_cats_per_key = model.adata_manager.get_state_registry(
            scvi.REGISTRY_KEYS.CAT_COVS_KEY
        ).n_cats_per_key
        sample_wise_cats = np.stack(
            [np.random.randint(0, n_cat, size=n_samples) for n_cat in n_cats_per_key], axis=1
        )
        cat_vector = sample_wise_cats[span_adata.obs["sample_id"]]
    else:
        cat_vector = None

    if model.adata_manager.get_state_registry(scvi.REGISTRY_KEYS.CONT_COVS_KEY):
        raise NotImplementedError(
            "Interpretability of models with continuous covariates are not implemented yet."
        )

    # lib size
    lib_vector = np.ones(n_samples) * 1e4
    lib_vector = lib_vector[span_adata.obs["sample_id"]]

    # Control and effect latent data
    control_data = noise_vector
    effect_data = noise_vector + span_adata.X.toarray()

    print("traversing latent ...")

    print(f"Input latent shape: control: {control_data.shape}, effect: {effect_data.shape}")
    # control and effect in mean parameter space
    control_mean_param = model.decode_latent_samples(
        control_data, lib=lib_vector, cat_values=cat_vector, cont_values=None, **kwargs
    )
    effect_mean_param = model.decode_latent_samples(
        effect_data, lib=lib_vector, cat_values=cat_vector, cont_values=None, **kwargs
    )
    print(
        f"Output mean param shape: control: {control_mean_param.shape}, effect: {effect_mean_param.shape}"
    )

    if copy_adata_var_info:
        traverse_adata_var = model.adata.var.copy()
    else:
        traverse_adata_var = pd.DataFrame(index=model.adata.var_names)
    traverse_adata_var["original_order"] = np.arange(effect_mean_param.shape[1])

    traverse_adata = AnnData(
        X=effect_mean_param - control_mean_param,
        obs=span_adata.obs,
        var=traverse_adata_var,
    )
    traverse_adata.layers["control"] = control_mean_param
    traverse_adata.layers["effect"] = effect_mean_param
    traverse_adata.obsm["control_latent"] = control_data
    traverse_adata.obsm["effect_latent"] = effect_data
    if cat_vector is not None:
        traverse_adata.obsm["cat_covs"] = cat_vector
    traverse_adata.obs["lib_size"] = lib_vector

    return traverse_adata


def get_dimensions_of_traverse_data(traverse_adata: AnnData) -> tuple[int, int, int, int]:
    """Get the dimensions of traversal data.

    This function extracts the key dimensions from a traversal AnnData object
    created by `make_traverse_adata` or `traverse_latent`.

    Parameters
    ----------
    traverse_adata
        AnnData object created by traversal functions.

    Returns
    -------
    tuple[int, int, int, int]
        A tuple containing:
        - n_latent: Number of latent dimensions traversed
        - n_steps: Number of steps in each traversal
        - n_samples: Number of samples per step
        - n_vars: Number of variables (genes) in the data

    Raises
    ------
    KeyError
        If required columns are missing from `traverse_adata.obs`.

    Notes
    -----
    The function expects the following columns in `traverse_adata.obs`:
    - `dim_id`: Latent dimension identifier
    - `step_id`: Step identifier within each dimension
    - `sample_id`: Sample identifier within each step

    Examples
    --------
    >>> # Get dimensions of traversal data
    >>> n_latent, n_steps, n_samples, n_vars = get_dimensions_of_traverse_data(traverse_adata)
    >>> print(f"Traversed {n_latent} dimensions with {n_steps} steps and {n_samples} samples")
    >>> print(f"Data has {n_vars} variables")
    """
    # Get the number of latent dimensions, steps, samples, and vars
    n_latent = traverse_adata.obs["dim_id"].nunique()
    n_steps = traverse_adata.obs["step_id"].nunique()
    n_samples = traverse_adata.obs["sample_id"].nunique()
    n_vars = traverse_adata.n_vars
    return n_latent, n_steps, n_samples, n_vars


def traverse_latent(
    model: DRVI,
    embed: AnnData,
    n_steps: int = 10 * 2,
    n_samples: int = 100,
    copy_adata_var_info: bool = True,
    **kwargs,
) -> AnnData:
    """Perform latent space traversal and enrich with metadata.

    This function generates systematic traversals through the latent space
    and decodes them to observe the effects on gene expression. It creates
    both control (baseline) and effect (traversal) conditions.
    Additionally, it enriches the data with dimension-specific metadata
    like titles, vanished status, and ordering.

    Parameters
    ----------
    model
        Trained DRVI model for decoding latent representations.
    embed
        AnnData object containing latent dimension statistics in `.var`.
        Must have columns: `original_dim_id`, `min`, `max`, `std`, `title`, `vanished`, `order`.
    n_steps
        Number of steps in the traversal. Must be even.
    n_samples
        Number of samples to generate for each step.
    copy_adata_var_info
        Whether to copy variable information from the original model's AnnData.
    **kwargs
        Additional keyword arguments passed to `make_traverse_adata`.

    Returns
    -------
    AnnData
        AnnData object containing the traversal results with enriched metadata.
        In addition to the structure returned by `make_traverse_adata`, the `.obs`
        also contains:
        - `title`: Dimension titles from `embed.var['title']`
        - `vanished`: Vanished status from `embed.var['vanished']`
        - `order`: Dimension ordering from `embed.var['order']`

    Raises
    ------
    ValueError
        If required columns are missing from `embed.var`.

    Notes
    -----
    The function performs the following steps:
    1. Generates systematic traversal vectors using `iterate_dimensions`
    2. Creates baseline noise for each sample based on dimension std values
    3. Generates random categorical covariates for individual samples if the model uses them
    4. Decodes both control (noise only) and effect (noise + traversal) conditions
    5. Returns the difference between effect and control as the main data
    6. Enriches the traversal data with dimension-specific information like titles, vanished status, and ordering.

    The function expects the following columns in `embed.var`:
    - `original_dim_id`: Original dimension indices
    - `min`, `max`, `std`: Dimension statistics
    - `title`: Human-readable dimension titles
    - `vanished`: Boolean indicating vanished dimensions
    - `order`: Dimension ordering

    Examples
    --------
    >>> # Basic traversal
    >>> traverse_data = traverse_latent(model, embed)
    >>> # Traversal with custom parameters
    >>> traverse_data = traverse_latent(model, embed, n_steps=30, n_samples=50)
    """
    if "original_dim_id" not in embed.var:
        raise ValueError(
            'Column "original_dim_id" not found in `embed.var`. Please run `set_latent_dimension_stats` to set vanished status.'
        )

    traverse_adata = make_traverse_adata(
        model=model,
        embed=embed,
        n_steps=n_steps,
        n_samples=n_samples,
        copy_adata_var_info=copy_adata_var_info,
        **kwargs,
    )

    # enrich traverse_adata with the additional info
    for col in ["title", "vanished", "order"]:
        mapping = dict(
            zip(embed.var["original_dim_id"].values, embed.var[col].values, strict=False)
        )
        traverse_adata.obs[col] = traverse_adata.obs["dim_id"].map(mapping)

    return traverse_adata
