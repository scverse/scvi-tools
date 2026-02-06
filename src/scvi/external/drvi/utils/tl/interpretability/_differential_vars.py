from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import scipy

from scvi.external.drvi import DRVI
from scvi.external.drvi.utils.tl.interpretability._latent_traverse import (
    get_dimensions_of_traverse_data,
    traverse_latent,
)

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData


def find_differential_effects(
    traverse_adata: AnnData,
    method: Literal["max_possible", "min_possible"] = "max_possible",
    key_added: str = "effect",
    add_to_counts: float = 0.1,
    relax_max_by: float = 0.0,
) -> None:
    """Find differential effects in latent space traversal data.

    This function analyzes the differential effects between control and effect
    conditions in traversal data to identify genes that respond to latent
    dimension changes. It supports two methods for calculating differential effects.

    The "max_possible" method is the simple log-fold-change effect between effect and control.
    The "min_possible" method is more conservative and considers the maximum possible effect of
    other dimensions to normalize the log-fold-change effect of the current dimension.

    Parameters
    ----------
    traverse_adata
        AnnData object created by `traverse_latent` or `make_traverse_adata`.
        Must contain `.layers['control']` and `.layers['effect']`.
    method
        Method for calculating differential effects:
        - "max_possible": Simple log-fold-change between effect and control conditions.
          This is the direct difference in log-space and represents the maximum
          possible effect a latent dimension can have on gene expression.
        - "min_possible": Conservative estimate that normalizes the effect by considering
          the maximum possible effects from other dimensions. Although the effect in count space
          is deterministic, the effect in log-space is not. This accounts for that, so changes in
          one dimension may constrain the possible change in other dimensions.
    key_added
        Prefix for the keys added to `traverse_adata.uns` and `traverse_adata.varm`.
        Results will be stored with keys like `{key_added}_traverse_effect_stepwise`.
    add_to_counts
        Small value added to counts to avoid log(0) issues in log-space calculations.
        This pseudo-count ensures numerical stability when computing log-fold-changes.
    relax_max_by
        Relaxation factor for the maximum possible effect calculation
        (only used with "min_possible" method).

    Returns
    -------
    None
        Results are stored in `traverse_adata`:
        - `.uns[f"{key_added}_traverse_effect_stepwise"]`: Stepwise effects for each
          latent dimension, step, and gene (shape: n_latent × n_steps × n_vars)
        - `.varm[f"{key_added}_traverse_effect_pos"]`: Maximum positive effects per
          dimension and gene (DataFrame with genes as rows, dimensions as columns)
        - `.varm[f"{key_added}_traverse_effect_neg"]`: Maximum negative effects per
          dimension and gene (DataFrame with genes as rows, dimensions as columns)
        - `.uns[f"{key_added}_traverse_effect_pos_dim_ids"]`: Array of dimension IDs
          for positive effects
        - `.uns[f"{key_added}_traverse_effect_neg_dim_ids"]`: Array of dimension IDs
          for negative effects

    Raises
    ------
    AssertionError
        If `method` is not one of the allowed values.
    ValueError
        If required data is missing from `traverse_adata` (e.g., missing layers).
    KeyError
        If required columns are missing from `traverse_adata.obs`.

    Notes
    -----
    The function performs the following steps:
    1. Calculates differential effects using the specified method:
       - "max_possible": Direct difference between effect and control conditions
       - "min_possible": Normalized difference that considers the maximum possible
         effect from other dimensions
    2. Identifies maximum effects in positive and negative directions
    3. Stores results in the AnnData object for further analysis

    **Method Details:**

    - **max_possible**: Computes the direct log-fold-change between effect and control
      conditions. This represents the maximum possible effect a latent dimension can
      have on gene expression, assuming no constraints from other dimensions.

    - **min_possible**: More conservative approach that normalizes the effect by
      considering the maximum possible effects from other dimensions:
      ```
      normalized_effect = log(exp(effect) + pseudo_count + baseline) - baseline
      ```
      where baseline is the maximum possible effect from other dimensions.

    Examples
    --------
    >>> # Using max_possible method (default)
    >>> find_differential_effects(traverse_adata, method="max_possible")
    >>>
    >>> # Using min_possible method with custom parameters
    >>> find_differential_effects(
    ...     traverse_adata,
    ...     method="min_possible",
    ...     key_added="conservative",
    ...     add_to_counts=0.05,
    ...     relax_max_by=0.1,
    ... )
    >>> # Access results
    >>> stepwise_effects = traverse_adata.uns["effect_traverse_effect_stepwise"]
    >>> positive_effects = traverse_adata.varm["effect_traverse_effect_pos"]
    >>> negative_effects = traverse_adata.varm["effect_traverse_effect_neg"]
    """
    assert method in ["max_possible", "min_possible"]

    # Reorder the traverse_adata to original order
    original_traverse_adata = traverse_adata
    traverse_adata = traverse_adata[
        traverse_adata.obs.sort_values(["original_order"]).index
    ].copy()
    traverse_adata = traverse_adata[
        :, traverse_adata.var.sort_values(["original_order"]).index
    ].copy()

    # Get the number of latent dimensions, steps, samples, and vars
    n_latent, n_steps, n_samples, n_vars = get_dimensions_of_traverse_data(traverse_adata)

    # Get the dim_id and span values
    span_values = traverse_adata.obs["span_value"].values.reshape(n_latent, n_steps, n_samples)
    assert np.allclose(span_values, span_values.max(axis=-1, keepdims=True))
    span_values = span_values[:, :, 0]  # n_latent x n_steps
    dim_ids = traverse_adata.obs["dim_id"].values.reshape(n_latent, n_steps, n_samples)[
        :, 0, 0
    ]  # n_latent

    # Get the output mean parameters in 4D format
    control_mean_param = traverse_adata.layers["control"].reshape(
        n_latent, n_steps, n_samples, n_vars
    )
    effect_mean_param = traverse_adata.layers["effect"].reshape(
        n_latent, n_steps, n_samples, n_vars
    )

    # Helper functions
    average_over_samples = lambda x: x.mean(axis=2)
    add_eps_in_count_space = lambda x: scipy.special.logsumexp(
        np.stack([x, np.log(add_to_counts) * np.ones_like(x)]), axis=0
    )
    find_relative_effect = lambda x, baseline: (
        scipy.special.logsumexp(
            np.stack([x, np.log(add_to_counts) * np.ones_like(x), baseline]), axis=0
        )
        - baseline
    )

    # Find DE for each sample and average over samples
    if method == "max_possible":
        diff_considering_small_values = average_over_samples(
            add_eps_in_count_space(effect_mean_param) - add_eps_in_count_space(control_mean_param)
        )  # n_latent x n_steps x n_vars
    elif method == "min_possible":
        reduce_dims = (
            1,
            2,
        )
        max_of_two = np.maximum(
            effect_mean_param.max(axis=reduce_dims, keepdims=True),
            control_mean_param.max(axis=reduce_dims, keepdims=True),
        )  # n_latent x 1 x n_samples, n_vars
        max_cumulative_possible_all = scipy.special.logsumexp(
            max_of_two, axis=0, keepdims=True
        )  # 1 x 1 x n_samples, n_vars
        max_cumulative_possible_other_dims = (
            np.log(np.exp(max_cumulative_possible_all) - np.exp(max_of_two)) - relax_max_by
        )  # n_latent x 1 x n_samples, n_vars
        max_cumulative_possible_other_dims = max_cumulative_possible_other_dims + np.zeros_like(
            effect_mean_param
        )  # n_latent x n_steps x n_samples, n_vars
        normalized_effect_mean_param = find_relative_effect(
            effect_mean_param, max_cumulative_possible_other_dims
        )  # n_latent x n_steps x n_samples, n_vars
        normalized_control_mean_param = find_relative_effect(
            control_mean_param, max_cumulative_possible_other_dims
        )  # n_latent x n_steps x n_samples, n_vars
        diff_considering_small_values = average_over_samples(
            normalized_effect_mean_param - normalized_control_mean_param
        )  # n_latent x n_steps x n_vars
    else:
        raise NotImplementedError()

    original_traverse_adata.uns[f"{key_added}_traverse_effect_stepwise"] = (
        diff_considering_small_values
    )

    # Find DE vars in positive and negative directions
    for effect_sign in ["pos", "neg"]:
        mask = (
            np.where(span_values >= 0, 1, 0)
            if effect_sign == "pos"
            else np.where(span_values <= 0, 1, 0)
        )
        max_effect = np.max(
            np.expand_dims(mask, axis=-1) * diff_considering_small_values, axis=1
        )  # n_latent x n_vars
        max_effect = pd.DataFrame(max_effect.T, index=traverse_adata.var_names, columns=dim_ids)
        original_traverse_adata.varm[f"{key_added}_traverse_effect_{effect_sign}"] = (
            max_effect.loc[original_traverse_adata.var_names].copy()
        )
        original_traverse_adata.uns[f"{key_added}_traverse_effect_{effect_sign}_dim_ids"] = (
            original_traverse_adata.varm[
                f"{key_added}_traverse_effect_{effect_sign}"
            ].columns.values
        )
        original_traverse_adata.varm[
            f"{key_added}_traverse_effect_{effect_sign}"
        ].columns = original_traverse_adata.varm[
            f"{key_added}_traverse_effect_{effect_sign}"
        ].columns.astype(str)


def combine_differential_effects(
    traverse_adata: AnnData,
    keys: list[str],
    key_added: str,
    combine_function: callable,
) -> None:
    """Combine differential effects from multiple analyses.

    This function combines differential effects calculated using different methods
    or parameters to create a unified score. It applies a custom combination
    function to merge the stepwise effects and recalculates the positive/negative
    direction effects.

    Parameters
    ----------
    traverse_adata
        AnnData object containing differential effects from `find_differential_effects`.
        Must have `.uns[f"{key}_traverse_effect_stepwise"]` for each key in `keys`.
    keys
        List of keys corresponding to existing differential effect analyses.
        Each key should have been used in a previous call to `find_differential_effects`.
    key_added
        Prefix for the keys added to store the combined results.
        Results will be stored with keys like `{key_added}_traverse_effect_stepwise`.
    combine_function
        Function to combine the stepwise effects. Should take multiple arrays
        (one for each key) as positional arguments and return a single combined
        array with the same shape as the input arrays.

    Returns
    -------
    None
        Combined results are stored in `traverse_adata`:
        - `.uns[f"{key_added}_traverse_effect_stepwise"]`: Combined stepwise effects
          (shape: n_latent × n_steps × n_vars)
        - `.varm[f"{key_added}_traverse_effect_pos"]`: Combined positive direction effects
          (DataFrame with genes as rows, dimensions as columns)
        - `.varm[f"{key_added}_traverse_effect_neg"]`: Combined negative direction effects
          (DataFrame with genes as rows, dimensions as columns)
        - `.uns[f"{key_added}_traverse_effect_pos_dim_ids"]`: Array of dimension IDs
          for positive effects
        - `.uns[f"{key_added}_traverse_effect_neg_dim_ids"]`: Array of dimension IDs
          for negative effects

    Raises
    ------
    KeyError
        If any key in `keys` is missing from `traverse_adata.uns`.
    ValueError
        If the combine_function returns an array with unexpected shape.

    Notes
    -----
    The function performs the following steps:
    1. Reorders the data to original order for consistent processing
    2. Extracts stepwise effects for each key in `keys`
    3. Applies the `combine_function` to merge the effects
    4. Recalculates positive and negative direction effects from the combined data
    5. Stores the results with the new `key_added` prefix

    The `combine_function` should be designed to work with the specific types of
    differential effects being combined. Common approaches include:
    - Taking the maximum/minimum of multiple analyses
    - Computing weighted averages
    - Applying logical operations (AND/OR) for binary effects

    Examples
    --------
    >>> # Combine max_possible and min_possible effects using maximum
    >>> def combine_max_min(max_effects, min_effects):
    ...     return np.maximum(max_effects, min_effects)
    >>> combine_differential_effects(
    ...     traverse_adata,
    ...     keys=["max_possible", "min_possible"],
    ...     key_added="combined",
    ...     combine_function=combine_max_min,
    ... )
    """
    # Reorder the traverse_adata to original order
    original_traverse_adata = traverse_adata
    traverse_adata = traverse_adata[
        traverse_adata.obs.sort_values(["original_order"]).index
    ].copy()
    traverse_adata = traverse_adata[
        :, traverse_adata.var.sort_values(["original_order"]).index
    ].copy()

    # Get the number of latent dimensions, steps, samples, and vars
    n_latent, n_steps, n_samples, n_vars = get_dimensions_of_traverse_data(traverse_adata)

    # Get the dim_id and span values
    span_values = traverse_adata.obs["span_value"].values.reshape(n_latent, n_steps, n_samples)
    assert np.allclose(span_values, span_values.max(axis=-1, keepdims=True))
    span_values = span_values[:, :, 0]  # n_latent x n_steps
    dim_ids = traverse_adata.obs["dim_id"].values.reshape(n_latent, n_steps, n_samples)[
        :, 0, 0
    ]  # n_latent

    # Combine effects
    combined_traverse_effect_stepwise = combine_function(
        *[traverse_adata.uns[f"{key}_traverse_effect_stepwise"] for key in keys]
    )

    original_traverse_adata.uns[f"{key_added}_traverse_effect_stepwise"] = (
        combined_traverse_effect_stepwise
    )

    # Find DE vars in positive and negative directions
    for effect_sign in ["pos", "neg"]:
        mask = (
            np.where(span_values >= 0, 1, 0)
            if effect_sign == "pos"
            else np.where(span_values <= 0, 1, 0)
        )
        max_effect = np.max(
            np.expand_dims(mask, axis=-1) * combined_traverse_effect_stepwise, axis=1
        )  # n_latent x n_vars
        max_effect = pd.DataFrame(max_effect.T, index=traverse_adata.var_names, columns=dim_ids)
        original_traverse_adata.varm[f"{key_added}_traverse_effect_{effect_sign}"] = (
            max_effect.loc[original_traverse_adata.var_names].copy()
        )
        original_traverse_adata.uns[f"{key_added}_traverse_effect_{effect_sign}_dim_ids"] = (
            original_traverse_adata.varm[
                f"{key_added}_traverse_effect_{effect_sign}"
            ].columns.values
        )
        original_traverse_adata.varm[
            f"{key_added}_traverse_effect_{effect_sign}"
        ].columns = original_traverse_adata.varm[
            f"{key_added}_traverse_effect_{effect_sign}"
        ].columns.astype(str)


def calculate_differential_vars(traverse_adata: AnnData, **kwargs) -> None:
    """Calculate differential variables based on a combination of max_possible and min_possible effects.

    This function performs a comprehensive differential variable analysis by
    calculating effects using both "max_possible" and "min_possible" methods,
    then combining them into a unified score that considers both approaches.

    Parameters
    ----------
    traverse_adata
        AnnData object created by `traverse_latent` or `make_traverse_adata`.
        Must contain `.layers['control']` and `.layers['effect']`.
    **kwargs
        Additional keyword arguments passed to `find_differential_effects`.

    Returns
    -------
    None
        Results are stored in `traverse_adata`:
        - `.uns["max_possible_traverse_effect_stepwise"]`: Max possible stepwise effects
        - `.uns["min_possible_traverse_effect_stepwise"]`: Min possible stepwise effects
        - `.uns["combined_score_traverse_effect_stepwise"]`: Combined stepwise effects
        - `.varm["max_possible_traverse_effect_pos/neg"]`: Max possible direction effects
        - `.varm["min_possible_traverse_effect_pos/neg"]`: Min possible direction effects
        - `.varm["combined_score_traverse_effect_pos/neg"]`: Combined direction effects

    Notes
    -----
    The function performs the following steps:
    1. Calculates differential effects using "max_possible" method
    2. Calculates differential effects using "min_possible" method
    3. Combines the results using a custom scoring function

    **Scoring Logic:**

    The combined score is designed to identify genes that show consistent
    differential effects across both calculation methods. The filtering criteria are:

    1. **Base threshold**: max_possible >= 1.0 (ensures substantial effect)
    2. **Relative thresholds**: Either:
       - max_possible > 50% of its maximum across all dimensions and steps, OR
       - min_possible > 10% of its maximum across all dimensions and steps
    3. **Final score**: For genes passing the filters, score = max_possible × min_possible

    This approach provides a robust way to identify genes that show consistent
    differential effects across both calculation methods, reducing false positives
    while maintaining sensitivity to real biological effects.

    **Biological Interpretation:**

    - **High combined scores**: Genes with strong specific effects
    - **High max_possible, low min_possible**: These genes are usually shared with another program,
      with the other program having a larger effect.
    - **Low max_possible, high min_possible**: Not possible.
    - **Low scores in both**: Genes with weak or inconsistent effects

    Examples
    --------
    >>> # Basic differential variable calculation
    >>> calculate_differential_vars(traverse_adata)
    >>> # With custom parameters
    >>> calculate_differential_vars(traverse_adata, add_to_counts=0.05, relax_max_by=0.1)
    >>> # Access results
    >>> max_effects = traverse_adata.varm["max_possible_traverse_effect_pos"]
    >>> min_effects = traverse_adata.varm["min_possible_traverse_effect_pos"]
    >>> combined_effects = traverse_adata.varm["combined_score_traverse_effect_pos"]
    """
    print("Finding differential variables per latent dimension ...")
    find_differential_effects(
        traverse_adata, method="max_possible", key_added="max_possible", **kwargs
    )
    find_differential_effects(
        traverse_adata, method="min_possible", key_added="min_possible", **kwargs
    )

    def combine_function(min_possible, max_possible):
        # min_possible and max_possible dimensions: n_latent x n_steps x n_vars
        keep = (max_possible >= 1.0) & (
            (max_possible > max_possible.max(axis=(0, 1), keepdims=True) / 2)
            | (min_possible > min_possible.max(axis=(0, 1), keepdims=True) / 10)
        )
        score = np.where(keep, max_possible * min_possible, 0)
        return score

    # Combine scores with product
    combine_differential_effects(
        traverse_adata,
        keys=["min_possible", "max_possible"],
        key_added="combined_score",
        combine_function=combine_function,
    )


def get_split_effects(
    model: DRVI,
    embed: AnnData,
    n_steps: int = 20,
    n_samples: int = 100,
    traverse_kwargs: dict | None = None,
    de_kwargs: dict | None = None,
) -> AnnData:
    """Get split effects by performing latent traversal and differential analysis.

    This is a high-level function that combines latent space traversal with
    differential variable analysis. It performs the complete pipeline from
    generating traversal data to calculating differential effects.

    Parameters
    ----------
    model
        Trained DRVI model for decoding latent representations.
    embed
        AnnData object containing latent dimension statistics in `.var`.
        Must have columns: `original_dim_id`, `min`, `max`, `std`, `title`, `vanished`, `order`.
    n_steps
        Number of steps in the traversal. Must be even (half negative, half positive).
    n_samples
        Number of samples to generate for each step.
    traverse_kwargs
        Additional arguments passed to `traverse_latent`. Common options include:
        - `copy_adata_var_info`: Whether to copy variable information
        - `noise_formula`: Custom noise generation function
        - `max_noise_std`: Maximum noise standard deviation.
    de_kwargs
        Additional arguments passed to `calculate_differential_vars`. Common options include:
        - `add_to_counts`: Pseudo-count for log calculations
        - `relax_max_by`: Relaxation factor for min_possible method.

    Returns
    -------
    AnnData
        AnnData object containing both traversal data and differential analysis results.
        Includes all outputs from `traverse_latent` and `calculate_differential_vars`:

        **Traversal Data:**
        - `.X`: Difference between effect and control conditions
        - `.layers['control']`: Control condition gene expression
        - `.layers['effect']`: Effect condition gene expression
        - `.obs`: Metadata including dim_id, sample_id, step_id, span_value, title, vanished, order

        **Differential Analysis Results:**
        - `.uns["max_possible_traverse_effect_stepwise"]`: Max possible stepwise effects
        - `.uns["min_possible_traverse_effect_stepwise"]`: Min possible stepwise effects
        - `.uns["combined_score_traverse_effect_stepwise"]`: Combined stepwise effects
        - `.varm["*_traverse_effect_pos/neg"]`: Direction-specific effects for all methods

    Raises
    ------
    ValueError
        If required columns are missing from `embed.var`.
    KeyError
        If required data is missing from the model or embed objects.

    Notes
    -----
    This function is a convenience wrapper that performs the complete analysis
    pipeline in one call. It's equivalent to:

    ```python
    traverse_adata = traverse_latent(model, embed, n_steps, n_samples, **traverse_kwargs)
    calculate_differential_vars(traverse_adata, **de_kwargs)
    return traverse_adata
    ```

    **Functionality:**

    1. **Traversal Generation**: Creates systematic traversals through latent space
    2. **Decoding**: Converts latent traversals to gene expression predictions
    3. **Differential Analysis**: Calculates effects using max_possible and min_possible methods
    4. **Combination**: Creates unified scores for robust gene identification

    **Use Cases:**

    - **Gene Discovery**: Identify genes associated with specific latent dimensions
    - **Pathway Analysis**: Understand biological processes captured by latent factors
    - **Model Validation**: Verify that latent dimensions have interpretable biological meaning
    - **Comparative Analysis**: Compare effects across different models or conditions

    Examples
    --------
    >>> # Basic split effects analysis
    >>> split_effects = get_split_effects(model, embed)
    >>> # With custom parameters
    >>> split_effects = get_split_effects(model, embed, n_steps=30, n_samples=50)
    >>> # Access results
    >>> combined_effects = split_effects.varm["combined_score_traverse_effect_pos"]
    >>> stepwise_effects = split_effects.uns["combined_score_traverse_effect_stepwise"]
    """
    if traverse_kwargs is None:
        traverse_kwargs = {}
    if de_kwargs is None:
        de_kwargs = {}

    traverse_adata = traverse_latent(
        model, embed, n_steps=n_steps, n_samples=n_samples, **traverse_kwargs
    )
    calculate_differential_vars(traverse_adata, **de_kwargs)

    return traverse_adata


def iterate_on_top_differential_vars(
    traverse_adata: AnnData,
    key: str,
    title_col: str = "title",
    order_col: str = "order",
    gene_symbols: str | None = None,
    score_threshold: float = 0.0,
) -> list[tuple[str, pd.Series]]:
    """Create an iterator of top differential variables per latent dimension.

    This function processes differential analysis results to create an organized
    list of top differentially expressed genes for each latent dimension,
    sorted by their effect scores and organized by dimension.

    Parameters
    ----------
    traverse_adata
        AnnData object with differential analysis results from `calculate_differential_vars`.
        Must contain differential effect data for the specified `key`.
    key
        Key prefix for the differential variables in `traverse_adata`.
        Should correspond to a key used in `find_differential_effects` or `calculate_differential_vars`.
        Common value: "combined_score".
    title_col
        Column name in `traverse_adata.obs` containing dimension titles.
        These titles will be used in the output dimension names.
    order_col
        Column name in `traverse_adata.obs` containing dimension ordering.
    gene_symbols
        Column name in `traverse_adata.var` containing gene symbols.
        If None, uses the index of `traverse_adata.var` (usually gene IDs).
        Useful for converting between gene IDs and readable gene names.
    score_threshold
        Minimum score threshold to include genes in the results.
        Only genes with scores above this threshold will be included.

    Returns
    -------
    list[tuple[str, pd.Series]]
        List of tuples, where each tuple contains:
        - str: Dimension title with direction indicator (e.g., "Cell Cycle+", "Cell Cycle-")
        - pd.Series: Series of gene scores for that dimension/direction, sorted descending

        The list is sorted by dimension order, with each dimension appearing at most twice
        (once for positive effects, once for negative effects).

    Raises
    ------
    KeyError
        If required columns or differential effect data are missing.
    ValueError
        If the specified key doesn't exist in the AnnData object.

    Notes
    -----
    The function performs the following steps:
    1. Extracts positive and negative differential effects for the specified key
    2. Maps gene names to symbols if `gene_symbols` is provided
    3. Filters genes by score threshold
    4. Organizes results by dimension and direction (positive/negative)
    5. Returns a list sorted by dimension order

    **Output Structure:**

    Each dimension appears twice in the results - once for positive effects
    and once for negative effects. The direction is indicated by "+" or "-"
    appended to the dimension title.

    Only dimensions with at least one gene above the threshold are included.

    Examples
    --------
    >>> # Basic iteration over top differential variables
    >>> top_vars = iterate_on_top_differential_vars(traverse_adata, "combined_score")
    >>> for dim_title, gene_scores in top_vars:
    ...     print(f"{dim_title}: {len(gene_scores)} genes")
    ...     print(f"Top genes: {gene_scores.head().index.tolist()}")
    >>> # With custom parameters and gene symbols
    >>> top_vars = iterate_on_top_differential_vars(
    ...     traverse_adata, "max_possible", gene_symbols="gene_symbol", score_threshold=1.0
    ... )
    >>> # Create a summary of results
    >>> for dim_title, gene_scores in top_vars:
    ...     print(f"{dim_title}: {gene_scores.head().index.tolist()}")
    """
    df_pos = traverse_adata.varm[f"{key}_traverse_effect_pos"].copy()
    df_neg = traverse_adata.varm[f"{key}_traverse_effect_neg"].copy()

    if gene_symbols is not None:
        gene_name_mapping = dict(
            zip(traverse_adata.var.index, traverse_adata.var[gene_symbols], strict=False)
        )
    else:
        gene_name_mapping = dict(
            zip(traverse_adata.var.index, traverse_adata.var.index, strict=False)
        )
    df_pos.index = df_pos.index.map(gene_name_mapping)
    df_neg.index = df_neg.index.map(gene_name_mapping)

    de_info = dict(
        **{
            (str(k) + "+"): v.sort_values(ascending=False)
            .where(lambda x: x > score_threshold)
            .dropna()
            for k, v in df_pos.to_dict(orient="series").items()
        },
        **{
            (str(k) + "-"): v.sort_values(ascending=False)
            .where(lambda x: x > score_threshold)
            .dropna()
            for k, v in df_neg.to_dict(orient="series").items()
        },
    )

    return [
        (f"{row[title_col]}{direction}", de_info[f"{row['dim_id']}{direction}"])
        for i, row in traverse_adata.obs[["dim_id", order_col, title_col]]
        .drop_duplicates()
        .sort_values(order_col)
        .iterrows()
        for direction in ["-", "+"]
        if len(de_info[f"{row['dim_id']}{direction}"]) > 0
    ]
