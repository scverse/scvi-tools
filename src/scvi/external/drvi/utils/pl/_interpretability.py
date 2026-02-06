from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from scvi.external.drvi.utils.tl import iterate_on_top_differential_vars

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    from anndata import AnnData


def make_heatmap_groups(ordered_list: list) -> tuple[list[tuple[int, int]], list[Any]]:
    """Create group positions and labels for scanpy heatmap visualization of marker genes.

    This helper function processes an ordered list to identify groups of
    consecutive identical elements and returns their positions and labels.
    It's used to create group annotations for scanpy heatmap plots.

    Parameters
    ----------
    ordered_list
        List of elements where consecutive identical elements form groups.

    Returns
    -------
    tuple[list[tuple[int, int]], list]
        A tuple containing:
        - List of tuples with (start_index, end_index) for each group
        - List of group labels (unique values from ordered_list)

    Notes
    -----
    The function uses `itertools.groupby` to identify consecutive groups
    of identical elements. Each group is represented by its start and end
    indices (inclusive).

    Examples
    --------
    >>> # Simple example
    >>> groups, labels = make_heatmap_groups(["A", "A", "B", "B", "B", "A"])
    >>> print(f"Groups: {groups}")  # [(0, 1), (2, 4), (5, 5)]
    >>> print(f"Labels: {labels}")  # ['A', 'B', 'A']
    """
    n_groups, group_names = zip(
        *[(len(list(group)), key) for (key, group) in itertools.groupby(ordered_list)],
        strict=False,
    )
    group_positions = [0] + list(itertools.accumulate(n_groups))
    group_positions = list(
        zip(group_positions[:-1], [c - 1 for c in group_positions[1:]], strict=False)
    )
    return group_positions, group_names


def _bar_plot_top_differential_vars(
    plot_info: Sequence[tuple[str, pd.Series]],
    dim_subset: Sequence[str] | None = None,
    n_top_genes: int = 10,
    ncols: int = 5,
    show: bool = True,
):
    """Plot the top differential variables in a group of bar plots.

    This internal function creates horizontal bar plots showing the top genes
    for each latent dimension based on their differential effect scores.

    Parameters
    ----------
    plot_info
        Sequence of tuples containing dimension titles and corresponding gene data.
    dim_subset
        Subset of dimensions to plot. If None, all dimensions are plotted.
    n_top_genes
        Number of top genes to show in each plot.
    ncols
        Number of columns in the subplot grid.
    show
        Whether to display the plot. If False, returns the figure object.

    Returns
    -------
    matplotlib.figure.Figure or None
        The figure object if `show=False`, otherwise None.

    Notes
    -----
    The function creates a grid of horizontal bar plots, with each subplot
    showing the top genes for one latent dimension. Genes are sorted by
    their effect scores in descending order.

    **Plot Features:**

    - **Horizontal bars**: Gene names on y-axis, scores on x-axis
    - **Color**: Sky blue bars for all genes
    - **Grid**: No grid lines for cleaner appearance
    - **Layout**: Automatic grid layout based on number of dimensions

    Examples
    --------
    >>> # Basic bar plot
    >>> _bar_plot_top_differential_vars(plot_info)
    >>> # Custom layout
    >>> _bar_plot_top_differential_vars(plot_info, n_top_genes=15, ncols=3, show=False)
    """
    if dim_subset is not None:
        plot_info = dict(plot_info)
        plot_info = [(dim_id, plot_info[dim_id]) for dim_id in dim_subset]

    n_row = int(np.ceil(len(plot_info) / ncols))
    fig, axes = plt.subplots(n_row, ncols, figsize=(3 * ncols, int(1 + 0.2 * n_top_genes) * n_row))

    for ax, info in zip(axes.flatten(), plot_info, strict=False):
        dim_title = info[0]

        top_indices = info[1].sort_values(ascending=False)[:n_top_genes]
        genes = top_indices.index
        values = top_indices.values

        # Create a horizontal bar plot
        ax.barh(genes, values, color="skyblue")
        ax.set_xlabel("Gene Score")
        ax.set_title(dim_title)
        ax.invert_yaxis()
        ax.grid(False)

    for ax in axes.flatten()[len(plot_info) :]:
        fig.delaxes(ax)

    plt.tight_layout()
    if show:
        plt.show()
    else:
        return fig


def show_top_differential_vars(
    traverse_adata: AnnData,
    key: str,
    title_col: str = "title",
    order_col: str = "order",
    dim_subset: Sequence[str] | None = None,
    gene_symbols: str | None = None,
    score_threshold: float = 0.0,
    n_top_genes: int = 10,
    ncols: int = 5,
    show: bool = True,
):
    """Show top differential variables in a bar plot.

    This function creates a comprehensive visualization of the top differentially
    expressed genes for each latent dimension. It generates horizontal bar plots
    showing the genes with the highest effect scores for each dimension.

    Parameters
    ----------
    traverse_adata
        AnnData object containing the differential analysis results from
        `calculate_differential_vars`. Must contain differential effect data
        for the specified key.
    key
        Key prefix for the differential variables in `traverse_adata.varm`.
        Should correspond to a key used in `find_differential_effects` or
        `calculate_differential_vars` (e.g., "max_possible", "min_possible", "combined_score").
    title_col
        Column name in `traverse_adata.obs` that contains the titles for each dimension.
        These titles will be used as subplot titles.
    order_col
        Column name in `traverse_adata.obs` that specifies the order of dimensions.
        Results will be sorted by this column. Ignored if `dim_subset` is provided.
    dim_subset
        List of dimensions to plot in the bar plot. If None, all dimensions
        with significant effects are plotted.
    gene_symbols
        Column name in `traverse_adata.var` that contains gene symbols.
        If provided, gene symbols will be used in the plot instead of gene indices.
        Useful for converting between gene IDs and readable gene names.
    score_threshold
        Threshold value for gene scores. Only genes with scores above this
        threshold will be plotted.
    n_top_genes
        Number of top genes to plot for each dimension.
    ncols
        Number of columns in the plot grid.
    show
        Whether to display the plot. If False, returns the figure object.

    Returns
    -------
    matplotlib.figure.Figure or None
        The figure object if `show=False`, otherwise None.

    Raises
    ------
    KeyError
        If required data is missing from `traverse_adata`.
    ValueError
        If the specified key doesn't exist in the AnnData object.

    Notes
    -----
    The function performs the following steps:
    1. Extracts top differential variables using `iterate_on_top_differential_vars`
    2. Filters dimensions based on `dim_subset` if provided
    3. Creates horizontal bar plots for each dimension
    4. Displays top `n_top_genes` genes sorted by their effect scores

    **Visualization Features:**

    - **Gene symbols**: If provided, gene symbols will be used instead of gene indices.
    - **Grid layout**: Automatic grid based on number of dimensions and `ncols`
    - **Horizontal bars**: Gene names on y-axis, scores on x-axis
    - **Color coding**: Sky blue bars for all genes
    - **Dimension titles**: Each subplot shows the dimension title
    - **Gene ordering**: Genes sorted by effect score (highest first)

    **Interpretation:**

    - **Bar length**: Represents the magnitude of the differential effect
    - **Gene position**: Higher bars indicate stronger effects
    - **Dimension separation**: Each subplot shows effects for one latent dimension
    - **Direction indicators**: Dimension titles include "+" or "-" to indicate effect direction

    Examples
    --------
    >>> # Basic visualization with combined scores
    >>> show_top_differential_vars(traverse_adata, "combined_score")
    >>> # Custom parameters with gene symbols
    >>> show_top_differential_vars(
    ...     traverse_adata,
    ...     "max_possible",
    ...     gene_symbols="gene_symbol",
    ...     score_threshold=1.0,
    ...     n_top_genes=15,
    ...     ncols=3,
    ... )
    >>> # Subset of dimensions
    >>> show_top_differential_vars(
    ...     traverse_adata, "combined_score", dim_subset=["DR 5+", "DR 12+", "DR 14+"]
    ... )
    """
    plot_info = iterate_on_top_differential_vars(
        traverse_adata, key, title_col, order_col, gene_symbols, score_threshold
    )

    return _bar_plot_top_differential_vars(plot_info, dim_subset, n_top_genes, ncols, show)
