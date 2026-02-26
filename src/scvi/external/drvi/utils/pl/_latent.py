from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import scanpy as sc
from matplotlib import pyplot as plt

from scvi.external.drvi.utils.pl import cmap

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

    from anndata import AnnData


def make_balanced_subsample(adata: AnnData, col: str, min_count: int = 10) -> AnnData:
    """Create a balanced subsample of AnnData based on a categorical column.

    This function creates a balanced subsample by sampling an equal number of cells
    from each category in the specified column, ensuring balanced representation.

    Parameters
    ----------
    adata
        Annotated data object to subsample.
    col
        Column name in `adata.obs` containing categorical labels for balancing.
    min_count
        Minimum number of samples per category. If a category has fewer samples
        than this, sampling will be done with replacement.

    Returns
    -------
    AnnData
        Balanced subsample of the input AnnData object.

    Notes
    -----
    The function uses a fixed random state (0) for reproducible results.
    If a category has fewer samples than `min_count`, sampling is done with replacement.
    """
    n_sample_per_cond = adata.obs[col].value_counts().min()
    balanced_sample_index = (
        adata.obs.groupby(col)
        .sample(
            n=max(min_count, n_sample_per_cond),
            random_state=0,
            replace=n_sample_per_cond < min_count,
        )
        .index
    )
    adata = adata[balanced_sample_index].copy()
    return adata


def plot_latent_dimension_stats(
    embed: AnnData,
    figsize: tuple[int, int] = (5, 3),
    log_scale: bool | Literal["try"] = "try",
    ncols: int = 5,
    columns: Sequence[str] = ("reconstruction_effect", "max_value", "mean", "std"),
    titles: dict[str, str] | None = None,
    remove_vanished: bool = False,
    show: bool = True,
):
    """Plot the statistics of latent dimensions.

    This function creates line plots showing various statistics of latent dimensions
    across their ranking order. It can optionally distinguish between vanished and
    non-vanished dimensions.

    Parameters
    ----------
    embed
        Annotated data object containing the latent dimensions and their statistics
        in the `.var` attribute.
    figsize
        The size of each subplot (width, height) in inches.
    log_scale
        Whether to use a log scale for the y-axis. If "try", log scale is used
        only if the minimum value is greater than 0.
    ncols
        The maximum number of columns in the subplot grid.
    columns
        The columns from `embed.var` to plot. These should be numeric columns
        containing dimension statistics.
    titles
        Custom titles for each column in the plot. If None, default titles are used.
    remove_vanished
        Whether to exclude vanished dimensions from the plot.
    show
        Whether to display the plot. If False, returns the figure object.

    Returns
    -------
    matplotlib.figure.Figure or None
        The matplotlib figure object if `show=False`, otherwise None.

    Notes
    -----
    The function expects the following columns in `embed.var`:
    - `order`: Ranking of dimensions
    - `vanished`: Boolean indicating vanished dimensions
    - The columns specified in the `columns` parameter

    If `remove_vanished=False`, a legend is added to distinguish between
    vanished (black dots) and non-vanished (blue dots) dimensions.

    Examples
    --------
    >>> # Default plot
    >>> plot_latent_dimension_stats(embed)
    >>>
    >>> # Plot basic statistics
    >>> plot_latent_dimension_stats(embed, columns=["reconstruction_effect", "max_value"])
    >>> # Plot with custom titles and log scale
    >>> titles = {"reconstruction_effect": "Reconstruction Impact", "max_value": "Max Activation"}
    >>> plot_latent_dimension_stats(embed, titles=titles, log_scale=True)
    """
    if titles is None:
        titles = {
            "reconstruction_effect": "Reconstruction effect",
            "max_value": "Max value",
            "mean": "Mean",
            "std": "Standard Deviation",
        }
    nrows = int(np.ceil(len(columns) / ncols))
    if nrows == 1:
        ncols = len(columns)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(figsize[0] * ncols, figsize[1] * nrows),
        sharey=False,
        sharex=False,
        squeeze=False,
    )

    # Iterate through columns and plot the data
    for ax, col in zip(axes.flatten(), columns, strict=False):
        df = embed.var
        if remove_vanished:
            df = df.query("vanished == False")
        df = df.sort_values("order")
        ranks = df["order"]
        x = df[col]

        ax.plot(ranks, x, linestyle="-", color="grey", label="Line")  # Solid line plot
        for vanished_status_to_plot in [True, False]:
            indices = df["vanished"] == vanished_status_to_plot
            ax.plot(
                ranks[indices],
                x[indices],
                "o",
                markersize=3,
                color="black" if vanished_status_to_plot else "blue",
                label="Data Points",
            )

        # Adding labels and title
        ax.set_xlabel("Rank based on Explanation Share")
        ax.set_ylabel(titles[col] if col in titles else col)
        if isinstance(log_scale, str):
            if log_scale == "try":
                if x.min() > 0:
                    ax.set_yscale("log")
        else:
            if log_scale:
                ax.set_yscale("log")

        # Removing the legend
        ax.legend().remove()

        # Adding grid
        ax.grid(axis="x")

    if not remove_vanished:
        # Create custom legend entries
        handles = []
        for vanished_status_to_plot in [False, True]:
            color = "black" if vanished_status_to_plot else "blue"
            label = "Vanished" if vanished_status_to_plot else "Non-vanished"
            handles.append(
                plt.Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    markerfacecolor=color,
                    markersize=5,
                    label=label,
                )
            )

        # Add the legend to the first subplot or the entire figure
        fig.legend(
            handles=handles,
            labels=[handle.get_label() for handle in handles],
            loc="center left",
            bbox_to_anchor=(1, 0.5),
            title=None,
        )

    for ax in axes.flatten()[len(columns) :]:
        fig.delaxes(ax)

    plt.tight_layout()
    if show:
        plt.show()
    else:
        return fig


def plot_latent_dims_in_heatmap(
    embed: AnnData,
    categorical_column: str,
    title_col: str | None = "title",
    sort_by_categorical: bool = False,
    make_balanced: bool = True,
    order_col: str | None = "order",
    remove_vanished: bool = True,
    figsize: tuple[int, int] | None = None,
    show: bool = True,
    **kwargs,
):
    """Plot the latent dimensions in a heatmap.

    This function creates a heatmap showing the values of latent dimensions
    across different categories. It can optionally create balanced subsamples
    and sort dimensions based on categorical differences.

    Parameters
    ----------
    embed
        Annotated data object containing the latent dimensions in `.X`
        and categorical metadata in `.obs`.
    categorical_column
        The column in `embed.obs` that represents the categorical variable
        for grouping cells.
    title_col
        The column in `embed.var` to use as titles for each dimension.
        If None, uses the dimension indices.
    sort_by_categorical
        Whether to sort dimensions based on their maximum absolute values
        within each category. If True, `order_col` is ignored.
    make_balanced
        Whether to create a balanced subsample of the data based on the
        categorical variable using `make_balanced_subsample`.
    order_col
        The column in `embed.var` to use for ordering the dimensions.
        Ignored if `sort_by_categorical=True`.
    remove_vanished
        Whether to remove vanished dimensions from the plot.
    figsize
        The size of the figure (width, height) in inches.
        If None, automatically calculated based on number of categories.
    show
        Whether to display the plot. If False, returns the plot object.
    **kwargs
        Additional keyword arguments passed to `sc.pl.heatmap`.

    Returns
    -------
    matplotlib.axes.Axes or None
        The heatmap axes if `show=False`, otherwise None.

    Raises
    ------
    ValueError
        If required columns (`order_col` or "vanished") are not found in `embed.var`.

    Notes
    -----
    The function expects the following columns in `embed.var`:
    - `order_col`: For ordering dimensions (if `sort_by_categorical=False`)
    - `title_col`: For dimension titles
    - `vanished`: Boolean indicating vanished dimensions (if `remove_vanished=True`)

    If `figsize=None`, the figure height is automatically calculated as
    `len(unique_categories) / 6` to accommodate all categories.

    The heatmap uses a red-blue color map centered at 0, with no dendrogram.

    Examples
    --------
    >>> # Basic heatmap of latent dimensions by cell type
    >>> plot_latent_dims_in_heatmap(embed, categorical_column="cell_type")
    >>> # Heatmap with balanced sampling and custom sorting
    >>> plot_latent_dims_in_heatmap(
    ...     embed, categorical_column="condition", sort_by_categorical=True, make_balanced=True
    ... )
    >>> # Heatmap with custom figure size
    >>> plot_latent_dims_in_heatmap(embed, categorical_column="batch", figsize=(12, 8))
    """
    if order_col is not None and order_col not in embed.var:
        raise ValueError(
            f'Column "{order_col}" not found in `embed.var`. '
            "Please run `set_latent_dimension_stats` to set order."
        )
    if remove_vanished:
        if "vanished" not in embed.var:
            raise ValueError(
                'Column "vanished" not found in `embed.var`. '
                "Please run `set_latent_dimension_stats` to set vanished status."
            )
        embed = embed[:, ~embed.var["vanished"]]

    if make_balanced:
        embed = make_balanced_subsample(embed, categorical_column)

    if sort_by_categorical:
        dim_order = np.abs(embed.X).argmax(axis=0).argsort().tolist()
    elif order_col is None:
        dim_order = np.arange(embed.n_vars)
    else:
        dim_order = embed.var[order_col].argsort().tolist()

    if title_col is None:
        vars_to_show = embed.var.iloc[dim_order].index
    else:
        vars_to_show = embed.var.iloc[dim_order][title_col]

    if figsize is None:
        figsize = (10, len(embed.obs[categorical_column].unique()) / 6)

    kwargs = {
        **dict(  # noqa: C408
            vcenter=0,
            cmap=cmap.saturated_red_blue_cmap,
            dendrogram=False,
        ),
        **kwargs,
    }

    return sc.pl.heatmap(
        embed,
        vars_to_show,
        categorical_column,
        gene_symbols=title_col,
        figsize=figsize,
        show_gene_labels=True,
        show=show,
        **kwargs,
    )
