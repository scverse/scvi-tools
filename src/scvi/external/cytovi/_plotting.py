import math

import anndata as ad
from anndata import AnnData

from scvi.utils import dependencies

from ._preprocessing import subsample
from ._utils import validate_layer_key, validate_marker, validate_obs_keys


@dependencies("seaborn", "matplotlib")
def plot_histogram(
    adata: ad.AnnData,
    marker: str | list[str] = "all",
    groupby: str = None,
    layer_key: str = "raw",
    downsample: bool = True,
    n_obs: int = 10000,
    col_wrap: int = None,
    tight_layout: bool = True,
    save: bool | str = None,
    return_plot: bool = False,
    kde_kwargs: dict = None,
    **kwargs,
):
    """
    Create a FacetGrid of histograms for specified markers in an AnnData object.

    Parameters
    ----------
    adata
        Annotated data matrix.
    marker
        Marker name or list of marker names to plot. Use ``'all'`` to plot all markers.
        Default: ``'all'`` uses all variables in ``adata.var_names``.
    groupby
        Key in ``adata.obs`` used to group/color the data (e.g., batch or condition).
        Default: no grouping.
    layer_key
        Layer key in ``adata.layers`` to draw values from. Default: use ``X``.
    downsample
        Whether to downsample when the number of observations is large. Default: True.
    n_obs
        Number of observations to keep if ``downsample`` is True. Default: 10_000.
    col_wrap
        Maximum number of columns before wrapping to a new row. If ``None``, chosen
        automatically based on the number of markers. Default: None.
    tight_layout
        Apply a tight layout to the figure. Default: True.
    save
        If ``True``, save as ``"marker_histogram.png"``; if a string is provided,
        use it as the filename. Default: False.
    return_plot
        If ``True``, return the ``seaborn.FacetGrid`` object. Default: False.
    kde_kwargs
        Additional keyword arguments passed to ``seaborn.kdeplot`` if a density
        overlay is enabled. Default: ``{}``.
    **kwargs
        Additional keyword arguments forwarded to ``seaborn.FacetGrid``.

    Returns
    -------
    If ``return_plot`` is ``True``, returns the ``seaborn.FacetGrid``; otherwise, nothing.

    Examples
    --------
    >>> cytovi.plot_histogram(adata, marker=["CD3", "CD4"], groupby="condition")
    >>> cytovi.plot_histogram(adata, marker="all", groupby="batch")
    """
    import seaborn as sns

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

    data_plot_melt = data_plot.melt(id_vars=groupby, var_name="variable", value_name="value")

    # generate the plot
    g = sns.FacetGrid(
        data_plot_melt,
        col="variable",
        hue=groupby,
        col_wrap=col_wrap,
        sharey=False,
        sharex=False,
        **kwargs,
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


@dependencies("seaborn", "matplotlib")
def plot_biaxial(
    adata: AnnData,
    marker_x: str | list[str] = None,
    marker_y: str | list[str] = None,
    color: str = None,
    n_bins: int = 10,
    layer_key: str = "raw",
    downsample: bool = True,
    n_obs: int = 10000,
    sample_color_groups: bool = False,
    save: bool | str = None,
    kde: bool = True,
    kde_kwargs: dict = None,
    scatter_kwargs: dict = None,
    **kwargs,
):
    """
    Create a PairGrid of biaxial plots for specified markers in an AnnData object.

    Parameters
    ----------
    adata
        Annotated data matrix.
    marker_x
        Variable name(s) to plot on the x-axis. A string or a list of marker names.
    marker_y
        Variable name(s) to plot on the y-axis. A string or a list of marker names.
    color
        Key in ``adata.obs`` used to color points (e.g., batch or condition).
        Default: no coloring.
    n_bins
        Number of contour levels for the KDE density. Default: 10.
    layer_key
        Layer key in ``adata.layers`` to draw values from. Default: use ``X``.
    downsample
        Whether to downsample when the number of observations is large. Default: True.
    n_obs
        Number of observations to keep if ``downsample`` is True. Default: 10_000.
    sample_color_groups
        If True and ``color`` is set, sample within each color group. Default: False.
    save
        If ``True``, save as ``"marker_biaxial.png"``;
        if a string is provided, use it as the filename. Default: False.
    kde
        Whether to overlay KDE density contours. Default: True.
    kde_kwargs
        Additional keyword arguments forwarded to ``seaborn.kdeplot``.
    scatter_kwargs
        Additional keyword arguments forwarded to ``seaborn.scatterplot``.
    **kwargs
        Additional keyword arguments forwarded to ``seaborn.PairGrid``.

    Returns
    -------
    No return value; the figure is displayed and optionally saved.

    Examples
    --------
    >>> cytovi.pl.biaxial(adata, marker_x="CD3", marker_y="CD4", color="condition")
    >>> cytovi.pl.biaxial(adata, marker_x=["CD8", "CD20"], marker_y="CD56", color="batch")
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

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
