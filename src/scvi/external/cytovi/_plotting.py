import math

import anndata as ad
from anndata import AnnData

from scvi.utils import dependencies

from ._preprocessing import subsample
from ._utils import validate_layer_key, validate_marker, validate_obs_keys


@dependencies("seaborn", "matplotlib")
def histogram(
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
        Number of columns to wrap the plots. If None, it is calculated based on the
        number of markers.

    tight_layout : bool, optional
        Whether to use tight layout for the plot.

    save : Union[bool, str], optional
        If True, the plot is saved as "marker_histogram.png". If a string is provided,
        the plot is saved with the given filename.

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
def biaxial(
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
        If True, save the plot as "marker_histogram.png". If a string is provided, save the
          plot with the given filename.

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
