import math

import anndata as ad

from scvi.utils import dependencies

from ._utils import validate_layer_key, validate_marker, validate_obs_keys, subsample


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
    """Create a FacetGrid of histograms for specified markers in an AnnData object.
    
    This function is adapted from CytoVI's plotting utilities.

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

    # Subsample if too many observations
    if downsample and adata.n_obs > 10000:
        adata = subsample(adata, n_obs=n_obs, groupby=groupby)

    num_plots = len(marker)

    if col_wrap is None:
        col_wrap = math.ceil(math.sqrt(num_plots))

    data_plot = adata[:, marker].to_df(layer=layer_key)

    if groupby is not None:
        data_plot[groupby] = adata.obs[groupby]

    data_plot_melt = data_plot.melt(id_vars=groupby, var_name="variable", value_name="value")

    # Generate the plot
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
