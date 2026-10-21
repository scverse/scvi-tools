from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from scvi.utils import dependencies

if TYPE_CHECKING:
    import xarray as xr


@dependencies("plotnine")
def plot_hier_importance(
    gene_results: xr.Dataset,
    feature: int | str = 0,
    base_resolution: int | None = None,
    significance_threshold: float = 0.1,
    color_by: str | None = None,
    plot_fig: bool = True,
    theme_kwargs: dict | None = None,
):
    """Visualize hierarchical CRT results as a multi-resolution significance dendrogram.

    Parameters
    ----------
    gene_results
        Output of :meth:`~scvi.external.VIVS.get_hier_importance`.
    feature
        Which response/feature to plot (index into ``gene_results.feature``, or one of
        its values). Required because significance is inherently per-response; call this
        once per feature of interest rather than aggregating across responses.
    base_resolution
        Resolution used to pick which genes to plot (defaults to the coarsest resolution).
    significance_threshold
        BH-adjusted p-value threshold for calling a gene/cluster significant.
    color_by
        Optional data variable in ``gene_results`` to color an extra annotation row by.
    plot_fig
        If ``True``, returns a ``plotnine`` figure; otherwise returns the raw
        ``(plot_df, labels, breaks)`` used to build it.
    theme_kwargs
        Extra kwargs forwarded to ``plotnine.theme()``.
    """
    if isinstance(feature, int):
        gene_results = gene_results.isel(feature=feature)
    else:
        gene_results = gene_results.sel(feature=feature)

    if base_resolution is None:
        base_resolution = gene_results.resolution.values.min()
    padjs = gene_results.loc[{"resolution": base_resolution}]["padj"].to_pandas()
    genes_to_plot = padjs.loc[lambda x: x < significance_threshold].index.tolist()
    res_subset = gene_results.loc[{"gene_name": genes_to_plot}].assign(
        gene_index=("gene_name", np.arange(len(genes_to_plot)))
    )
    plot_df = []
    for resolution_idx, resolution in enumerate(res_subset.resolution.values):
        res_ = res_subset.loc[{"resolution": resolution}]
        unique_clusters = np.unique(res_["cluster_assignment"].values)
        for cluster in unique_clusters:
            gene_is_in_cluster = res_["cluster_assignment"] == cluster
            res_cluster = res_.loc[{"gene_name": gene_is_in_cluster}]
            xmin = res_cluster.gene_index.values.min()
            xmax = res_cluster.gene_index.values.max()
            are_indices_contiguous = xmax - xmin + 1 == res_cluster.gene_name.shape[0]
            if not are_indices_contiguous:
                raise ValueError("Gene indices are not contiguous")
            is_cluster_detected = (res_cluster["padj"] < significance_threshold).all()
            if is_cluster_detected:
                plot_df.append(
                    {
                        "resolution_idx": resolution_idx,
                        "resolution": resolution,
                        "xmin": xmin - 0.5,
                        "xmax": xmax + 0.5,
                        "ymin": resolution_idx,
                        "ymax": resolution_idx + 1,
                    }
                )
    plot_df = pd.DataFrame(
        plot_df, columns=["resolution_idx", "resolution", "xmin", "xmax", "ymin", "ymax"]
    )

    plot_df_color = None
    if color_by is not None:
        plot_df_color = []
        for gene in res_subset.gene_name.values:
            res_cluster = res_subset.loc[{"gene_name": gene}]
            xmin = res_cluster.gene_index.values.min()
            xmax = res_cluster.gene_index.values.max()
            plot_df_color.append(
                {
                    "xmin": xmin - 0.5,
                    "xmax": xmax + 0.5,
                    "ymin": -2,
                    "ymax": -1,
                    "color": res_cluster[color_by].values.item(),
                }
            )
        plot_df_color = pd.DataFrame(plot_df_color)

    labels = list(res_subset.gene_name.values)
    breaks = list(res_subset.gene_index.values)
    if not plot_fig:
        return plot_df, labels, breaks

    import plotnine as p9

    theme_kwargs = theme_kwargs if theme_kwargs is not None else {"figure_size": (15, 2)}
    fig = (
        p9.ggplot(plot_df)
        + p9.geom_rect(
            plot_df,
            p9.aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"),
            inherit_aes=False,
            fill="#ededed",
            color="#080808",
        )
        + p9.scale_x_continuous(labels=labels, breaks=breaks)
        + p9.theme_classic()
        + p9.theme(
            axis_text_x=p9.element_text(rotation=90),
            axis_line_y=p9.element_blank(),
            axis_text_y=p9.element_blank(),
            axis_ticks_major_y=p9.element_blank(),
        )
        + p9.labs(x="", y="")
        + p9.theme(**theme_kwargs)
    )
    if plot_df_color is not None:
        fig = fig + p9.geom_rect(
            plot_df_color,
            p9.aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill="color"),
        )
    return fig
