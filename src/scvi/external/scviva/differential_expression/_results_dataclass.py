from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Union

import numpy as np
import pandas as pd
from rich import print
from sklearn.gaussian_process import GaussianProcessClassifier

from scvi.utils import dependencies


@dataclass
class DifferentialExpressionResults:
    """Dataclass for storing the results of the differential expression analysis,
    including the GP classifier
    """  # noqa: D205

    gpc: GaussianProcessClassifier
    g1_g2: pd.DataFrame
    g1_n1: pd.DataFrame
    n1_g2: pd.DataFrame
    n1_n2: Union[pd.DataFrame, None] = field(default=None)  # noqa: UP007
    n1_index: Union[np.array, None] = field(default=None)  # noqa: UP007
    n2_index: Union[np.array, None] = field(default=None)  # noqa: UP007

    def gpc_info(self):
        """Print the log marginal likelihood value and the kernel
        of the Gaussian Process Classifier

        """  # noqa: D205
        print("Training score: ", self.gpc.train_score_)
        print("Marginal likelihood: ", self.gpc.log_marginal_likelihood_value_)
        print("Kernel: ", self.gpc.kernel_)

    @dependencies("matplotlib")
    def plot(
        self,
        X: pd.DataFrame | None = None,
        y: pd.Series | None = None,
        filter: Iterable | None = None,
        background_filter: Iterable | None = None,
        markersize: int = 50,
        fontsize: int = 10,
        chosen_colormap: str = "seismic",
        path_to_save: str | None = None,
        dpi: int = 1000,
        margin: float = 0.1,
        manual_limits: tuple | None = None,
        legend_loc: str = "upper right",
    ) -> None:
        """Plot the results of the differential expression analysis.

        Parameters
        ----------
        X
            Design matrix: log fold changes for group 1 vs group 2 and neighbors 1 vs group 2.
            If None, the design matrix from the classifier is used.
        y
            Labels: FDR values for group 1 vs neighbors 1.
            If None, the labels from the classifier are used.
        filter
            Genes to include in the plot.
        background_filter
            Genes to include in the plot as background.
        markersize
            Size of the markers.
        fontsize
            Font size for the gene names.
        chosen_colormap
            Colormap to use for the decision boundary.
        path_to_save
            Path to save the plot.
        dpi
            Resolution of the saved plot.
        margin
            Margin to add to the plot limits.
        manual_limits
            Manual limits for the plot.
        legend_loc
            Location of the legend.

        Returns
        -------
        None
        """
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.font_manager import FontProperties
        from sklearn.inspection import DecisionBoundaryDisplay

        # Set SVG font type to 'none' to keep text as text in SVG files
        plt.rcParams["svg.fonttype"] = "none"
        plt.rcParams["font.size"] = fontsize
        italic_font = FontProperties(style="italic", size=fontsize - 2)

        if X is None:
            X = self.gpc.X_
            y = self.gpc.y_

        if filter is None:
            filter = X.index

        X_display = X.loc[filter].values
        fdr_g1_n1_display = y.loc[filter]
        lfc_g1_g2_display = X.iloc[:, 0].loc[filter]
        lfc_n1_g2_display = X.iloc[:, 1].loc[filter]

        self.gpc.confident_genes = fdr_g1_n1_display.index[fdr_g1_n1_display]

        if background_filter is not None:
            X_background = X.loc[background_filter].values
            lfc_g1_g2_bg = X.iloc[:, 0].loc[background_filter]
            lfc_n1_g2_bg = X.iloc[:, 1].loc[background_filter]

            # Determine plot limits using both main and background data, with adjustable margin
            x_min = min(lfc_g1_g2_display.min(), lfc_g1_g2_bg.min()) - margin
            x_max = max(lfc_g1_g2_display.max(), lfc_g1_g2_bg.max()) + margin
            y_min = min(lfc_n1_g2_display.min(), lfc_n1_g2_bg.min()) - margin
            y_max = max(lfc_n1_g2_display.max(), lfc_n1_g2_bg.max()) + margin

        else:
            x_min, x_max = lfc_g1_g2_display.min() - margin, lfc_g1_g2_display.max() + margin
            y_min, y_max = lfc_n1_g2_display.min() - margin, lfc_n1_g2_display.max() + margin

        if manual_limits is not None:
            x_min, x_max, y_min, y_max = manual_limits

        # Create a larger figure and axes explicitly
        fig, ax = plt.subplots(figsize=(8, 8))

        # Manually create mesh grid for the decision boundary
        xx, yy = np.meshgrid(
            np.linspace(x_min, x_max, 100),
            np.linspace(y_min, y_max, 100),
        )

        # Create a custom colormap with alpha transparency
        cmap = plt.get_cmap(chosen_colormap)
        cmap_with_alpha = cmap(np.arange(cmap.N))
        cmap_with_alpha[:, -1] = 0.6  # Set the alpha channel
        chosen_colormap = mcolors.ListedColormap(cmap_with_alpha)

        # Generate the decision boundary with custom limits
        disp = DecisionBoundaryDisplay.from_estimator(
            self.gpc,
            np.c_[xx.ravel(), yy.ravel()],
            response_method="predict_proba",
            xlabel="LFC group1-group2",
            ylabel="LFC neighbors1-group2",
            ax=ax,
            cmap=chosen_colormap,
        )

        # Set the limits to ensure the decision boundary display matches the full range of data
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        # Plot the identity line using the retrieved limits with np.linspace
        line_points = np.linspace(min(x_min, y_min), max(x_max, y_max), 100)
        ax.plot(line_points, line_points, "k--", alpha=0.85, zorder=0)

        # Plot background points in light grey if background_filter is provided
        if background_filter is not None:
            disp.ax_.scatter(
                X_background[:, 0],
                X_background[:, 1],
                c="lightgrey",
                edgecolor="none",
                s=markersize,
                alpha=0.3,
            )
            for _i, gene in enumerate(background_filter):
                ax.annotate(
                    gene,
                    (lfc_g1_g2_bg[gene], lfc_n1_g2_bg[gene]),
                    xytext=(0, 5),
                    textcoords="offset points",
                    fontsize=fontsize - 2,
                    fontproperties=italic_font,
                )

        # Scatter plot with fixed colors for True (yellow) and False (green) DE genes
        colors = np.where(fdr_g1_n1_display, "yellow", "green")
        disp.ax_.scatter(X_display[:, 0], X_display[:, 1], c=colors, edgecolor="k", s=markersize)

        # Manually add a colorbar for the decision boundary
        sm = plt.cm.ScalarMappable(cmap=chosen_colormap, norm=mcolors.Normalize(vmin=0, vmax=1))
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label("Decision Boundary Probability", rotation=270, labelpad=15)

        # Add legend for scatter plot
        legend_elements = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label="DE g1_n1 True",
                markerfacecolor="yellow",
                markersize=8,
            ),
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label="DE g1_n1 False",
                markerfacecolor="green",
                markersize=8,
            ),
        ]

        # Add legend for background points if background_filter is provided
        if background_filter is not None:
            legend_elements += [
                plt.Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    label="DE g1_g2 in background",
                    markerfacecolor="lightgrey",
                    markersize=8,
                )
            ]
        disp.ax_.legend(handles=legend_elements, loc=legend_loc)

        # Annotate gene names
        for _i, gene in enumerate(fdr_g1_n1_display.index):
            ax.annotate(
                gene,
                (
                    lfc_g1_g2_display[gene],
                    lfc_n1_g2_display[gene],
                ),
                xytext=(0, 5),
                textcoords="offset points",
                fontsize=fontsize,
            )

        ax.set_aspect("auto")  # Allow the plot to adjust freely

        # Adjust layout and save/show plot
        plt.tight_layout()
        if path_to_save is not None:
            plt.savefig(
                path_to_save,
                bbox_inches="tight",
                dpi=dpi,
            )

        plt.show()
