import umap
import seaborn as sns
import matplotlib.pyplot as plt

import numpy as np
from scipy.stats import spearmanr


def plot_umap(trainer):
    latent_seq, latent_fish = trainer.get_latent()
    latent2d = umap.UMAP().fit_transform(np.concatenate([latent_seq, latent_fish]))
    latent2d_seq = latent2d[: latent_seq.shape[0]]
    latent2d_fish = latent2d[latent_seq.shape[0] :]

    data_seq, data_fish = [p.gene_dataset for p in trainer.all_dataset]

    colors = sns.color_palette(n_colors=30)
    plt.figure(figsize=(25, 10))
    ax = plt.subplot(1, 3, 1)
    ax.scatter(*latent2d_seq.T, color="r", label="seq", alpha=0.5, s=0.5)
    ax.scatter(*latent2d_fish.T, color="b", label="osm", alpha=0.5, s=0.5)
    ax.legend()

    ax = plt.subplot(1, 3, 2)
    labels = data_seq.labels.ravel()
    for i, label in enumerate(data_seq.cell_types):
        ax.scatter(
            *latent2d_seq[labels == i].T,
            color=colors[i],
            label=label[:12],
            alpha=0.5,
            s=5
        )
    ax.legend()
    ax.set_title("Seq cells")

    ax = plt.subplot(1, 3, 3)
    labels = data_fish.labels.ravel()
    for i, label in enumerate(data_fish.cell_types):
        ax.scatter(
            *latent2d_fish[labels == i].T, color=colors[i], label=label, alpha=0.5, s=5
        )
    ax.legend()
    ax.set_title("Spatial cells")


def imputation_score(trainer_both, data_spatial, gene_ids_test, normalized=True):
    _, fish_imputation = trainer_both.get_imputed_values(normalized=normalized)
    original, imputed = (
        data_spatial.X[:, gene_ids_test],
        fish_imputation[:, gene_ids_test],
    )

    if normalized:
        original /= data_spatial.X.sum(axis=1).reshape(-1, 1)

    spearman_gene = []
    for g in range(imputed.shape[1]):
        if np.all(imputed[:, g] == 0):
            correlation = 0
        else:
            correlation = spearmanr(original[:, g], imputed[:, g])[0]
        spearman_gene.append(correlation)
    return np.median(np.array(spearman_gene))


def plot_gene_spatial(trainer, data_spatial, gene):
    data_seq, _ = [p.gene_dataset for p in trainer.all_dataset]
    data_fish = data_spatial

    fig, (ax_gt, ax) = plt.subplots(1, 2)

    if type(gene) == str:
        gene_id = list(data_seq.gene_names).index(gene)
    else:
        gene_id = gene

    x_coord = data_fish.x_coord.ravel()
    y_coord = data_fish.y_coord.ravel()

    def order_by_strenght(x, y, z):
        ind = np.argsort(z)
        return x[ind], y[ind], z[ind]

    s = 20

    def transform(data):
        return np.log(1 + 100 * data)

    # Plot groundtruth
    x, y, z = order_by_strenght(
        x_coord, y_coord, data_fish.X[:, gene_id] / (data_fish.X.sum(axis=1) + 1)
    )
    ax_gt.scatter(x, y, c=transform(z), s=s, edgecolors="none", marker="s", cmap="Reds")
    ax_gt.set_title("Groundtruth")
    ax_gt.axis("off")

    _, imputed = trainer.get_imputed_values(normalized=True)
    x, y, z = order_by_strenght(x_coord, y_coord, imputed[:, gene_id])
    ax.scatter(x, y, c=transform(z), s=s, edgecolors="none", marker="s", cmap="Reds")
    ax.set_title("Imputed")
    ax.axis("off")
    plt.tight_layout()
    plt.show()
