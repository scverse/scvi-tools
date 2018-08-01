import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def show_cell_types(latent_seq, labels_seq, latent_fish, labels_fish, title='latent.svg', x_lim=[-85, 85],
                    y_lim=[-85, 85]):
    plt.figure(figsize=(15, 15))

    gs = gridspec.GridSpec(4, 5)
    ax_common = plt.subplot(gs[:, 2:])
    ax_seq = plt.subplot(gs[:2, :2])
    ax_fish = plt.subplot(gs[2:, :2])

    # Plotting fish data
    ax_common.scatter(latent_fish[:, 0][labels_fish == 0], latent_fish[:, 1][labels_fish == 0],
                      c='g', edgecolors='none', label='astrocytes')
    ax_common.scatter(latent_fish[:, 0][labels_fish == 2], latent_fish[:, 1][labels_fish == 2],
                      c='r', edgecolors='none', label='interneurons')
    ax_common.scatter(latent_fish[:, 0][labels_fish == 4], latent_fish[:, 1][labels_fish == 4],
                      c='b', edgecolors='none', label='oligodendrocytes')
    ax_common.scatter(latent_fish[:, 0][labels_fish == 1], latent_fish[:, 1][labels_fish == 1],
                      c='k', edgecolors='none', label='endothelial')
    ax_common.scatter(latent_fish[:, 0][labels_fish == 3], latent_fish[:, 1][labels_fish == 3],
                      c='m', edgecolors='none', label='microglia')
    ax_common.scatter(latent_fish[:, 0][labels_fish == 5], latent_fish[:, 1][labels_fish == 5],
                      c='y', edgecolors='none', label='pyramidal')
    ax_fish.scatter(latent_fish[:, 0][labels_fish == 0], latent_fish[:, 1][labels_fish == 0],
                    c='g', edgecolors='none', label='astrocytes')
    ax_fish.scatter(latent_fish[:, 0][labels_fish == 2], latent_fish[:, 1][labels_fish == 2],
                    c='r', edgecolors='none', label='interneurons')
    ax_fish.scatter(latent_fish[:, 0][labels_fish == 4], latent_fish[:, 1][labels_fish == 4],
                    c='b', edgecolors='none', label='oligodendrocytes')
    ax_fish.scatter(latent_fish[:, 0][labels_fish == 1], latent_fish[:, 1][labels_fish == 1],
                    c='k', edgecolors='none', label='endothelial')
    ax_fish.scatter(latent_fish[:, 0][labels_fish == 3], latent_fish[:, 1][labels_fish == 3],
                    c='m', edgecolors='none', label='microglia')
    ax_fish.scatter(latent_fish[:, 0][labels_fish == 5], latent_fish[:, 1][labels_fish == 5],
                    c='y', edgecolors='none', label='pyramidal')

    # Plotting seq data
    ax_common.scatter(latent_seq[:, 0][labels_seq == 0], latent_seq[:, 1][labels_seq == 0], c='g', edgecolors='none')
    ax_common.scatter(latent_seq[:, 0][labels_seq == 2], latent_seq[:, 1][labels_seq == 2], c='r', edgecolors='none')
    ax_common.scatter(latent_seq[:, 0][labels_seq == 4], latent_seq[:, 1][labels_seq == 4], c='b', edgecolors='none')
    ax_common.scatter(latent_seq[:, 0][labels_seq == 1], latent_seq[:, 1][labels_seq == 1], c='k', edgecolors='none')
    ax_common.scatter(latent_seq[:, 0][labels_seq == 3], latent_seq[:, 1][labels_seq == 3], c='m', edgecolors='none')
    ax_common.scatter(latent_seq[:, 0][labels_seq == 5], latent_seq[:, 1][labels_seq == 5], c='y', edgecolors='none')
    ax_common.scatter(latent_seq[:, 0][labels_seq == 6], latent_seq[:, 1][labels_seq == 6], c='y', edgecolors='none')
    ax_seq.scatter(latent_seq[:, 0][labels_seq == 0], latent_seq[:, 1][labels_seq == 0], c='g', edgecolors='none',
                   label='astrocytes')
    ax_seq.scatter(latent_seq[:, 0][labels_seq == 2], latent_seq[:, 1][labels_seq == 2], c='r', edgecolors='none',
                   label='interneurons')
    ax_seq.scatter(latent_seq[:, 0][labels_seq == 4], latent_seq[:, 1][labels_seq == 4], c='b', edgecolors='none',
                   label='oligodendrocytes')
    ax_seq.scatter(latent_seq[:, 0][labels_seq == 1], latent_seq[:, 1][labels_seq == 1], c='k', edgecolors='none',
                   label='endothelial')
    ax_seq.scatter(latent_seq[:, 0][labels_seq == 3], latent_seq[:, 1][labels_seq == 3], c='m', edgecolors='none',
                   label='microglia')
    ax_seq.scatter(latent_seq[:, 0][labels_seq == 5], latent_seq[:, 1][labels_seq == 5], c='y', edgecolors='none',
                   label='pyramidal')
    ax_seq.scatter(latent_seq[:, 0][labels_seq == 6], latent_seq[:, 1][labels_seq == 6], c='y', edgecolors='none',
                   label='pyramidal')
    ax_common.legend(loc="best")
    ax_seq.legend(loc="best")
    ax_fish.legend(loc="best")
    ax_seq.set_title('Latent embedding of scRNA cells')
    ax_fish.set_title('Latent embedding of smFISH_utils cells')
    ax_common.set_title('Common latent embedding')
    if x_lim is not None:
        ax_seq.set_xlim([-85, 85])
        ax_fish.set_xlim([-85, 85])
        ax_common.set_xlim([-85, 85])
        ax_seq.set_ylim([-85, 85])
        ax_fish.set_ylim([-85, 85])
        ax_common.set_ylim([-85, 85])
    plt.show()
    plt.savefig(title)


def show_mixing(latent_seq, latent_fish, title='mixing.svg', x_lim=[-85, 85], y_lim=[-85, 85]):
    plt.figure(figsize=(15, 15))
    gs = gridspec.GridSpec(4, 5)
    ax_common = plt.subplot(gs[:, 2:])
    ax_seq = plt.subplot(gs[:2, :2])
    ax_fish = plt.subplot(gs[2:, :2])
    ax_seq.scatter(latent_seq[:, 0], latent_seq[:, 1], c='m', label="scRNA cells")
    ax_common.scatter(latent_seq[:, 0], latent_seq[:, 1], c='m', label="scRNA cells")
    ax_fish.scatter(latent_fish[:, 0], latent_fish[:, 1], c='y', label="smFISH_utils cells")
    ax_common.scatter(latent_fish[:, 0], latent_fish[:, 1], c='y', label="smFISH_utils cells")
    ax_common.legend(loc="best")
    ax_seq.legend(loc="best")
    ax_fish.legend(loc="best")
    ax_seq.set_title('Latent embedding of scRNA cells')
    ax_fish.set_title('Latent embedding of smFISH_utils cells')
    ax_common.set_title('Common latent embedding')
    if x_lim is not None:
        ax_seq.set_xlim([-85, 85])
        ax_fish.set_xlim([-85, 85])
        ax_common.set_xlim([-85, 85])
        ax_seq.set_ylim([-85, 85])
        ax_fish.set_ylim([-85, 85])
        ax_common.set_ylim([-85, 85])
    plt.show()
    plt.savefig(title)


def compare_cell_types(latent, labels, inferred_labels, title='labelling.svg', x_lim=[-85, 85], y_lim=[-85, 85]):
    plt.figure(figsize=(10, 15))
    gs = gridspec.GridSpec(2, 2)
    ax_real = plt.subplot(gs[:1, :2])
    ax_inferred = plt.subplot(gs[1:, :2])

    # Plotting fish data
    ax_real.scatter(latent[:, 0][labels == 0], latent[:, 1][labels == 0],
                    c='g', edgecolors='none', label='astrocytes')
    ax_real.scatter(latent[:, 0][labels == 2], latent[:, 1][labels == 2],
                    c='r', edgecolors='none', label='interneurons')
    ax_real.scatter(latent[:, 0][labels == 4], latent[:, 1][labels == 4],
                    c='b', edgecolors='none', label='oligodendrocytes')
    ax_real.scatter(latent[:, 0][labels == 1], latent[:, 1][labels == 1],
                    c='k', edgecolors='none', label='endothelial')
    ax_real.scatter(latent[:, 0][labels == 3], latent[:, 1][labels == 3],
                    c='m', edgecolors='none', label='microglia')
    ax_real.scatter(latent[:, 0][labels == 5], latent[:, 1][labels == 5],
                    c='y', edgecolors='none', label='pyramidal')
    ax_inferred.scatter(latent[:, 0][inferred_labels == 0], latent[:, 1][inferred_labels == 0],
                        c='g', edgecolors='none', label='astrocytes')
    ax_inferred.scatter(latent[:, 0][inferred_labels == 2], latent[:, 1][inferred_labels == 2],
                        c='r', edgecolors='none', label='interneurons')
    ax_inferred.scatter(latent[:, 0][inferred_labels == 4], latent[:, 1][inferred_labels == 4],
                        c='b', edgecolors='none', label='oligodendrocytes')
    ax_inferred.scatter(latent[:, 0][inferred_labels == 1], latent[:, 1][inferred_labels == 1],
                        c='k', edgecolors='none', label='endothelial')
    ax_inferred.scatter(latent[:, 0][inferred_labels == 3], latent[:, 1][inferred_labels == 3],
                        c='m', edgecolors='none', label='microglia')
    ax_inferred.scatter(latent[:, 0][inferred_labels == 5], latent[:, 1][inferred_labels == 5],
                        c='y', edgecolors='none', label='pyramidal')

    ax_real.legend(loc="best")
    ax_inferred.legend(loc="best")
    ax_real.set_title('Real cell types')
    ax_inferred.set_title('Inferred cell types')
    if x_lim is not None:
        ax_real.set_xlim([-85, 85])
        ax_inferred.set_xlim([-85, 85])
        ax_real.set_ylim([-85, 85])
        ax_inferred.set_ylim([-85, 85])
    plt.show()
    plt.savefig(title)


def show_gene_exp(latent, gene_exp, labels=None, title='latent_expression.svg', gene_name=None, x_lim=[-85, 85],
                  y_lim=[-85, 85]):
    if labels is None:
        plt.figure(figsize=(15, 15))
        plt.scatter(latent[:, 0], latent[:, 1], c=gene_exp.ravel(), label=gene_name)
        plt.legend()
        plt.savefig(title)
    else:
        plt.figure(figsize=(10, 15))
        gs = gridspec.GridSpec(2, 2)
        ax_real = plt.subplot(gs[:1, :2])
        ax_inferred = plt.subplot(gs[1:, :2])

        # Plotting fish data
        ax_real.scatter(latent[:, 0], latent[:, 1], c=gene_exp.ravel(), label=gene_name)
        ax_inferred.scatter(latent[:, 0][labels == 0], latent[:, 1][labels == 0],
                            c='g', edgecolors='none', label='astrocytes')
        ax_inferred.scatter(latent[:, 0][labels == 2], latent[:, 1][labels == 2],
                            c='r', edgecolors='none', label='interneurons')
        ax_inferred.scatter(latent[:, 0][labels == 4], latent[:, 1][labels == 4],
                            c='b', edgecolors='none', label='oligodendrocytes')
        ax_inferred.scatter(latent[:, 0][labels == 1], latent[:, 1][labels == 1],
                            c='k', edgecolors='none', label='endothelial')
        ax_inferred.scatter(latent[:, 0][labels == 3], latent[:, 1][labels == 3],
                            c='m', edgecolors='none', label='microglia')
        ax_inferred.scatter(latent[:, 0][labels == 5], latent[:, 1][labels == 5],
                            c='y', edgecolors='none', label='pyramidal')

        ax_real.legend(loc="best")
        ax_inferred.legend(loc="best")
        ax_real.set_title("Gene" + gene_name + " expression")
        ax_inferred.set_title('Cell types')
        if x_lim is not None:
            ax_real.set_xlim([-85, 85])
            ax_inferred.set_xlim([-85, 85])
            ax_real.set_ylim([-85, 85])
            ax_inferred.set_ylim([-85, 85])
        plt.show()
        plt.savefig(title)


def show_spatial_expression(x_coords, y_coords, gene_exp, labels=None, title='spatial_expression.svg', gene_name=None):
    if labels is None:
        plt.figure(figsize=(15, 15))
        x_coords = x_coords.ravel()
        y_coords = y_coords.ravel()
        plt.scatter(x_coords, y_coords, c=gene_exp.ravel(), label=gene_name)
        plt.legend()
        plt.savefig(title)
    else:
        plt.figure(figsize=(10, 15))
        gs = gridspec.GridSpec(2, 2)
        ax_real = plt.subplot(gs[:1, :2])
        ax_inferred = plt.subplot(gs[1:, :2])

        # Plotting fish data
        x_coords = x_coords.ravel()
        y_coords = y_coords.ravel()
        ax_real.scatter(x_coords, y_coords, c=gene_exp.ravel(), label=gene_name)
        ax_inferred.scatter(x_coords[labels == 0], y_coords[labels == 0],
                            c='g', edgecolors='none', label='astrocytes')
        ax_inferred.scatter(x_coords[labels == 2], y_coords[labels == 2],
                            c='r', edgecolors='none', label='interneurons')
        ax_inferred.scatter(x_coords[labels == 4], y_coords[labels == 4],
                            c='b', edgecolors='none', label='oligodendrocytes')
        ax_inferred.scatter(x_coords[labels == 1], y_coords[labels == 1],
                            c='k', edgecolors='none', label='endothelial')
        ax_inferred.scatter(x_coords[labels == 3], y_coords[labels == 3],
                            c='m', edgecolors='none', label='microglia')
        ax_inferred.scatter(x_coords[labels == 5], y_coords[labels == 5],
                            c='y', edgecolors='none', label='pyramidal')

        ax_real.legend(loc="best")
        ax_inferred.legend(loc="best")
        ax_real.set_title("Gene" + gene_name + " expression")
        ax_inferred.set_title('Cell types')
        plt.show()
        plt.savefig(title)
