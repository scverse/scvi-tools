from scvi.dataset.scanorama import DatasetSCANORAMA
from scvi.dataset.dataset import GeneExpressionDataset
import numpy as np
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns
from scvi.metrics.clustering import select_indices_evenly
from sklearn.manifold import TSNE
dirs = (
        open('/data/scanorama/conf/panorama.txt')
        .read().rstrip().split()
)
# ['data/293t_jurkat/293t', 'data/293t_jurkat/jurkat', 'data/293t_jurkat/jurkat_293t_50_50', 'data/293t_jurkat/jurkat_293t_99_1',
#  'data/brain/neuron_9k',
#  'data/macrophage/infected', 'data/macrophage/mixed_infected', 'data/macrophage/uninfected', 'data/macrophage/uninfected_donor2',
#  'data/hsc/hsc_mars', 'data/hsc/hsc_ss2',
#  'data/pancreas/pancreas_inDrop', 'data/pancreas/pancreas_multi_celseq2_expression_matrix', 'data/pancreas/pancreas_multi_celseq_expression_matrix', 'data/pancreas/pancreas_multi_fluidigmc1_expression_matrix', 'data/pancreas/pancreas_multi_smartseq2_expression_matrix',
#  'data/pbmc/10x/68k_pbmc', 'data/pbmc/10x/b_cells', 'data/pbmc/10x/cd14_monocytes', 'data/pbmc/10x/cd4_t_helper', 'data/pbmc/10x/cd56_nk', 'data/pbmc/10x/cytotoxic_t', 'data/pbmc/10x/memory_t', 'data/pbmc/10x/regulatory_t', 'data/pbmc/pbmc_kang', 'data/pbmc/pbmc_10X']

datasets = [DatasetSCANORAMA(d) for d in dirs]

all_dataset = GeneExpressionDataset.concat_datasets(*datasets)
# Keeping 5216 genes

labels = (open('/data/scanorama/data/cell_labels/all.txt').read().rstrip().split())
all_dataset.cell_types,all_dataset.labels = np.unique(labels,return_inverse=True)
all_dataset.labels = all_dataset.labels.reshape(len(all_dataset.labels),1)
all_dataset.batch_indices = all_dataset.batch_indices.astype('int')

from scvi.harmonization.utils_chenling import trainVAE,VAEstats
# full = trainVAE(all_dataset, 'scanorama', 1, nlayers=3,n_hidden=256)
full = trainVAE(all_dataset, 'scanorama', 0) #  nlayers=2,n_hidden=128
# full = trainVAE(all_dataset, 'scanorama', 2, nlayers=3,n_hidden=128)

# for 250 iterations, takes 45:14 to train VAE
latent, batch_indices, labels, stats = VAEstats(full)
from scvi.metrics.clustering import clustering_scores
res_knn = clustering_scores(np.asarray(latent), labels, 'knn')
res_kmeans = clustering_scores(np.asarray(latent), labels, 'KMeans')

celltype_names = ['293t', 'CD19+ B cells', 'Hsc CMP', 'Hsc GMP', 'PBMC CD14+ monocyte', 'Hsc MEP', 'PBMC CD56+ NK', 'PBMC T cells', 'PBMC CD4 T cells', 'PBMC CD8 T cells','Hsc Unknown', 'Pancreas Acinar', 'Pancreas Alpha', 'Pancreas Beta', 'Pancreas Beta(ER stress)', 'Pancreas Delta', 'Pancreas Ductal','Pancreas Endothelial', 'Pancreas Gamma', 'Jurkat', 'Macrophage', 'Neuron', 'Pancreas Unknown','Pancreas  Stellate']
plt.figure(figsize=(5, 4))
order = np.argsort(res_knn['clusteracc'])
plt.bar(np.arange(len(all_dataset.cell_types)), np.asarray(res_knn['clusteracc'])[order])
plt.xticks(np.arange(len(all_dataset.cell_types)), np.asarray(celltype_names)[order],rotation='vertical')
plt.tight_layout()
plt.savefig("/data/scanorama/acc_barplot.pdf")
# # 1:05:58 for the more complex model
#     # , clustering_scores,clustering_accuracy
# plotname='scanorama'
# sample = select_indices_evenly(2000, batch_indices)
# colors = sns.color_palette('bright') + \
#          sns.color_palette('muted') + \
#          sns.color_palette('pastel') + \
#          sns.color_palette('dark') + \
#          sns.color_palette('colorblind')
#
#
# from umap import UMAP
#
# latent_s = latent[sample, :]
# label_s = labels[sample]
# batch_s = batch_indices[sample]
# if latent_s.shape[1] != 2:
#     latent_s = UMAP().fit_transform(latent_s)
#
# keys = all_dataset.cell_types
# fig, ax = plt.subplots(figsize=(20, 18))
# key_order = np.argsort(keys)
# for i, k in enumerate(key_order):
#     ax.scatter(latent_s[label_s == k, 0], latent_s[label_s == k, 1], c=colors[i % 30], label=keys[k],
#                edgecolors='none')
#     # ax.legend(bbox_to_anchor=(1.1, 0.5), borderaxespad=0, fontsize='x-large')
# fig.patch.set_visible(False)
# ax.axis('off')
# fig.tight_layout()
# plt.savefig('../' + plotname + '.2l_128.umap.labels.pdf')
#
#
# keys = dirs
# fig, ax = plt.subplots(figsize=(20, 14))
# key_order = np.argsort(keys)
#
# colors = sns.light_palette("navy",4, reverse=True) + \
# sns.light_palette("gray",1, reverse=True) + \
# sns.light_palette("orange",4, reverse=True) + \
# sns.light_palette("purple",2,reverse=True) + \
# sns.light_palette("green",5,reverse=True) + \
# sns.light_palette("red",10,reverse=True)
#
# for k in key_order:
#     ax.scatter(latent_s[batch_s == k, 0], latent_s[batch_s == k, 1], c=colors[k], label=keys[k],
#                edgecolors='none')
#     # ax.legend(bbox_to_anchor=(1.1, 0.5), borderaxespad=0, fontsize='x-large')
# fig.patch.set_visible(False)
# ax.axis('off')
# fig.tight_layout()
# plt.savefig('../' + plotname + '.2l_128.umap.batches.pdf')
