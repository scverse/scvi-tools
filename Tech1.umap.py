use_cuda = True
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import run_model
from scvi.harmonization.utils_chenling import entropy_batch_mixing
from scvi.metrics.clustering import select_indices_evenly
from scvi.dataset.dataset import SubsetGenes

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from umap import UMAP
plotname = 'Tech1'

from scvi.dataset.muris_tabula import TabulaMuris
dataset1 = TabulaMuris('facs')
dataset2 = TabulaMuris('droplet')
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)

gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)

seurat, batch_indices, labels, keys, stats = run_model('readSeurat', gene_dataset, dataset1, dataset2,
                                                       filename=plotname)

dataset1, dataset2, gene_dataset = SubsetGenes(dataset1, dataset2, gene_dataset, plotname)
vae,  _, _, _, _ = run_model('vae', gene_dataset, dataset1, dataset2,filename=plotname, rep='0')

sorted_key = [
    'hematopoietic precursor cell','Slamf1-positive multipotent progenitor cell','Slamf1-negative multipotent progenitor cell',
    'common lymphoid progenitor',
    'immature T cell','T cell','regulatory T cell',
    'immature NK T cell','immature natural killer cell','pre-natural killer cell','mature natural killer cell',
    'Fraction A pre-pro B cell','early pro-B cell','pro-B cell','late pro-B cell','immature B cell','naive B cell','B cell',
    'granulocyte monocyte progenitor cell','granulocytopoietic cell','granulocyte',
    'promonocyte','monocyte',
    'megakaryocyte-erythroid progenitor cell','proerythroblast', 'erythroblast',
    'basophil',
    'macrophage',
    'nan']

key_order = [list(keys).index(x) for x in sorted_key]


colors = sns.cubehelix_palette(3, start=0,rot=0,light=0.65,dark=0.3,hue=1) + \
sns.cubehelix_palette(1,start=2.7,rot=0,light=0.65, dark=0.3,hue=1) + \
sns.cubehelix_palette(3, start=0.3,rot=0, light=0.65, dark=0.3,hue=1) + \
sns.cubehelix_palette(4, start=2.4,rot=0, light=0.65, dark=0.3,hue=1) + \
sns.cubehelix_palette(7,start=0.6,rot=0, light=0.65, dark=0.3,hue=1) + \
sns.cubehelix_palette(3,start=2.1,rot=0, light=0.65, dark=0.3,hue=1) + \
sns.cubehelix_palette(2,start=0.9,rot=0, light=0.65, dark=0.3,hue=1) + \
sns.cubehelix_palette(3,start=1.8,rot=0, light=0.65, dark=0.3,hue=1) + \
sns.cubehelix_palette(1,start=1.2,rot=0, light=0.65, dark=0.3,hue=1) + \
sns.cubehelix_palette(1,start=1.5,rot=0, light=0.65, dark=0.3,hue=1) + \
sns.cubehelix_palette(1,start=3,rot=0, light=0.65, dark=0.3,hue=1)

latent = vae
sample = select_indices_evenly(3000, batch_indices)
latent_s = latent[sample, :]
label_s = labels[sample]
batch_s = batch_indices[sample]
if latent_s.shape[1] != 2:
    latent_s = UMAP().fit_transform(latent_s)

fig, ax = plt.subplots(figsize=(18, 18))
for i, k in enumerate(key_order):
    ax.scatter(latent_s[label_s == k, 0], latent_s[label_s == k, 1], c=colors[i], label=keys[k],
               edgecolors='none')
    ax.legend(bbox_to_anchor=(1.1, 0.5), borderaxespad=0, fontsize='x-large')

fig.patch.set_visible(False)
ax.axis('off')
fig.tight_layout()
plt.savefig('../' + plotname + '.vae.pretty.labels.pdf')

plt.figure(figsize=(18, 18))
plt.scatter(latent_s[:, 0], latent_s[:, 1], c=batch_s, edgecolors='none')
plt.axis("off")
plt.tight_layout()
plt.savefig('../' + plotname + '.vae.pretty.batchid.pdf')


latent = seurat
sample = select_indices_evenly(2000, labels)
if plotname is not None:
    latent_s = latent[sample, :]
    label_s = labels[sample]
    batch_s = batch_indices[sample]
    if latent_s.shape[1] != 2:
        latent_s = UMAP().fit_transform(latent_s)
    fig, ax = plt.subplots(figsize=(18, 18))
    key_order = np.argsort(keys)
    for i, k in enumerate(key_order):
        ax.scatter(latent_s[label_s == k, 0], latent_s[label_s == k, 1], c=colors[i % 30], label=keys[k],
                   edgecolors='none')
    fig.patch.set_visible(False)
    ax.axis('off')
    fig.tight_layout()
    plt.savefig('../' + plotname + '.cca.pretty.labels.pdf')
    plt.figure(figsize=(18, 18))
    plt.scatter(latent_s[:, 0], latent_s[:, 1], c=batch_s, edgecolors='none')
    plt.axis("off")
    plt.tight_layout()
    plt.savefig('../' + plotname + '.cca.pretty.batchid.pdf')
