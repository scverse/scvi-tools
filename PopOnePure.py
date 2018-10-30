from scvi.dataset import Dataset10X

use_cuda = True
from scvi.harmonization.utils_chenling import get_matrix_from_dir
from scvi.dataset.pbmc import PbmcDataset
from scvi.harmonization.utils_chenling import assign_label
import numpy as np
from scipy.stats import entropy
from sklearn.neighbors import NearestNeighbors

from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import run_model
from copy import deepcopy

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns

from scvi.metrics.clustering import select_indices_evenly,clustering_scores
from sklearn.manifold import TSNE



plotname = 'PopOne'

cell_types = np.array(["cd4_t_helper", "regulatory_t", "naive_t", "memory_t", "cytotoxic_t", "naive_cytotoxic",
                       "b_cells", "cd34", "cd56_nk", "cd14_monocytes"])
cell_type_name = np.array(["CD4 T cells", "CD4 T cells Regulatory", "CD4 T cells Naive", "CD4 Memory T cells", "CD8 T cells", "CD8 T cells Naive",
                       "B cells", "CD34 cells", "NK cells", "CD14+ Monocytes"])

datasets = []
for i,cell_type in enumerate(cell_types):
    dataset = Dataset10X(cell_type, save_path='data/')
    dataset.cell_types = np.array([cell_type_name[i]])
    dataset.labels = dataset.labels.astype('int')
    dataset.subsample_genes(dataset.nb_genes)
    dataset.gene_names = dataset.gene_symbols
    datasets += [dataset]

dataset1 = GeneExpressionDataset.concat_datasets(*datasets, shared_batches=True)


dataset2 = PbmcDataset(filter_out_de_genes=False)
dataset2.update_cells(dataset2.batch_indices.ravel()==0)
dataset2.subsample_genes(dataset2.nb_genes)


prop2 = dict( zip(dataset2.cell_types[np.unique(dataset2.labels,return_counts=True)[0]],
                  np.unique(dataset2.labels,return_counts=True)[1]/len(dataset2)))

prop =dict(zip(dataset1.cell_types, np.unique(dataset1.labels,return_counts=True)[1]/len(dataset1)))
for i,name in enumerate(list(set(dataset2.cell_types).intersection(set(dataset2.cell_types)))):
    print("%s & %.3f & %.3f\\\\ \n" % (name,prop[name],prop2[name]))


dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)


def entropy_from_indices(indices):
    return entropy(np.array(np.unique(indices, return_counts=True)[1].astype(np.int32)))


def entropy_batch_mixing_subsampled(latent, batches, labels, removed_type, sampled_batch=0, n_neighbors=50, n_pools=50, n_samples_per_pool=100):
    X = latent[labels == removed_type,:]
    nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(latent)
    indices = nbrs.kneighbors(X, return_distance=False)[:, 1:]
    batch_indices = np.vectorize(lambda i: batches[i])(indices)
    entropies = np.apply_along_axis(entropy_from_indices, axis=1, arr=batch_indices)
    if n_pools == 1:
        res = np.mean(entropies)
    else:
        res = np.mean([
            np.mean(entropies[np.random.choice(len(entropies), size=n_samples_per_pool)])
            for _ in range(n_pools)
        ])
    return res

from scvi.metrics.clustering import clustering_scores

f = open('../'+plotname+'/res.txt', "w+")

from scvi.dataset.dataset import SubsetGenes

# scp chenlingantelope@s128.millennium.berkeley.edu:/data/yosef2/users/chenling/harmonization/Seurat_data/PopOneCD8*

shared = set(dataset1.cell_types).intersection(set(dataset2.cell_types))
shared = list(shared)
for KeepCellTypes in shared:
    pbmc = deepcopy(dataset1)
    pbmc.filter_cell_types([KeepCellTypes])
    gene_dataset = GeneExpressionDataset.concat_datasets(pbmc, dataset2)
    pbmc = deepcopy(gene_dataset)
    pbmc.update_cells(pbmc.batch_indices.ravel() == 0)
    pbmc.subsample_genes(pbmc.nb_genes)
    pbmc2 = deepcopy(gene_dataset)
    pbmc2.update_cells(gene_dataset.batch_indices.ravel() == 1)
    pbmc2.subsample_genes(dataset2.nb_genes)
    latent, batch_indices, labels, keys, stats = run_model(
        'writedata', gene_dataset, pbmc, pbmc2,filename=plotname+KeepCellTypes.replace(' ',''))
    latent, batch_indices, labels, keys, stats = run_model(
        'readSeurat', gene_dataset, pbmc, pbmc2,filename=plotname+KeepCellTypes.replace(' ',''))

    otheridx = np.arange(len(keys))[keys == 'Other'][0]
    latent = latent[labels!=otheridx,:]
    batch_indices = batch_indices[labels!=otheridx]
    labels = labels[labels!=otheridx]
    map = dict(zip(np.unique(labels),np.argsort(np.unique(labels))))
    labels = np.asarray([map[x] for x in labels])
    keys = keys[keys!='Other']

    rm_idx = np.arange(len(keys))[keys == KeepCellTypes][0]
    other_idx = np.arange(len(keys))[keys != KeepCellTypes]
    cell_type = [keys[rm_idx]] + list(keys[other_idx])
    BE1 = entropy_batch_mixing_subsampled(latent, batch_indices, labels, removed_type=rm_idx)
    BE2 = [entropy_batch_mixing_subsampled(latent, batch_indices, labels, removed_type= i ) for i in other_idx]
    res_knn = clustering_scores(np.asarray(latent), labels, 'knn')
    res = [BE1] + BE2 + [res_knn[x] for x in res_knn]
    f.write('Seurat' + '\t' + KeepCellTypes + ("\t%.4f" * len(res) + "\t%s"*len(cell_type) + "\n") % tuple(res+cell_type))

    colors = sns.color_palette('tab20')
    sample = select_indices_evenly(2000, labels)
    latent_s = latent[sample, :]
    label_s = labels[sample]
    batch_s = batch_indices[sample]
    if latent_s.shape[1] != 2:
        latent_s = TSNE().fit_transform(latent_s)

    fig, ax = plt.subplots(figsize=(13, 10))
    key_order = np.argsort(keys)
    for i, k in enumerate(key_order):
        ax.scatter(latent_s[label_s == k, 0], latent_s[label_s == k, 1], c=colors[i % 20], label=keys[k],
                   edgecolors='none')
        ax.legend(bbox_to_anchor=(1.1, 0.5), borderaxespad=0, fontsize='x-large')

    fig.tight_layout()
    plt.savefig('../' + plotname+ '/Seurat' + '.' + KeepCellTypes.replace(' ', '') + '.labels.png')
    plt.figure(figsize=(10, 10))
    plt.scatter(latent_s[:, 0], latent_s[:, 1], c=batch_s, edgecolors='none')
    plt.axis("off")
    plt.tight_layout()
    plt.savefig('../' + plotname + '/Seurat' + '.' + KeepCellTypes.replace(' ', '') + '.batchid.png')

    pbmc, pbmc2, gene_dataset = SubsetGenes(pbmc, pbmc2, gene_dataset, plotname+KeepCellTypes.replace(' ',''))
    latent, batch_indices, labels, keys, stats = run_model(
        'vae', gene_dataset, pbmc, pbmc2,filename=plotname, rep = KeepCellTypes.replace(' ',''))

    otheridx = np.arange(len(keys))[keys == 'Other'][0]
    latent = latent[labels!=otheridx,:]
    batch_indices = batch_indices[labels!=otheridx]
    labels = labels[labels!=otheridx]
    map = dict(zip(np.unique(labels),np.argsort(np.unique(labels))))
    labels = np.asarray([map[x] for x in labels])
    keys = keys[keys!='Other']

    rm_idx = np.arange(len(keys))[keys == KeepCellTypes][0]
    other_idx = np.arange(len(keys))[keys != KeepCellTypes]
    cell_type = [keys[rm_idx]] + list(keys[other_idx])
    BE1 = entropy_batch_mixing_subsampled(latent, batch_indices, labels, removed_type=rm_idx)
    BE2 = [ entropy_batch_mixing_subsampled(latent, batch_indices, labels, removed_type= i ) for i in other_idx]
    res_knn = clustering_scores(np.asarray(latent), labels, 'knn')
    res = [BE1] + BE2 + [res_knn[x] for x in res_knn]
    f.write('vae' + '\t' + KeepCellTypes + ("\t%.4f" * len(res) + "\t%s"* len(cell_type) + "\n") % tuple(res+cell_type))


    colors = sns.color_palette('tab20')
    sample = select_indices_evenly(2000, labels)
    latent_s = latent[sample, :]
    label_s = labels[sample]
    batch_s = batch_indices[sample]
    if latent_s.shape[1] != 2:
        latent_s = TSNE().fit_transform(latent_s)

    fig, ax = plt.subplots(figsize=(13, 10))
    key_order = np.argsort(keys)
    for i, k in enumerate(key_order):
        ax.scatter(latent_s[label_s == k, 0], latent_s[label_s == k, 1], c=colors[i % 20], label=keys[k],
                   edgecolors='none')
        ax.legend(bbox_to_anchor=(1.1, 0.5), borderaxespad=0, fontsize='x-large')

    fig.tight_layout()
    plt.savefig('../'+ plotname + '/vae.' + KeepCellTypes.replace(' ', '') + '.labels.png')
    plt.figure(figsize=(10, 10))
    plt.scatter(latent_s[:, 0], latent_s[:, 1], c=batch_s, edgecolors='none')
    plt.axis("off")
    plt.tight_layout()
    plt.savefig('../'+ plotname + '/vae.' + KeepCellTypes.replace(' ', '') + '.batchid.png')

f.close()
