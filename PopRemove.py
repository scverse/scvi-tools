use_cuda = True
from sklearn.neighbors import KNeighborsClassifier
from scvi.harmonization.utils_chenling import get_matrix_from_dir,SubsetGenes
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
from scvi.dataset.dataset import SubsetGenes

from scvi.metrics.clustering import select_indices_evenly
from umap import UMAP

from scvi.models.vae import VAE
from scvi.models.scanvi import SCANVI
from scvi.inference import UnsupervisedTrainer, SemiSupervisedTrainer
import torch
import os

plotname = 'PopRemove'

dataset1 = PbmcDataset(filter_out_de_genes=False)
dataset1.update_cells(dataset1.batch_indices.ravel()==0)
newCellType = [k for i, k in enumerate(dataset1.cell_types) if k not in ['Other']]
dataset1.filter_cell_types(newCellType)
dataset1.subsample_genes(dataset1.nb_genes)

count, geneid, cellid = get_matrix_from_dir('cite')
count = count.T.tocsr()
seurat = np.genfromtxt('../cite/cite.seurat.labels', dtype='str', delimiter=',')
cellid = np.asarray([x.split('-')[0] for x in cellid])
labels_map = [0, 0, 1, 2, 3, 4, 5, 6]
labels = seurat[1:, 4]
cell_type = ['CD4 T cells', 'NK cells', 'CD14+ Monocytes', 'B cells','CD8 T cells', 'FCGR3A+ Monocytes', 'Other']
dataset2 = assign_label(cellid, geneid, labels_map, count, cell_type, seurat)
newCellType = [k for i, k in enumerate(dataset2.cell_types) if k not in ['Other']]
dataset2.filter_cell_types(newCellType)


dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)


def entropy_from_indices(indices):
    return entropy(np.array(np.unique(indices, return_counts=True)[1].astype(np.int32)))


def entropy_batch_mixing_subsampled(latent, batches, labels, removed_type, n_neighbors=30, n_pools=1, n_samples_per_pool=100):
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


def plotUMAP(latent,plotname,method,rmCellTypes):
    colors = sns.color_palette('tab20')
    sample = select_indices_evenly(2000, labels)
    latent_s = latent[sample, :]
    label_s = labels[sample]
    batch_s = batch_indices[sample]
    if latent_s.shape[1] != 2:
        latent_s = UMAP(spread=2).fit_transform(latent_s)
    fig, ax = plt.subplots(figsize=(10, 10))
    key_order = np.argsort(keys)
    for i, k in enumerate(key_order):
        ax.scatter(latent_s[label_s == k, 0], latent_s[label_s == k, 1], c=colors[i % 20], label=keys[k],
                   edgecolors='none')
        # ax.legend(bbox_to_anchor=(1.1, 0.5), borderaxespad=0, fontsize='x-large')
    ax.axis("off")
    fig.tight_layout()
    plt.savefig('../%s/%s.%s.%s.labels.pdf' % (plotname, plotname, method, rmCellTypes))
    plt.figure(figsize=(10, 10))
    plt.scatter(latent_s[:, 0], latent_s[:, 1], c=batch_s, edgecolors='none')
    plt.axis("off")
    plt.tight_layout()
    plt.savefig('../%s/%s.%s.%s.batches.pdf' % (plotname, plotname, method, rmCellTypes))


def KNNacc(latent, labels, keys):
    neigh = KNeighborsClassifier(n_neighbors=10)
    neigh = neigh.fit(latent, labels)
    labels_pred = neigh.predict(latent)
    acc = [np.mean(labels[labels_pred == i] == i) for i in np.unique(labels)]
    cell_type = keys[np.unique(labels)]
    return acc, cell_type

def BEbyType(keys,latent,labels,batch_indices):
    rm_idx = np.arange(len(keys))[keys == rmCellTypes][0]
    other_idx = np.arange(len(keys))[keys != rmCellTypes]
    cell_type = [keys[rm_idx]] + list(keys[other_idx])
    BE1 = entropy_batch_mixing_subsampled(latent, batch_indices, labels, removed_type=rm_idx)
    BE2 = [entropy_batch_mixing_subsampled(latent, batch_indices, labels, removed_type=i) for i in other_idx]
    res = [BE1] + BE2
    return(res, cell_type)



f = open('../PopRemove/'+plotname+'.acc.res.txt', "w+")
g = open('../PopRemove/'+plotname+'.res.txt', "w+")
# f.write('model_type\tcell_type\tBE_removed\tBE_kept\tasw\tnmi\tari\tca\twca\n')


# scp chenlingantelope@s128.millennium.berkeley.edu:/data/yosef2/users/chenling/harmonization/Seurat_data/PopRemove* .
for rmCellTypes in dataset2.cell_types[:6]:
    pbmc = deepcopy(dataset1)
    newCellType = [k for i, k in enumerate(dataset1.cell_types) if k not in [rmCellTypes,'Others']]
    pbmc.filter_cell_types(newCellType)
    gene_dataset = GeneExpressionDataset.concat_datasets(pbmc, dataset2)
    pbmc = deepcopy(gene_dataset)
    pbmc.update_cells(pbmc.batch_indices.ravel() == 0)
    pbmc.subsample_genes(pbmc.nb_genes)
    pbmc2 = deepcopy(gene_dataset)
    pbmc2.update_cells(gene_dataset.batch_indices.ravel() == 1)
    pbmc2.subsample_genes(dataset2.nb_genes)
    # _,_,_,_,_ = run_model('writedata', gene_dataset, pbmc, pbmc2,filename=plotname+rmCellTypes.replace(' ',''))

    latent, batch_indices, labels, keys, stats = run_model(
        'readSeurat', gene_dataset, pbmc, pbmc2, filename=plotname + rmCellTypes.replace(' ', ''))
    acc, cell_type = KNNacc(latent, labels, keys)
    f.write('Seurat' + '\t' + rmCellTypes + ("\t%.4f" * 8 + "\t%s" * 8 + "\n") % tuple(acc + list(cell_type)))
    be, cell_type2 = BEbyType(keys, latent, labels, batch_indices)
    g.write('Seurat' + '\t' + rmCellTypes + ("\t%.4f" * 8 + "\t%s" * 8 + "\n") % tuple(be + list(cell_type2)))
    # plotUMAP(latent, plotname, 'Seurat', rmCellTypes)

    pbmc, pbmc2, gene_dataset = SubsetGenes(pbmc, pbmc2, gene_dataset, plotname + rmCellTypes.replace(' ', ''))
    vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches, n_labels=gene_dataset.n_labels,
              n_hidden=128, n_latent=10, n_layers=2, dispersion='gene')
    trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=1.0)
    if os.path.isfile('../PopRemove/vae.%s.pkl' % rmCellTypes):
        trainer.model.load_state_dict(torch.load('../PopRemove/vae.%s.pkl' % rmCellTypes))
        trainer.model.eval()
    else:
        trainer.train(n_epochs=150)
        torch.save(trainer.model.state_dict(), '../PopRemove/vae.%s.pkl' % rmCellTypes)
    full = trainer.create_posterior(trainer.model, gene_dataset, indices=np.arange(len(gene_dataset)))
    latent, batch_indices, labels = full.sequential().get_latent()
    batch_indices = batch_indices.ravel()
    acc, cell_type = KNNacc(latent, labels, keys)
    f.write('vae' + '\t' + rmCellTypes + ("\t%.4f" * 8 + "\t%s" * 8 + "\n") % tuple(acc + list(cell_type)))
    be, cell_type2 = BEbyType(keys, latent, labels, batch_indices)
    g.write('vae' + '\t' + rmCellTypes + ("\t%.4f" * 8 + "\t%s" * 8 + "\n") % tuple(be + list(cell_type2)))
    # plotUMAP(latent, plotname, 'vae', rmCellTypes)

    scanvi = SCANVI(gene_dataset.nb_genes, gene_dataset.n_batches, (gene_dataset.n_labels), n_layers=2)
    scanvi.load_state_dict(trainer.model.state_dict(), strict=False)
    trainer_scanvi = SemiSupervisedTrainer(scanvi, gene_dataset, classification_ratio=50,
                                           n_epochs_classifier=1, lr_classification=5 * 1e-3)

    trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=gene_dataset.batch_indices.ravel() == 0)
    trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(indices=gene_dataset.batch_indices.ravel() == 1)
    if os.path.isfile('../PopRemove/scanvi.%s.pkl' % rmCellTypes):
        trainer_scanvi.model.load_state_dict(torch.load('../PopRemove/scanvi.%s.pkl' % rmCellTypes))
        trainer_scanvi.model.eval()
    else:
        trainer_scanvi.train(n_epochs=3)
        torch.save(trainer_scanvi.model.state_dict(), '../PopRemove/scanvi.%s.pkl' % rmCellTypes)
    scanvi_full = trainer_scanvi.create_posterior(trainer_scanvi.model, gene_dataset, indices=np.arange(len(gene_dataset)))
    latent, _, _ = scanvi_full.sequential().get_latent()
    batch_indices = batch_indices.ravel()
    acc, cell_type = KNNacc(latent, labels, keys)
    f.write('scanvi' + '\t' + rmCellTypes + ("\t%.4f" * 8 + "\t%s" * 8 + "\n") % tuple(acc + list(cell_type)))
    be, cell_type2 = BEbyType(keys, latent, labels, batch_indices)
    g.write('scanvi' + '\t' + rmCellTypes + ("\t%.4f" * 8 + "\t%s" * 8 + "\n") % tuple(be + list(cell_type2)))
    # plotUMAP(latent, plotname, 'scanvi', rmCellTypes)

f.close()
g.close()
