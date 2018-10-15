from scvi.harmonization.utils_chenling import get_matrix_from_dir, assign_label
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import eval_latent, run_model
from copy import deepcopy
from scvi.inference import UnsupervisedTrainer, AdversarialTrainerVAE
from scvi.inference.posterior import *
from scvi.metrics.clustering import select_indices_evenly, clustering_scores, entropy_batch_mixing
from scvi.models.vae import VAE
from sklearn.neighbors import NearestNeighbors
from scvi.dataset import GeneExpressionDataset, Dataset10X, CiteSeqDataset

use_cuda = True

# EASY1
# count, geneid, cellid = get_matrix_from_dir('pbmc8k')
# geneid = geneid[:, 1]
# count = count.T.tocsr()
# seurat = np.genfromtxt('../pbmc8k/pbmc8k.seurat.labels', dtype='str', delimiter=',')
# cellid = np.asarray([x.split('-')[0] for x in cellid])
# labels_map = [0, 2, 4, 4, 0, 3, 3, 1, 5, 6]
# cell_type = ["CD4+ T Helper2", "CD56+ NK", "CD14+ Monocyte", "CD19+ B", "CD8+ Cytotoxic T", "FCGR3A Monocyte",
#              "dendritic"]
# dataset1 = assign_label(cellid, geneid, labels_map, count, cell_type, seurat)
# rmCellTypes = {'na', 'dendritic'}
# newCellType = [k for i, k in enumerate(dataset1.cell_types) if k not in rmCellTypes]
# dataset1.filter_cell_types(newCellType)
# count, geneid, cellid = get_matrix_from_dir('cite')
# count = count.T.tocsr()
# seurat = np.genfromtxt('../cite/cite.seurat.labels', dtype='str', delimiter=',')
# cellid = np.asarray([x.split('-')[0] for x in cellid])
# labels_map = [0, 0, 1, 2, 3, 4, 5, 6]
# labels = seurat[1:, 4]
# cell_type = ["CD4+ T Helper2", "CD56+ NK", "CD14+ Monocyte", "CD19+ B", "CD8+ Cytotoxic T", "FCGR3A Monocyte", "na"]
# dataset2 = assign_label(cellid, geneid, labels_map, count, cell_type, seurat)
# rmCellTypes = {'na', 'dendritic'}
# newCellType = [k for i, k in enumerate(dataset2.cell_types) if k not in rmCellTypes]
# dataset2.filter_cell_types(newCellType)
# dataset1.subsample_genes(dataset1.nb_genes)
# dataset2.subsample_genes(dataset2.nb_genes)
# gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)
# genes = np.genfromtxt('../Seurat_data/Easy1.CCA.genes.txt')
# genes = genes.astype('int')
# gene_dataset.X = gene_dataset.X[:, genes]
# gene_dataset.update_genes(genes)

from scvi.dataset.muris_tabula import TabulaMuris
dataset1 = TabulaMuris('facs')
dataset2 = TabulaMuris('droplet')
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)
genes = np.genfromtxt('../Seurat_data/Tech1.CCA.genes.txt')
genes = genes.astype('int')
gene_dataset.X = gene_dataset.X[:, genes]
gene_dataset.update_genes(genes)
cells = np.genfromtxt('../Seurat_data/Tech1.CCA.cells.txt')
print(cells.shape)
print(gene_dataset.X.shape)

vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches, n_labels=gene_dataset.n_labels,
          n_hidden=128, n_latent=10, n_layers=1, n_layers_decoder=1, dispersion='gene')
# trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=1.0)
trainer = AdversarialTrainerVAE(vae, gene_dataset, train_size=1.0, scale=50)
trainer.train(n_epochs=250)
batch_entropy = trainer.train_set.entropy_batch_mixing()
full = trainer.create_posterior(vae, gene_dataset, indices=np.arange(len(gene_dataset)))
ll = full.ll(verbose=True)
latent, batch_indices, labels = full.sequential().get_latent()
batch_indices = batch_indices.ravel()
labels = labels.ravel()
res = clustering_scores(np.asarray(latent), labels, 'knn', len(np.unique(labels)))
res["batch_entropy"] = batch_entropy
res["ll"] = ll
print("layer" + str(1))
print(res)

# # TODO: run scVI
# for i in (1, 2, 3):
#     vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches, n_labels=gene_dataset.n_labels,
#               n_hidden=128, n_latent=10, n_layers=1, n_layers_decoder=i, dispersion='gene')
#     trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=1.0)
#     trainer.train(n_epochs=250)
#     batch_entropy = trainer.train_set.entropy_batch_mixing()
#     full = trainer.create_posterior(vae, gene_dataset, indices=np.arange(len(gene_dataset)))
#     ll = full.ll(verbose=True)
#     latent, batch_indices, labels = full.sequential().get_latent()
#     batch_indices = batch_indices.ravel()
#     labels = labels.ravel()
#     res = clustering_scores(np.asarray(latent), labels, 'knn', len(np.unique(labels)))
#     res["batch_entropy"] = batch_entropy
#     res["ll"] = ll
#     print("layer" + str(i))
#     print(res)
