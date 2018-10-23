from scvi.dataset.dataset10X import Dataset10X
from scvi.dataset.pbmc import PbmcDataset
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset
import pandas as pd

from scvi.models.vae import VAE
from scvi.models.scanvi import SCANVI
from scvi.inference import UnsupervisedTrainer, SemiSupervisedTrainer
from sklearn.metrics import roc_auc_score
from scvi.inference.posterior import get_bayes_factors
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset
import sys
from scvi.metrics.clustering import select_indices_evenly
import torch

models = str(sys.argv[1])
use_cuda = True
# pbmc.vae.model.pkl  vae.model.allgenes.pkl  vae.model.rmdeT.pkl


def auc_score_threshold(gene_set, bayes_factor, gene_symbols):
    # put ones on the genes from the gene_set
    true_labels = np.array([g in gene_set for g in gene_symbols])
    estimated_score = np.abs(bayes_factor)
    indices = np.isfinite(estimated_score)
    return roc_auc_score(true_labels[indices], estimated_score[indices])


pbmc = PbmcDataset()
de_data  = pbmc.de_metadata
pbmc.update_cells(pbmc.batch_indices.ravel()==0)
# pbmc.labels = pbmc.labels.reshape(len(pbmc),1)

donor = Dataset10X('fresh_68k_pbmc_donor_a')
donor.gene_names = donor.gene_symbols
import os
if not os.path.isfile('data/10X/fresh_68k_pbmc_donor_a/68k_pbmc_barcodes_annotation.tsv'):
    import urllib.request
    annotation_url = 'https://raw.githubusercontent.com/10XGenomics/single-cell-3prime-paper/master/pbmc68k_analysis/68k_pbmc_barcodes_annotation.tsv'
    urllib.request.urlretrieve(annotation_url, 'data/10X/fresh_68k_pbmc_donor_a/68k_pbmc_barcodes_annotation.tsv')

annotation = pd.read_csv('data/10X/fresh_68k_pbmc_donor_a/68k_pbmc_barcodes_annotation.tsv',sep='\t')
cellid1 = donor.barcodes
temp = cellid1.join(annotation)
assert all(temp[0]==temp['barcodes'])

donor.cell_types,donor.labels = np.unique(temp['celltype'],return_inverse=True)
donor.labels = donor.labels.reshape(len(donor.labels),1)
donor.cell_types = np.array([ 'CD14+ Monocytes','B cells','CD34 cells', 'CD4 T cells','CD4 T cells Regulatory',
                             'CD4 T cells Naive','CD4 Memory T cells','NK cells',
                            'CD8 T cells',  'CD8 T cells Naive', 'Dendritic'])


all_dataset = GeneExpressionDataset.concat_datasets(pbmc, donor)
all_dataset.subsample_genes(5000)
# Now resolve the Gene symbols to properly work with the DE
all_gene_symbols = donor.gene_symbols[
    np.array(
        [np.where(donor.gene_names == x)[0][0] for x in list(all_dataset.gene_names)]
    )]


#####################################################################
# Gene sets
############################################################################
# all_gene_symbols = pbmc.gene_names
# For that, let is import the genesets and intersect with the largest gene set from scRNA-seq
path_geneset = "Additional_Scripts/genesets.txt"
geneset_matrix = np.loadtxt(path_geneset, dtype=np.str)[:, 2:]
CD4_TCELL_VS_BCELL_NAIVE, CD8_TCELL_VS_BCELL_NAIVE, CD8_VS_CD4_NAIVE_TCELL, NAIVE_CD8_TCELL_VS_NKCELL \
    = [set(geneset_matrix[i:i + 2, :].flatten()) & set(all_gene_symbols) for i in [0, 2, 4, 6]]

# these are the length of the positive gene sets for the DE
print((len(CD4_TCELL_VS_BCELL_NAIVE), len(CD8_TCELL_VS_BCELL_NAIVE),
       len(CD8_VS_CD4_NAIVE_TCELL), len(NAIVE_CD8_TCELL_VS_NKCELL)))

print(all_dataset.cell_types)

comparisons = [
    ['CD4 T cells', 'B cells'],
    ['CD8 T cells', 'B cells'],
    ['CD8 T cells', 'CD4 T cells'],
    ['CD8 T cells', 'NK cells']
               ]

# comparisons = [['CD4+ T Helper2','CD19+ B'],
#     ['CD8+ Cytotoxic T', 'CD19+ B'],
#     ['CD8+ Cytotoxic T', 'CD4+ T Helper2'],
#     ['CD8+ Cytotoxic T', 'CD56+ NK']]
#

gene_sets = [CD4_TCELL_VS_BCELL_NAIVE,
             CD8_TCELL_VS_BCELL_NAIVE,
             CD8_VS_CD4_NAIVE_TCELL,
             NAIVE_CD8_TCELL_VS_NKCELL]


path_geneset = "Additional_Scripts/genesets.txt"
geneset_matrix = np.loadtxt(path_geneset, dtype=np.str)[:, 2:]

CD4_TCELL_VS_BCELL_NAIVE, CD8_TCELL_VS_BCELL_NAIVE, CD8_VS_CD4_NAIVE_TCELL, NAIVE_CD8_TCELL_VS_NKCELL, CD8_VS_CD4, B_VS_MDC \
    = [set(geneset_matrix[i:i + 2, :].flatten()) & set(all_gene_symbols) for i in [0, 2, 4, 6, 8, 10]]
  # these are the length of the positive gene sets for the DE

print((len(CD4_TCELL_VS_BCELL_NAIVE), len(CD8_TCELL_VS_BCELL_NAIVE),
       len(CD8_VS_CD4_NAIVE_TCELL), len(NAIVE_CD8_TCELL_VS_NKCELL),
      len(CD8_VS_CD4),len(B_VS_MDC)))

print(all_dataset.cell_types)

gene_sets = [CD4_TCELL_VS_BCELL_NAIVE,
             CD8_TCELL_VS_BCELL_NAIVE,
             CD8_VS_CD4_NAIVE_TCELL,
             NAIVE_CD8_TCELL_VS_NKCELL,
             CD8_VS_CD4,
             B_VS_MDC]


comparisons = [
    ['CD4 T cells', 'B cells'],
    ['CD8 T cells', 'B cells'],
    ['CD8 T cells', 'CD4 T cells'],
    ['CD8 T cells', 'NK cells'],
    ['CD8 T cells', 'CD4 T cells'],
    ['B cells', 'Dendritic Cells']
               ]


print(de_data.columns.values)
CD = de_data['CD_adj.P.Val']
BDC = de_data['BDC_adj.P.Val']
BDC2 = de_data['BDC2_adj.P.Val']
CD = np.asarray(de_data['GS'][CD<0.05])
BDC = np.asarray(de_data['GS'][BDC<0.05])
BDC2 = np.asarray(de_data['GS'][BDC2<0.05])

gene_sets = [set(CD) & set(all_gene_symbols),
             set(BDC)& set(all_gene_symbols),
             set(BDC2) &  set(all_gene_symbols)]

comparisons = [
    ['CD8 T cells', 'CD4 T cells'],
    ['B cells', 'Dendritic Cells'],
    ['B cells', 'Dendritic Cells']
               ]

############################################################################################
# pbmc only
############################################################################################


vae = VAE(pbmc.nb_genes, n_batch=pbmc.n_batches, n_labels=pbmc.n_labels,
          n_hidden=128, n_latent=10, n_layers=1, dispersion='gene')

import torch
trainer = UnsupervisedTrainer(vae, pbmc, train_size=1.0)
# trainer.train(n_epochs=200)
# torch.save(trainer.model,'../DE/pbmc.vae.model.pkl')
trainer.model = torch.load('../DE/pbmc.vae.model.pkl')
full = trainer.create_posterior(trainer.model, pbmc, indices=np.arange(len(pbmc)))
latent, batch_indices, labels = full.sequential().get_latent()
keys = pbmc.cell_types

cell_type_label = [[np.where(pbmc.cell_types == x[i])[0].astype('int')[0] for i in [0, 1]] for x in comparisons]
for t, comparison in enumerate(comparisons):
    # Now for each comparison, let us create a posterior object and compute a Bayes factor
    gene_set = gene_sets[t]
    cell_indices = np.where(np.logical_or(
        labels == cell_type_label[t][0],
        labels == cell_type_label[t][1]))[0]
    de_posterior = trainer.create_posterior(trainer.model, pbmc, indices=cell_indices)
    scale_pbmc = de_posterior.sequential().get_harmonized_scale(0)
    bayes_pbmc = get_bayes_factors(scale_pbmc,
                                   pbmc.labels.ravel()[cell_indices],
                                   cell_type_label[t][0],
                                   cell_type_label[t][1])
    print(auc_score_threshold(gene_set, bayes_pbmc, pbmc.gene_names))

############################################################################################
# all_dataset
############################################################################################



# Start building the models

vae = VAE(all_dataset.nb_genes, n_batch=all_dataset.n_batches, n_labels=all_dataset.n_labels,
          n_hidden=128, n_latent=10, n_layers=2, dispersion='gene')

import torch
trainer = UnsupervisedTrainer(vae, all_dataset, train_size=1.0)
# trainer.train(n_epochs=200)
# torch.save(trainer.model,'../DE/vae.model.rmdeT.pkl')
trainer.model = torch.load('../DE/vae.model.rmdeT.pkl')

trainer.train_set.entropy_batch_mixing()
full = trainer.create_posterior(trainer.model, all_dataset, indices=np.arange(len(all_dataset)))
latent, batch_indices, labels = full.sequential().get_latent()
keys = all_dataset.cell_types


from scvi.inference.posterior import entropy_batch_mixing
sample = select_indices_evenly(2000, batch_indices)
batch_entropy = entropy_batch_mixing(latent[sample, :], batch_indices[sample])

latent_labelled = latent[batch_indices.ravel()==0, :]
latent_unlabelled = latent[batch_indices.ravel()==1, :]
labels_labelled = labels[batch_indices.ravel()==0]
n_labels = np.sum(batch_indices.ravel()==1)
from sklearn.neighbors import KNeighborsClassifier
neigh = KNeighborsClassifier(n_neighbors=10)
neigh = neigh.fit(latent_labelled, labels_labelled)

pred = neigh.predict(latent_labelled)
np.mean(pred==labels_labelled)

pred = neigh.predict(latent_unlabelled)
pred = np.concatenate([labels_labelled,pred])

# from scvi.metrics.clustering import clustering_scores
# clustering_scores(np.asarray(latent), labels, 'knn')
# batch_indices = batch_indices.ravel()
#
#
# scanvi = SCANVI(all_dataset.nb_genes, all_dataset.n_batches, all_dataset.n_labels, n_layers=2)
# scanvi.load_state_dict(vae.state_dict(), strict=False)
# trainer_scanvi = SemiSupervisedTrainer(scanvi, all_dataset, classification_ratio=50,
#                                        n_epochs_classifier=1, lr_classification=5 * 1e-3)
# trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=(all_dataset.batch_indices.ravel() ==0))
# trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(
#     indices=(all_dataset.batch_indices.ravel() == 1)
# )
# trainer_scanvi.train(n_epochs=50)
# full = trainer_scanvi.create_posterior(scanvi, all_dataset, indices=np.arange(len(all_dataset)))
# latent, batch_indices, labels = full.sequential().get_latent()
# sample = select_indices_evenly(2000, batch_indices)
# batch_entropy = entropy_batch_mixing(latent[sample, :], batch_indices[sample])
#
# # import torch
# torch.save(trainer.model,'../DE/vae.model.pkl')
# torch.save(trainer_scanvi.model,'../DE/scanvi.model.pkl')
#
# keys = all_dataset.cell_types
# latent, batch_indices, labels = trainer_scanvi.labelled_set.get_latent()
# pred = trainer_scanvi.unlabelled_set.compute_predictions()
# np.mean(pred[0] == pred[1])
# trainer_scanvi.full_dataset.entropy_batch_mixing()
# # 0.6
#
# scanvi_posterior = trainer_scanvi.create_posterior(trainer_scanvi.model, all_dataset)
# # Extract the predicted labels from SCANVI
#
# pred = scanvi_posterior.compute_predictions()[1]
# batch_indices = all_dataset.batch_indices.ravel()
# np.mean(pred[batch_indices==0]==labels_labelled)
# # 0.84 instead of 0.96 with nn


cell_type_label = [[np.where(all_dataset.cell_types == x[i])[0].astype('int')[0] for i in [0, 1]] for x in comparisons]
for t, comparison in enumerate(comparisons):
    # Now for each comparison, let us create a posterior object and compute a Bayes factor
    gene_set = gene_sets[t]
    cell_indices = np.where(np.logical_or(
        pred == cell_type_label[t][0],
        pred == cell_type_label[t][1]))[0]
    de_posterior = trainer.create_posterior(trainer.model, all_dataset, indices=cell_indices)
    scale_pbmc = de_posterior.sequential().get_harmonized_scale(0)
    scale_68k = de_posterior.sequential().get_harmonized_scale(1)
    # For Chenling: I looked again at the number of cells,
    # if we use all of them, we are OK using just one sample from the posterior
    # first grab the original bayes factor by ignoring the unlabeled cells
    bayes_pbmc = get_bayes_factors(scale_pbmc,
                                   all_dataset.labels.ravel()[cell_indices],
                                   cell_type_label[t][0],
                                   cell_type_label[t][1])
    # second get them for all the predicted labels cross-datasets
    probs_all_imputed_pbmc = get_bayes_factors(scale_pbmc,
                                               pred[cell_indices],
                                               cell_type_label[t][0],
                                               cell_type_label[t][1], logit=False)
    probs_all_imputed_68k = get_bayes_factors(scale_68k,
                                              pred[cell_indices],
                                              cell_type_label[t][0],
                                              cell_type_label[t][1], logit=False)
    p_s = pbmc.labels.shape[0] / all_dataset.labels.shape[0]
    bayes_all_imputed = p_s * probs_all_imputed_pbmc + (1 - p_s) * probs_all_imputed_68k
    bayes_all_imputed = np.log(bayes_all_imputed + 1e-8) - np.log(1 - bayes_all_imputed + 1e-8)
    # COMMENT FOR CHENLING: the bayes variables (bayes_pbmc and bayes_all_imputed) should now span all values
    # (not only positive or negative).
    # Can you check that with np.min and np.max ?
    print(auc_score_threshold(gene_set, bayes_pbmc, all_gene_symbols))
    print(auc_score_threshold(gene_set, bayes_all_imputed, all_gene_symbols))

