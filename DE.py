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


def WeightedAccuracy(y,y_pred,cell_types):
    res = dict()
    for i in np.unique(y):
        res[cell_types[i]] = (np.mean(y_pred[y == i] == i), sum(y==i))
    return(res)



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
donor.cell_types = np.array([ 'CD14+ Monocytes','B cells','CD34 cells', 'CD4 T cells 2','CD4 T cells Regulatory',
                             'CD4 T cells Naive','CD4 T cells Memory','NK cells',
                            'CD8 T cells',  'CD8 T cells Naive', 'Dendritic Cells'])


all_dataset = GeneExpressionDataset.concat_datasets(pbmc, donor)
# Now resolve the Gene symbols to properly work with the DE
all_gene_symbols = donor.gene_symbols[
    np.array(
        [np.where(donor.gene_names == x)[0][0] for x in list(all_dataset.gene_names)]
    )]


#####################################################################
# Gene sets 1
############################################################################
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


gene_sets = [CD4_TCELL_VS_BCELL_NAIVE,
             CD8_TCELL_VS_BCELL_NAIVE,
             CD8_VS_CD4_NAIVE_TCELL,
             NAIVE_CD8_TCELL_VS_NKCELL]

#####################################################################
# Gene sets 2
############################################################################
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
torch.save(trainer.model,'../DE/pbmc.vae.model.pkl')
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
# VAE
############################################################################################
vae = VAE(all_dataset.nb_genes, n_batch=all_dataset.n_batches, n_labels=all_dataset.n_labels,
          n_hidden=128, n_latent=10, n_layers=2, dispersion='gene')

import torch
trainer = UnsupervisedTrainer(vae, all_dataset, train_size=1.0)
# trainer.train(n_epochs=200)
# torch.save(trainer.model,'../DE/vae.model.pkl')
trainer.model = torch.load('../DE/vae.model.pkl')

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
labels_unlabelled = labels[batch_indices.ravel()==1]
n_labels = np.sum(batch_indices.ravel()==1)
from sklearn.neighbors import KNeighborsClassifier
neigh = KNeighborsClassifier(n_neighbors=10)
neigh = neigh.fit(latent_labelled, labels_labelled)

vae_pred = neigh.predict(latent)
np.mean(vae_pred[batch_indices.ravel()==0]==labels[batch_indices.ravel()==0])
np.mean(vae_pred[batch_indices.ravel()==1]==labels[batch_indices.ravel()==1])
############################################################################################
#  SCANVI
############################################################################################
scanvi = SCANVI(all_dataset.nb_genes, all_dataset.n_batches, all_dataset.n_labels, n_layers=2)
scanvi.load_state_dict(trainer.model.state_dict(), strict=False)
trainer_scanvi = SemiSupervisedTrainer(scanvi, all_dataset, classification_ratio=50,
                                       n_epochs_classifier=1, lr_classification=5 * 1e-3)

trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(
    indices=(all_dataset.batch_indices.ravel() == 0)
)
trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(
    indices=(all_dataset.batch_indices.ravel() == 1)
)

# trainer_scanvi.train(n_epochs=50)
# torch.save(trainer_scanvi.model,'../DE/scanvi.model.pkl')
trainer_scanvi.model = torch.load('../DE/scanvi.model.pkl')
full = trainer_scanvi.create_posterior(trainer_scanvi.model, all_dataset, indices=np.arange(len(all_dataset)))
latent, batch_indices, labels = full.sequential().get_latent()
keys = all_dataset.cell_types

latent_labelled = latent[batch_indices.ravel()==0, :]
latent_unlabelled = latent[batch_indices.ravel()==1, :]
labels_labelled = labels[batch_indices.ravel()==0]
labels_unlabelled = labels[batch_indices.ravel()==1]
n_labels = np.sum(batch_indices.ravel()==1)
neigh = KNeighborsClassifier(n_neighbors=10)
neigh = neigh.fit(latent_labelled, labels_labelled)
vae_pred = neigh.predict(latent)
np.mean(vae_pred[batch_indices.ravel()==0]==labels[batch_indices.ravel()==0])
np.mean(vae_pred[batch_indices.ravel()==1]==labels[batch_indices.ravel()==1])

pred = full.sequential().compute_predictions()
np.mean(pred[1][batch_indices.ravel()==0]==labels[batch_indices.ravel()==0])
np.mean(pred[1][batch_indices.ravel()==1]==labels[batch_indices.ravel()==1])

vae_pred = neigh.predict(latent)
np.mean(vae_pred[batch_indices.ravel()==0]==labels[batch_indices.ravel()==0])
np.mean(vae_pred[batch_indices.ravel()==1]==labels[batch_indices.ravel()==1])

vae_pred = np.concatenate([all_dataset.labels.ravel()[batch_indices.ravel()==0],vae_pred[batch_indices.ravel()==1]])
#
# # #todo map the other labels to the same labels as pbmc8k labels
# i = np.arange(len(keys))[keys=='CD8 T cells Naive']
# j = np.arange(len(keys))[keys=='CD8 T cells']
# pred[0][pred[0]==i]=j
#
# for CD4sub in ['CD4 T cells Naive','CD4 T cells Regulatory', 'CD4 T cells Memory', 'CD4 T cells 2']:
#     i = np.arange(len(keys))[keys == CD4sub]
#     j = np.arange(len(keys))[keys == 'CD4 T cells']
#     pred[0][pred[0] == i] = j
#     pred[1][pred[1]==i]=j
#
# print(keys[np.unique(pred[0][batch_indices.ravel()==1])])
# print(keys[np.unique(pred[1][batch_indices.ravel()==1])])
#
#
# from sklearn.metrics import confusion_matrix
# shared = set(keys[np.unique(pred[0][batch_indices.ravel()==1])]).intersection(set(keys[np.unique(pred[0][batch_indices.ravel()==0])]))
# shared = list(shared)
# CM_VAE = confusion_matrix(keys[pred[0][batch_indices.ravel()==1]],keys[vae_pred[batch_indices.ravel()==1]],labels=shared)
# CM_VAE = pd.DataFrame(CM_VAE,index=shared,columns=shared)
# WeightedAccuracy(pred[0][batch_indices.ravel()==0],vae_pred[batch_indices.ravel()==0],all_dataset.cell_types)
# WeightedAccuracy(pred[0][batch_indices.ravel()==1],vae_pred[batch_indices.ravel()==1],all_dataset.cell_types)
#
# CM_SCANVI = confusion_matrix(keys[pred[0][batch_indices.ravel()==1]],keys[pred[1][batch_indices.ravel()==1]],labels=shared)
# CM_SCANVI = pd.DataFrame(CM_SCANVI,index=shared,columns=shared)
# WeightedAccuracy(pred[0][batch_indices.ravel()==0],pred[1][batch_indices.ravel()==0],all_dataset.cell_types)
# WeightedAccuracy(pred[0][batch_indices.ravel()==1],pred[1][batch_indices.ravel()==1],all_dataset.cell_types)
#
#
# sample = select_indices_evenly(2000, batch_indices)
# batch_entropy = entropy_batch_mixing(latent[sample, :], batch_indices[sample])
#
#
#

###########################################################################################
# Computing AUC
############################################################################################
from copy import deepcopy
batch2 = deepcopy(all_dataset)
batch2.update_cells(batch_indices.ravel()==1)
cell_type_label = [[np.where(all_dataset.cell_types == x[i])[0].astype('int')[0] for i in [0, 1]] for x in comparisons]
for t, comparison in enumerate(comparisons):
    # Now for each comparison, let us create a posterior object and compute a Bayes factor
    # t=0
    # comparison = comparisons[0]
    # gene_set = gene_sets[t]
    cell_indices = np.where(np.logical_or(
        vae_pred == cell_type_label[t][0],
        vae_pred == cell_type_label[t][1]))[0]
    cell_idx_68k = np.where(np.logical_or(
        vae_pred[batch_indices.ravel()==1] == cell_type_label[t][0],
        vae_pred[batch_indices.ravel()==1] == cell_type_label[t][1]))[0]
    joint_de_posterior = trainer.create_posterior(trainer.model, all_dataset, indices=cell_indices)
    scale_pbmc = joint_de_posterior.sequential().get_harmonized_scale(0)
    scale_68k = joint_de_posterior.sequential().get_harmonized_scale(1)
    questionable_de_posterior = trainer.create_posterior(trainer.model, batch2,indices=cell_idx_68k)
    questionable_scale_68k = questionable_de_posterior.sequential().get_harmonized_scale(0)
    # For Chenling: I looked again at the number of cells,
    # if we use all of them, we are OK using just one sample from the posterior
    # first grab the original bayes factor by ignoring the unlabeled cells
    res_pbmc=[]
    res_all = []
    res_questionable = []
    for rep in range(10):
        bayes_pbmc = get_bayes_factors(scale_pbmc,
                                       all_dataset.labels.ravel()[cell_indices],
                                       cell_type_label[t][0],
                                       cell_type_label[t][1])
        # second get them for all the predicted labels cross-datasets
        probs_all_imputed_pbmc = get_bayes_factors(scale_pbmc,
                                                   vae_pred[cell_indices],
                                                   cell_type_label[t][0],
                                                   cell_type_label[t][1], logit=False)
        probs_all_imputed_68k = get_bayes_factors(scale_68k,
                                                  vae_pred[cell_indices],
                                                  cell_type_label[t][0],
                                                  cell_type_label[t][1], logit=False)
        bayes_questionable = get_bayes_factors(questionable_scale_68k,
                                               vae_pred[batch_indices.ravel()==1][cell_idx_68k],
                                               cell_type_label[t][0],
                                               cell_type_label[t][1], logit=True)
        p_s = pbmc.labels.shape[0] / all_dataset.labels.shape[0]
        bayes_all_imputed = p_s * probs_all_imputed_pbmc + (1 - p_s) * probs_all_imputed_68k
        bayes_all_imputed = np.log(bayes_all_imputed + 1e-8) - np.log(1 - bayes_all_imputed + 1e-8)
        # COMMENT FOR CHENLING: the bayes variables (bayes_pbmc and bayes_all_imputed) should now span all values
        # (not only positive or negative).
        # Can you check that with np.min and np.max ?
        res_pbmc.append(auc_score_threshold(gene_set, bayes_pbmc, all_gene_symbols))
        res_all.append(auc_score_threshold(gene_set, bayes_all_imputed, all_gene_symbols))
        res_questionable.append(auc_score_threshold(gene_set, bayes_questionable, all_gene_symbols))
    print(str(t))
    print(res_pbmc)
    print(res_all)
    print(res_questionable)

