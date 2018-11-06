import sys
from scvi.dataset.dataset10X import Dataset10X
from scvi.dataset.pbmc import PbmcDataset
import pandas as pd
import matplotlib.pyplot as plt

from scvi.models.vae import VAE
from scvi.models.scanvi import SCANVI
from scvi.inference import UnsupervisedTrainer, SemiSupervisedTrainer
from sklearn.metrics import roc_auc_score
from scvi.inference.posterior import get_bayes_factors
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset
import os
from scvi.metrics.clustering import select_indices_evenly


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



subsampling_rate = float(sys.argv[1])
print("starting subsampling with rate" + str(subsampling_rate))

pbmc = PbmcDataset()
de_data  = pbmc.de_metadata
pbmc.update_cells(pbmc.batch_indices.ravel()==0)
# pbmc.labels = pbmc.labels.reshape(len(pbmc),1)

donor = Dataset10X('fresh_68k_pbmc_donor_a')
donor.gene_names = donor.gene_symbols
donor.labels = np.repeat(0,len(donor)).reshape(len(donor),1)
donor.cell_types = ['unlabelled']
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


vae = VAE(all_dataset.nb_genes, n_batch=all_dataset.n_batches, n_labels=all_dataset.n_labels,
          n_hidden=128, n_latent=10, n_layers=2, dispersion='gene')

import torch
trainer = UnsupervisedTrainer(vae, all_dataset, train_size=1.0)
# trainer.train(n_epochs=200)
# torch.save(trainer.model,'../DE/vae.model.pkl')
trainer.model = torch.load('DE/vae.model.pkl')

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


#######
#
#  CORRUPT THE DATASET
#
######


print("SUBSAMPLING DATA: "+str(subsampling_rate)) 

A_info = list(np.unique(all_dataset.labels.ravel(), return_counts=True))
A_info[0] = all_dataset.cell_types[A_info[0]]
AB_info = list(np.unique(vae_pred, return_counts=True))
AB_info[0] = all_dataset.cell_types[AB_info[0]]
A_info, AB_info

# JUST SUBSAMPLE THE CD4 AND CD8

type_set = [np.where(all_dataset.cell_types == comparisons[0][i])[0].astype('int')[0] for i in [0, 1]]
cell_set = np.array([False] * all_dataset.labels.ravel().shape[0])
for x in np.unique(all_dataset.labels.ravel()):
    cell_idx = np.where(vae_pred == x)[0]
    if x in type_set:
        cell_idx = np.random.choice(cell_idx, size=int(subsampling_rate * cell_idx.shape[0]))
    print(all_dataset.cell_types[x], " from ", np.sum(vae_pred == x), " cells to ", cell_idx.shape[0], " cells")
    cell_set[cell_idx] = True
   

all_dataset.update_cells(cell_set)
vae_pred = vae_pred[cell_set]
batch_indices = all_dataset.batch_indices

print("CD8 subsampled to"+ str(np.sum(vae_pred == type_set[0])/128.9)+" percent of the original data")

trainer = UnsupervisedTrainer(vae, all_dataset, train_size=1.0)
trainer.train(n_epochs=200)

#######
#
#  DIFFERENTIAL EXPRESSION CODE
#
######

from copy import deepcopy
batch2 = deepcopy(all_dataset)
batch2.update_cells(batch_indices.ravel()==1)
cell_type_label = \
    [[np.where(all_dataset.cell_types == x[i])[0].astype('int')[0] for i in [0, 1]] for x in comparisons]
    
from scipy.stats import kendalltau

import rpy2
from rpy2.robjects import r
import rpy2.robjects as robj
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()
from rpy2.robjects.lib import grid
from rpy2.robjects import r, Formula
py2ri_orig = rpy2.robjects.conversion.py2ri
from rpy2.rinterface import RRuntimeWarning

r["library"]("idr")

def idr(bayes, p_value, p_prior=0.1):
    p_val_1r = r.matrix(bayes[:, np.newaxis], nrow=3343, ncol=1)
    r.assign("p_val_1", p_val_1r)

    p_val_2r = r.matrix(p_value[:, np.newaxis], nrow=3343, ncol=1)
    r.assign("p_val_2", p_val_2r)

    r("x <- cbind(p_val_1[, 1], p_val_2[, 1])")
    r("mu = 1")
    r("sigma = 0.5")
    r("rho = 0.5")
    r.assign("p", 0.25)
    return kendalltau(bayes, p_value)[0]

    r("idr.out <- est.IDR(x, mu, sigma, rho, p, eps=0.001, max.ite=20)")
    return r("idr.out$para$p")[0]

r["library"]("edgeR")

import pandas 
def conversion_pydataframe(obj):
    """
    Convert pandas DataFrame or python object to an R dataframe/object.
    """
    if isinstance(obj, pandas.core.frame.DataFrame):
        od = OrderedDict()
        for name, values in obj.iteritems():
            if values.dtype.kind == 'O':
                od[name] = rpy2.robjects.vectors.StrVector(values)
            else:
                od[name] = rpy2.robjects.conversion.py2ri(values)
        return rpy2.robjects.vectors.DataFrame(od)
    else:
        return py2ri_orig(obj)
    
def run_edgeR(gene_expression, bio_assignment, gene_names, batch_info=None, batch=True):
    if batch_info is None:
        batch = False
    r_counts = conversion_pydataframe(gene_expression)
    r_bio_group = conversion_pydataframe(bio_assignment)
    r_dge = r.DGEList(counts=r.t(r_counts), genes=gene_names)
    r.assign("dge", r_dge)
    r.assign("bio_group", r.factor(r_bio_group))
    r("dge$samples$bio_group <- bio_group")

    if batch:
        r_batch_group = conversion_pydataframe(batch_info)
        r.assign("batch_group", r.factor(r_batch_group))
        r("dge$samples$batch_group <- batch_group")

    r("""dge <- suppressWarnings(edgeR::calcNormFactors(dge))""")
    
    if not batch:
        r("""design <- model.matrix(~bio_group, data = dge$samples)""")
        r("""colnames(design) <- c("Intercept", "bio")""")

    if batch:
        r("""design <- model.matrix(~bio_group+batch_group, data = dge$samples)""")
        r("""colnames(design) <- c("Intercept", "bio", "batch")""")

    r("""dge <- estimateDisp(dge, design)""")
    
    r("""fit <- glmFit(dge, design)""")
    if not batch:
        r("""lrt <- glmLRT(fit)""")
    if batch:
        r("""lrt <- glmLRT(fit, coef="bio")""")
    return r("lrt$table$PValue")


interest_list = ["CD", "BDC"]
p_prior_list = [0.03, 0.25]
de_data.index = de_data["GS"]

# TODO: CHANGE SAMPLING
def bayes(t, comparaison, T=1, n_perm = 10000):
    gene_set = gene_sets[t]
    
    cell_idx_8k = np.where(np.logical_or(
        all_dataset.labels.ravel() == cell_type_label[t][0],
        all_dataset.labels.ravel() == cell_type_label[t][1]))[0]
    
    cell_idx_68k = np.where(np.logical_or(
        vae_pred[batch_indices.ravel()==1] == cell_type_label[t][0],
        vae_pred[batch_indices.ravel()==1] == cell_type_label[t][1]))[0]
    
    cell_indices = np.where(np.logical_or(
        vae_pred == cell_type_label[t][0],
        vae_pred == cell_type_label[t][1]))[0]

    joint_de_posterior = trainer.create_posterior(trainer.model, all_dataset, indices=cell_indices)
    scale_pbmc = joint_de_posterior.sequential().get_harmonized_scale(0)
    scale_68k = joint_de_posterior.sequential().get_harmonized_scale(1)
    questionable_de_posterior = trainer.create_posterior(trainer.model, batch2, indices=cell_idx_68k)
    questionable_scale_68k = questionable_de_posterior.sequential().get_harmonized_scale(1)
    
    res_vi = np.zeros((3, 2, T)) # 3 datasets, 2 metrics, T indep runs
    res_eR = np.zeros((3, 2, T))
    p_value = de_data[interest_list[t] + "_adj.P.Val"][all_gene_symbols].values
    p_prior = p_prior_list[t]

    for rep in range(T):
        
        #PBMC8K only
        bayes_pbmc = get_bayes_factors(scale_pbmc,
                                       all_dataset.labels.ravel()[cell_indices],
                                       cell_type_label[t][0],
                                       cell_type_label[t][1], m_permutation=n_perm)
            
        res_vi[0, 0, rep] = auc_score_threshold(gene_set, bayes_pbmc, all_gene_symbols)
        res_vi[0, 1, rep] = idr(np.abs(bayes_pbmc), -np.log(p_value), p_prior=p_prior)
        
        ind_0 = np.random.choice(np.where(all_dataset.labels.ravel() == cell_type_label[t][0])[0], 100)
        ind_1 = np.random.choice(np.where(all_dataset.labels.ravel() == cell_type_label[t][1])[0], 100)
        expression_data = np.vstack((all_dataset.X[ind_0].A, all_dataset.X[ind_1].A))
        bio_data = np.hstack((all_dataset.labels.ravel()[ind_0], all_dataset.labels.ravel()[ind_1]))
        edgeR_pbmc = run_edgeR(expression_data, bio_data, all_dataset.gene_names)
        
        res_eR[0, 0, rep] = auc_score_threshold(gene_set, -np.log(edgeR_pbmc), all_gene_symbols)
        res_eR[0, 1, rep] = idr(-np.log(edgeR_pbmc), -np.log(p_value), p_prior=p_prior)
        
        # PBMC68K only        
        bayes_questionable = get_bayes_factors(questionable_scale_68k,
                                       vae_pred[batch_indices.ravel()==1][cell_idx_68k],
                                       cell_type_label[t][0],
                                       cell_type_label[t][1], logit=True, m_permutation=n_perm)
        
        res_vi[1, 0, rep] = auc_score_threshold(gene_set, bayes_questionable, all_gene_symbols)
        res_vi[1, 1, rep] = idr(np.abs(bayes_questionable), -np.log(p_value), p_prior=p_prior)
        
        ind_0 = np.random.choice(np.where(vae_pred[batch_indices.ravel()==1] == cell_type_label[t][0])[0], 100)
        ind_1 = np.random.choice(np.where(vae_pred[batch_indices.ravel()==1] == cell_type_label[t][1])[0], 100)
        expression_data = np.vstack((all_dataset.X[batch_indices.ravel()==1][ind_0].A, 
                                     all_dataset.X[batch_indices.ravel()==1][ind_1].A))
        bio_data = np.hstack((vae_pred[batch_indices.ravel()==1][ind_0], 
                                vae_pred[batch_indices.ravel()==1][ind_1]))
        edgeR_68k = run_edgeR(expression_data, bio_data, all_dataset.gene_names)
        
        
        res_eR[1, 0, rep] = auc_score_threshold(gene_set, -np.log(edgeR_68k), all_gene_symbols)
        res_eR[1, 1, rep] = idr(-np.log(edgeR_68k), -np.log(p_value), p_prior=p_prior)
        
        
        #WHOLE PBMC
        probs_all_imputed_pbmc = get_bayes_factors(scale_pbmc,
                                                   vae_pred[cell_indices],
                                                   cell_type_label[t][0],
                                                   cell_type_label[t][1], logit=False, m_permutation=n_perm)
        probs_all_imputed_68k = get_bayes_factors(scale_68k,
                                                  vae_pred[cell_indices],
                                                  cell_type_label[t][0],
                                                  cell_type_label[t][1], logit=False, m_permutation=n_perm)

        p_s = 0.5
        bayes_all_imputed = p_s * probs_all_imputed_pbmc + (1 - p_s) * probs_all_imputed_68k
        bayes_all_imputed = np.log(bayes_all_imputed + 1e-8) - np.log(1 - bayes_all_imputed + 1e-8)
        
        res_vi[2, 0, rep] = auc_score_threshold(gene_set, bayes_all_imputed, all_gene_symbols)
        res_vi[2, 1, rep] = idr(np.abs(bayes_all_imputed), -np.log(p_value), p_prior=p_prior)
        
        
        
        ind_0 = np.random.choice(np.where(vae_pred == cell_type_label[t][0])[0], 100)
        ind_1 = np.random.choice(np.where(vae_pred == cell_type_label[t][1])[0], 100)
        expression_data = np.vstack((all_dataset.X[ind_0].A, 
                                     all_dataset.X[ind_1].A))
        bio_data = np.hstack((vae_pred[ind_0], 
                                vae_pred[ind_1]))
        batch_data = np.hstack((batch_indices.ravel()[ind_0], 
                                batch_indices.ravel()[ind_1]))
        edgeR_all = run_edgeR(expression_data, bio_data, all_dataset.gene_names, batch_info=batch_data)
        
        res_eR[2, 0, rep] = auc_score_threshold(gene_set, -np.log(edgeR_all), all_gene_symbols)
        res_eR[2, 1, rep] = idr(-np.log(edgeR_all), -np.log(p_value), p_prior=p_prior)
            
    return res_vi, res_eR

for t, comparison in enumerate(comparisons):
    print(t, comparison)
    
res_vi_CD, res_eR_CD = bayes(0, ['CD8 T cells', 'CD4 T cells'], T=10)

print("\n\n\n\ RESULTS:     \n")
print("vi")
print(res_vi_CD.__repr__())
print("eR")
print(res_eR_CD.__repr__())

