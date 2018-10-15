from scvi.dataset.dataset10X import Dataset10X
from scvi.dataset.pbmc import PbmcDataset
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset

from scvi.models.vae import VAE
from scvi.models.scanvi import SCANVI
from scvi.inference import UnsupervisedTrainer, SemiSupervisedTrainer
from sklearn.metrics import roc_auc_score
from scvi.inference.posterior import get_bayes_factors
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset
import sys
from scvi.metrics.clustering import select_indices_evenly

models = str(sys.argv[1])
plotname = 'Easy1'
use_cuda = True


def auc_score_threshold(gene_set, bayes_factor, gene_symbols):
    # put ones on the genes from the gene_set
    true_labels = np.array([g in gene_set for g in gene_symbols])
    estimated_score = np.abs(bayes_factor)
    indices = np.isfinite(estimated_score)
    return roc_auc_score(true_labels[indices], estimated_score[indices])


# We need to modify this import to get all the genes
pbmc = PbmcDataset()
pbmc.update_cells(pbmc.batch_indices.ravel()==0)

from scvi.harmonization.utils_chenling import get_matrix_from_dir,assign_label
count, geneid, cellid = get_matrix_from_dir('cite')
count = count.T.tocsr()
seurat = np.genfromtxt('../cite/cite.seurat.labels', dtype='str', delimiter=',')
cellid = np.asarray([x.split('-')[0] for x in cellid])
labels_map = [0, 0, 1, 2, 3, 4, 5, 6]
labels = seurat[1:, 4]
cell_type = ["CD4+ T Helper2", "CD56+ NK", "CD14+ Monocyte", "CD19+ B", "CD8+ Cytotoxic T", "FCGR3A Monocyte", "na"]
dataset2 = assign_label(cellid, geneid, labels_map, count, cell_type, seurat)
rmCellTypes = {'na', 'dendritic'}
newCellType = [k for i, k in enumerate(dataset2.cell_types) if k not in rmCellTypes]
dataset2.filter_cell_types(newCellType)
pbmc = dataset2

pbmc68k = Dataset10X('fresh_68k_pbmc_donor_a')
pbmc68k.cell_types = ['unlabelled']
pbmc68k.labels = np.repeat(0, len(pbmc68k)).reshape(len(pbmc68k), 1)
pbmc68k.gene_names = pbmc68k.gene_symbols

all_dataset = GeneExpressionDataset.concat_datasets(pbmc, pbmc68k)
all_dataset.subsample_genes(5000)
# Now resolve the Gene symbols to properly work with the DE
all_gene_symbols = pbmc68k.gene_symbols[
    np.array(
        [np.where(pbmc68k.gene_names == x)[0][0] for x in list(all_dataset.gene_names)]
    )]

# For that, let is import the genesets and intersect with the largest gene set from scRNA-seq
path_geneset = "Additional_Scripts/genesets.txt"
geneset_matrix = np.loadtxt(path_geneset, dtype=np.str)[:, 2:]
CD4_TCELL_VS_BCELL_NAIVE, CD8_TCELL_VS_BCELL_NAIVE, CD8_VS_CD4_NAIVE_TCELL, NAIVE_CD8_TCELL_VS_NKCELL \
    = [set(geneset_matrix[i:i + 2, :].flatten()) & set(all_gene_symbols) for i in [0, 2, 4, 6]]

# these are the length of the positive gene sets for the DE
print((len(CD4_TCELL_VS_BCELL_NAIVE), len(CD8_TCELL_VS_BCELL_NAIVE),
       len(CD8_VS_CD4_NAIVE_TCELL), len(NAIVE_CD8_TCELL_VS_NKCELL)))

print(all_dataset.cell_types)

# Start building the models

vae = VAE(all_dataset.nb_genes, n_batch=all_dataset.n_batches, n_labels=all_dataset.n_labels,
          n_hidden=128, n_latent=10, n_layers=2, dispersion='gene')

trainer = UnsupervisedTrainer(vae, all_dataset, train_size=1.0)
trainer.train(n_epochs=100)
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

# pred = scanvi_posterior.compute_predictions()[1]

# comparisons = [
#     ['CD4 T cells', 'B cells'],
#     ['CD8 T cells', 'B cells'],
#     ['CD8 T cells', 'CD4 T cells'],
#     ['CD8 T cells', 'NK cells']
#                ]

comparisons = [['CD4+ T Helper2','CD19+ B'],
    ['CD8+ Cytotoxic T', 'CD19+ B'],
    ['CD8+ Cytotoxic T', 'CD4+ T Helper2'],
    ['CD8+ Cytotoxic T', 'CD56+ NK']]


gene_sets = [CD4_TCELL_VS_BCELL_NAIVE,
             CD8_TCELL_VS_BCELL_NAIVE,
             CD8_VS_CD4_NAIVE_TCELL,
             NAIVE_CD8_TCELL_VS_NKCELL]

cell_type_label = [[np.where(all_dataset.cell_types == x[i])[0].astype('int')[0] for i in [0, 1]] for x in comparisons]
for t, comparison in enumerate(comparisons):
    # Now for each comparison, let us create a posterior object and compute a Bayes factor
    gene_set = gene_sets[0]
    cell_indices = np.where(np.logical_or(
        pred == cell_type_label[t][0],
        pred == cell_type_label[t][1]))[0]
    de_posterior = trainer.create_posterior(vae, all_dataset, indices=cell_indices)
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

