use_cuda = True
from scvi.dataset.dataset10X import Dataset10X
from scvi.dataset.pbmc import PbmcDataset
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import eval_latent, run_model
import sys
from scvi.models.vae import VAE
from scvi.models.scanvi import SCANVI
from scvi.inference import UnsupervisedTrainer, SemiSupervisedTrainer

# We need to modify this import to get all the genes
pbmc = PbmcDataset(filter_out_de_genes=False, use_symbols=False)
pbmc68k = Dataset10X('fresh_68k_pbmc_donor_a')

pbmc68k.cell_types = ['unlabelled']
pbmc68k.labels = np.repeat(0, len(pbmc68k)).reshape(len(pbmc68k), 1)

all_dataset = GeneExpressionDataset.concat_datasets(pbmc, pbmc68k)
all_dataset.subsample_genes(5000)
# Now resolve the Gene symbols to properly work with the DE
all_gene_symbols = pbmc68k.gene_symbols[
    np.array(
        [np.where(pbmc68k.gene_names == x)[0][0] for x in list(all_dataset.gene_names)]
    )]

# Need to filter the genes for the DE before subsampling
# ALL GENE SET REFERENCES COME FROM GSE22886
# [CD4_TCELL_VS_BCELL_NAIVE_UP, CD4_TCELL_VS_BCELL_NAIVE_DN
# CD8_TCELL_VS_BCELL_NAIVE_UP, CD8_TCELL_VS_BCELL_NAIVE_DN
# CD8_VS_CD4_NAIVE_TCELL_UP, CD8_VS_CD4_NAIVE_TCELL_DN
# NAIVE_CD8_TCELL_VS_NKCELL_UP, NAIVE_CD8_TCELL_VS_NKCELL_DN]

# For that, let is import the genesets and intersect with the largest gene set from scRNA-seq
path_geneset = "../Additional_Scripts/genesets.txt"
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
trainer.train(n_epochs=250)

scanvi = SCANVI(all_dataset.nb_genes, all_dataset.n_batches, all_dataset.n_labels, n_layers=2)
scanvi.load_state_dict(vae.state_dict(), strict=False)
trainer_scanvi = SemiSupervisedTrainer(scanvi, all_dataset, classification_ratio=1,
                                       n_epochs_classifier=1, lr_classification=5 * 1e-3)
trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=(all_dataset.batch_indices != 2))
trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(
    indices=(all_dataset.batch_indices == 2)
)
trainer_scanvi.train(n_epochs=50)

keys = all_dataset.cell_types
latent, batch_indices, labels = trainer_scanvi.labelled_set.get_latent()
pred = trainer_scanvi.labelled_set.compute_predictions()
np.mean(pred[0] == pred[1])

scanvi_posterior = trainer_scanvi.create_posterior(trainer_scanvi.model, all_dataset)
# Extract the predicted labels from SCANVI
pred = scanvi_posterior.compute_predictions()[1]

# comparison = ['CD4 T cells', 'CD8 T cells']
comparison = ['CD8 T cells', 'B cells']
# Now for each comparison, let us create a posterior object and compute a Bayes factor
cell_type_label = [np.where(all_dataset.cell_types == comparison[i])[0][0] for i in [0, 1]]
cell_indices = np.where(np.logical_or(pred == cell_type_label[0], pred == cell_type_label[1]))[0]
de_posterior = trainer.create_posterior(vae, all_dataset, indices=cell_indices)

# TODO: clarify how many batches are used in the experiment.
# TODO: 0 and 1 might be from the same Michael PBMC dataset
scale_pbmc = de_posterior.sequential().get_harmonized_scale(0)
scale_68k = de_posterior.sequential().get_harmonized_scale(1)

# For Chenling: I looked again at the number of cells,
# if we use all of them, we are OK using just one sample from the posterior

from scvi.inference.posterior import get_bayes_factors

# first grab the original bayes factor by ignoring the unlabeled cells
bayes_pbmc = get_bayes_factors(scale_pbmc,
                               all_dataset.labels.ravel()[cell_indices],
                               cell_type_label[0],
                               cell_type_label[1])

# second get them for all the predicted labels cross-datasets
bayes_all_imputed_pbmc = get_bayes_factors(scale_pbmc,
                                           pred[cell_indices],
                                           cell_type_label[0],
                                           cell_type_label[1])
bayes_all_imputed_68k = get_bayes_factors(scale_68k,
                                          pred[cell_indices],
                                          cell_type_label[0],
                                          cell_type_label[1])

# TODO: We need to average the Bayes Factor in the exponential natural parametrization

# TODO: We need to compute the AUC scores
