use_cuda = True
from scvi.harmonization.utils_chenling import get_matrix_from_dir
from scvi.harmonization.benchmark import assign_label
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import eval_latent, run_model
import sys

model_type = str(sys.argv[1])
# option = str(sys.argv[2])
# plotname = 'Easy1'+option
plotname = 'Easy1'
print(model_type)

count, geneid, cellid = get_matrix_from_dir('pbmc8k')
geneid = geneid[:, 1]
count = count.T.tocsr()
seurat = np.genfromtxt('../pbmc8k/pbmc8k.seurat.labels', dtype='str', delimiter=',')
cellid = np.asarray([x.split('-')[0] for x in cellid])
labels_map = [0, 2, 4, 4, 0, 3, 3, 1, 5, 6]
cell_type = ["CD4+ T Helper2", "CD56+ NK", "CD14+ Monocyte", "CD19+ B", "CD8+ Cytotoxic T", "FCGR3A Monocyte",
             "dendritic"]
dataset1 = assign_label(cellid, geneid, labels_map, count, cell_type, seurat)
rmCellTypes = {'na', 'dendritic'}
newCellType = [k for i, k in enumerate(dataset1.cell_types) if k not in rmCellTypes]
dataset1.filter_cell_types(newCellType)

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

dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)

from scvi.models.vae import VAE
from scvi.inference import UnsupervisedTrainer
gene_dataset.subsample_genes(5000)
vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches, n_labels=gene_dataset.n_labels,
          n_hidden=128, n_latent=10, n_layers=2, dispersion='gene')
trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=1.0)
trainer.train(n_epochs=5)
full = trainer.create_posterior(vae, gene_dataset, indices=np.arange(len(gene_dataset)))
scale1 = full.sequential().get_harmonized_scale(0)
scale2 = full.sequential().get_harmonized_scale(1)
from scvi.inference.posterior import get_bayes_factors
bayes1 = get_bayes_factors(scale1,full.gene_dataset.labels.ravel(),0,4)
bayes2 = get_bayes_factors(scale2,full.gene_dataset.labels.ravel(),0,4)
def get_harmonized_scale(posterior, fixed_batch):
    px_scales = []
    for tensors in posterior:
        sample_batch, local_l_mean, local_l_var, batch_index, label = tensors
        px_scales += [posterior.model.scale_from_z(sample_batch, fixed_batch).cpu()]
    return np.concatenate(px_scales)


# if option == 'large':
#     from scvi.dataset.dataset10X import Dataset10X
#     dataset3 = Dataset10X('fresh_68k_pbmc_donor_a')
#     dataset3.cell_types = np.asarray(['unlabelled'])
#     dataset3.subsample_genes(dataset3.nb_genes)
#     dataset3.gene_names = dataset3.gene_symbols
#     gene_dataset = GeneExpressionDataset.concat_datasets(gene_dataset, dataset3)
#     latent, batch_indices, labels, keys = run_model(model_type, gene_dataset, dataset1, dataset2, ngenes=5000,filename=plotname)
# elif option == 'small':
latent, batch_indices, labels, keys = run_model(model_type, gene_dataset, dataset1, dataset2, ngenes=500,filename=plotname)
# if option=='large':
#     latent = latent[batch_indices!=2]
#     labels = labels[batch_indices!=2]
#     keys = gene_dataset.cell_types[np.unique(labels)]
#     _, labels = np.unique(labels,return_inverse=True)
#     batch_indices = batch_indices[batch_indices!=2]
#     batch_indices = batch_indices.reshape(len(batch_indices), 1)

eval_latent(batch_indices, labels, latent, keys, plotname + '.' + model_type)

labels = gene_dataset.labels

