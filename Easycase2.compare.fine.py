use_cuda = True
from scvi.dataset.BICCN import *
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import CompareModels
from scvi.inference import UnsupervisedTrainer, SemiSupervisedTrainer

import sys
plotname = 'Macosko_Regev_fine'

dataset1 = MacoskoDataset(coarse=False)
dataset2 = RegevDataset(coarse=False)
# dataset1.subsample_genes(20000)
# dataset2.subsample_genes(20000)
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)
from scvi.harmonization.utils_chenling import SubsetGenes
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)
dataset1, dataset2, gene_dataset = SubsetGenes(dataset1, dataset2, gene_dataset, 'Macosko_Regev')

groups = ['Pvalb', 'L2/3', 'Sst', 'L5 PT', 'L5 IT Tcap', 'L5 IT Aldh1a7', 'L5 IT Foxp2', 'L5 NP',
                      'L6 IT', 'L6 CT', 'L6 NP', 'L6b', 'Lamp5', 'Vip', 'Astro', 'OPC', 'VLMC', 'Oligo', 'Sncg', 'Endo',
                      'SMC', 'MICRO']
groups = np.asarray([x.upper() for x in groups])
cell_type_bygroup = np.concatenate([[x for x in gene_dataset.cell_types if x.startswith(y)] for y in groups])
new_labels_dict = dict(zip(cell_type_bygroup, np.arange(len(cell_type_bygroup))))
labels = np.asarray([gene_dataset.cell_types[x] for x in gene_dataset.labels.ravel()])
new_labels = np.asarray([new_labels_dict[x] for x in labels])
labels_groups = [[i for i, x in enumerate(groups) if y.startswith(x)][0] for y in cell_type_bygroup]
dict_merged = dict(zip(cell_type_bygroup, groups[labels_groups]))
cell_groups = [dict_merged[x] for x in gene_dataset.cell_types]
gene_dataset.groups, gene_dataset.labels_groups = np.unique(cell_groups,return_inverse=True)
for i in np.unique(labels_groups):
    print(cell_type_bygroup[labels_groups==i])

freq = np.unique(new_labels,return_counts=True)[1]
for i in np.unique(labels_groups):
    print(freq[labels_groups==i])

gene_dataset.labels = new_labels
gene_dataset.cell_types = cell_type_bygroup

from scvi.harmonization.utils_chenling import trainVAE, SCANVI,SemiSupervisedTrainer
vae_posterior = trainVAE(gene_dataset,filename='Macosko_Regev',rep='0')
##########################################################################
# hierarchical
##########################################################################
scanvi = SCANVI(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels, n_layers=2,
                labels_groups=gene_dataset.labels_groups,use_labels_groups=True)
scanvi.load_state_dict(vae_posterior.model.state_dict(), strict=False)
trainer_scanvi = SemiSupervisedTrainer(scanvi, gene_dataset, classification_ratio=50,
                                       n_epochs_classifier=1, lr_classification=5 * 1e-3)
trainer_scanvi.train(n_epochs=50)
import torch
torch.save(trainer_scanvi.model, '../Macosko_Regev/scanvi.hier.pkl')
full = trainer_scanvi.create_posterior(trainer_scanvi.model, gene_dataset, indices=np.arange(len(gene_dataset)))
all_y, all_y_pred = full.compute_predictions()
acc = np.mean(all_y == all_y_pred)
all_y_groups = np.array([gene_dataset.labels_groups[y] for y in all_y])
all_y_pred_groups = np.array([gene_dataset.labels_groups[y] for y in all_y_pred])
h_acc = np.mean(all_y_groups == all_y_pred_groups)

# acc: 0.5097011494252873
# h_acc: 0.5340383141762453
##########################################################################
# not hierarchical
##########################################################################
scanvi = SCANVI(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels, n_layers=2,use_labels_groups=False)
scanvi.load_state_dict(vae_posterior.model.state_dict(), strict=False)
trainer_scanvi = SemiSupervisedTrainer(scanvi, gene_dataset, classification_ratio=50,
                                       n_epochs_classifier=1, lr_classification=5 * 1e-3)
trainer_scanvi.train(n_epochs=50)
full = trainer_scanvi.create_posterior(trainer_scanvi.model, gene_dataset, indices=np.arange(len(gene_dataset)))
import torch
torch.save(trainer_scanvi.model, '../Macosko_Regev/scanvi.fine.pkl')
all_y, all_y_pred = full.compute_predictions()
acc = np.mean(all_y == all_y_pred)
all_y_groups = np.array([gene_dataset.labels_groups[y] for y in all_y])
all_y_pred_groups = np.array([gene_dataset.labels_groups[y] for y in all_y_pred])
h_acc = np.mean(all_y_groups == all_y_pred_groups)
# acc: 0.4495785440613027
# hacc: 0.4787432950191571
##########################################################################
# coarse
##########################################################################
dataset1 = MacoskoDataset()
dataset2 = RegevDataset()
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)
scanvi = SCANVI(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels, n_layers=2,use_labels_groups=False)
scanvi.load_state_dict(vae_posterior.model.state_dict(), strict=False)
trainer_scanvi = SemiSupervisedTrainer(scanvi, gene_dataset, classification_ratio=50,
                                       n_epochs_classifier=1, lr_classification=5 * 1e-3)
trainer_scanvi.train(n_epochs=50)
import torch
torch.save(trainer_scanvi.model, '../Macosko_Regev/scanvi.coarse.pkl')
full = trainer_scanvi.create_posterior(trainer_scanvi.model, gene_dataset, indices=np.arange(len(gene_dataset)))
all_y, all_y_pred = full.compute_predictions()
acc = np.mean(all_y == all_y_pred)
# acc: 0.5303601532567049

#
# if model_type=='scanvi1':
#     trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices == 0))
#     trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices == 1))
# elif model_type=='scanvi2':
#     trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices == 1))
#     trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices == 0))
# elif model_type=='scanvi0':
#     trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices <0))
#     trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices >= 0))

