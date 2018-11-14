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
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)
CompareModels(gene_dataset,dataset1,dataset1,plotname,'others')



dict1 = dict(zip(dataset1.cell_types, dataset1.groups[dataset1.labels_groups]))
dict2 = dict(zip(dataset2.cell_types, dataset2.groups[dataset2.labels_groups]))
dict_merged = {**dict1, **dict2}
cell_groups = [dict_merged[x] for x in gene_dataset.cell_types]
gene_dataset.groups, gene_dataset.groups_labels = np.unique(cell_groups,return_inverse=True)
for i in np.unique(gene_dataset.groups_labels):
    print(gene_dataset.cell_types[gene_dataset.groups_labels==i])

freq = np.unique(gene_dataset.labels,return_counts=True)[1]
for i in np.unique(gene_dataset.groups_labels):
    print(freq[gene_dataset.groups_labels==i])


from scvi.harmonization.utils_chenling import trainVAE, SCANVI,SemiSupervisedTrainer
vae_posterior = trainVAE(gene_dataset,filename='Zeng',rep='1')
##########################################################################
# hierarchical
##########################################################################
scanvi = SCANVI(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels, n_layers=2,
                labels_groups=gene_dataset.groups_labels,use_labels_groups=True)
scanvi.load_state_dict(vae_posterior.model.state_dict(), strict=False)
trainer_scanvi = SemiSupervisedTrainer(scanvi, gene_dataset, classification_ratio=50,
                                       n_epochs_classifier=1, lr_classification=5 * 1e-3)
trainer_scanvi.train(n_epochs=50)
import torch
torch.save(trainer_scanvi.model, '../Macosko_Regev/scanvi.hier.pkl')
full = trainer_scanvi.create_posterior(trainer_scanvi.model, gene_dataset, indices=np.arange(len(gene_dataset)))
acc = full.hierarchical_accuracy(verbose=True)
print(acc)

# acc: 0.40706513409961687
# h_acc: 0.569088122605364
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
all_y_groups = np.array([gene_dataset.groups_labels[y] for y in all_y])
all_y_pred_groups = np.array([gene_dataset.groups_labels[y] for y in all_y_pred])
h_acc = np.mean(all_y_groups == all_y_pred_groups)
# acc: 0.5251340996168582
# hacc: 0.6727816091954023
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
# acc 0.3812107279693487


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

