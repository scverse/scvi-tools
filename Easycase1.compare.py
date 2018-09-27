use_cuda = True
from scvi.harmonization.utils_chenling import get_matrix_from_dir,assign_label
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import CompareModels
import sys

models = str(sys.argv[1])
plotname = 'Easy1'

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

genes = np.genfromtxt('../Seurat_data/'+plotname+'.CCA.genes.txt')
genes = genes.astype('int')
gene_dataset.X = gene_dataset.X[:,genes]
gene_dataset.update_genes(genes)

CompareModels(gene_dataset, dataset1, dataset2, plotname, models)


from copy import deepcopy
from scvi.harmonization.utils_chenling import run_model
from scvi.harmonization.utils_chenling import eval_latent

dataset1 = deepcopy(gene_dataset)
dataset1.update_cells(gene_dataset.batch_indices.ravel() == 0)
dataset1.subsample_genes(dataset1.nb_genes)
latent1, _, _, _, _ = run_model('vae', dataset1, 0, 0, filename=plotname)
dataset2 = deepcopy(gene_dataset)
dataset2.update_cells(gene_dataset.batch_indices.ravel()  == 1)
dataset2.subsample_genes(dataset2.nb_genes)
latent2, _, _, _, _ = run_model('vae', dataset2, 0, 0, filename=plotname)
for model_type in ['vae','scanvi','scanvi1','scanvi2','scanvi0']:
print(model_type)
latent, batch_indices, labels, keys,stats = run_model(model_type, gene_dataset, dataset1, dataset2, filename=plotname)
unlabelled = stats[4]
latent1_id = unlabelled[unlabelled<latent1.shape[0]]
latent2_id = unlabelled[unlabelled>=latent1.shape[0]]-latent1.shape[0]
res_knn, res_kmeans, res_jaccard = eval_latent(batch_indices, labels, latent, latent1[latent1_id,:], latent2[latent2_id,:], keys, plotname + '.' + model_type, plotting=True)
for i in [1, 2, 3]:
    latent, batch_indices, labels, keys, stats = run_model(model_type, gene_dataset, dataset1, dataset2,filename=plotname)
    res_knn, res_kmeans, res_jaccard = eval_latent(batch_indices, labels, latent, latent1, latent2,
                                                   keys, plotname + '.' + model_type, plotting=True)
    res = [res_knn[x] for x in res_knn] + [res_kmeans[x] for x in res_kmeans] + [res_jaccard]
    f.write(model_type+str(i) + " %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f" % tuple(res))
