from scvi.harmonization.utils_chenling import get_matrix_from_dir,assign_label,KNNJaccardIndex
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import eval_latent, run_model
from copy import deepcopy
import sys
use_cuda = True

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


# cells = np.genfromtxt('../Seurat_data/'+plotname+'.CCA.cells.txt')
latent, batch_indices, labels, keys = run_model(model_type, gene_dataset, dataset1, dataset2,filename=plotname)
eval_latent(batch_indices, labels, latent, keys, plotname + '.' + model_type, plotting=True)

dataset1 = deepcopy(gene_dataset)
dataset1.update_cells(batch_indices==0)
latent1, _, _, _ = run_model(model_type, dataset1, 0, 0,filename=plotname)

dataset2 = deepcopy(gene_dataset)
dataset2.update_cells(batch_indices==1)
latent2, _, _, _ = run_model(model_type, dataset2, 0, 0,filename=plotname)

KNeighbors  = np.concatenate([np.arange(10,100,10),np.arange(100,500,50)])
res_vae = [KNNJaccardIndex(latent1,latent2,latent,batch_indices,i)[0] for i in KNeighbors]

for i in [1,2,3]:
    latent, batch_indices, labels, keys = run_model(model_type, gene_dataset, dataset1, dataset2,filename=plotname)
    eval_latent(batch_indices, labels, latent, keys, plotname + '.' + model_type,plotting=False)
