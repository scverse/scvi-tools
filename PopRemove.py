use_cuda = True
from scvi.harmonization.utils_chenling import get_matrix_from_dir
from scvi.dataset.pbmc import PbmcDataset
from scvi.harmonization.benchmark import assign_label
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import eval_latent, run_model
from copy import deepcopy

import sys
model_type = str(sys.argv[1])
print(model_type)
plotname = 'PopRemove'

dataset1 = PbmcDataset()
dataset1.update_cells(dataset1.batch_indices.ravel()==0)

count, geneid, cellid = get_matrix_from_dir('cite')
count = count.T.tocsr()
seurat = np.genfromtxt('../cite/cite.seurat.labels', dtype='str', delimiter=',')
cellid = np.asarray([x.split('-')[0] for x in cellid])
labels_map = [0, 0, 1, 2, 3, 4, 5, 6]
labels = seurat[1:, 4]
cell_type = ['CD4 T cells', 'NK cells', 'CD14+ Monocytes', 'B cells','CD8 T cells', 'FCGR3A+ Monocytes', 'Other']
dataset2 = assign_label(cellid, geneid, labels_map, count, cell_type, seurat)
set(dataset2.cell_types).intersection(set(dataset2.cell_types))

dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)

for rmCellTypes in dataset2.cell_types:
    pbmc = deepcopy(dataset1)
    newCellType = [k for i, k in enumerate(dataset1.cell_types) if k not in [rmCellTypes]]
    pbmc.filter_cell_types(newCellType)
    gene_dataset = GeneExpressionDataset.concat_datasets(pbmc, dataset2)
    latent, batch_indices, labels, keys = run_model(model_type, gene_dataset, pbmc, dataset2, ngenes=500,filename=plotname+rmCellTypes.replace(' ',''))
    eval_latent(batch_indices, labels, latent, keys, plotname+rmCellTypes.replace(' ','') + '.' + model_type)
