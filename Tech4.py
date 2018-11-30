from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import CompareModels
from scipy.sparse import csr_matrix
plotname = 'Tech4'
import sys
from copy import deepcopy
models = str(sys.argv[1])
import numpy as np

from scvi.dataset.scanorama import DatasetSCANORAMA
dirs = ['data/pancreas/pancreas_inDrop', 'data/pancreas/pancreas_multi_celseq2_expression_matrix', 'data/pancreas/pancreas_multi_celseq_expression_matrix', 'data/pancreas/pancreas_multi_fluidigmc1_expression_matrix', 'data/pancreas/pancreas_multi_smartseq2_expression_matrix']

datasets = [DatasetSCANORAMA(d) for d in dirs]

labels = (open('/data/scanorama/data/cell_labels/pancreas_cluster.txt').read().rstrip().split())

all_dataset = GeneExpressionDataset.concat_datasets(*datasets)
batch_id = all_dataset.batch_indices.ravel()

all_dataset = GeneExpressionDataset.concat_datasets(datasets[0],datasets[1])
all_dataset.cell_types,all_dataset.labels = np.unique(np.asarray(labels)[(batch_id==0)+(batch_id==1)],return_inverse=True)
all_dataset.labels = all_dataset.labels.reshape(len(all_dataset.labels),1)
all_dataset.n_labels = len(np.unique(all_dataset.labels))

all_dataset.subsample_genes(all_dataset.nb_genes)
dataset1 = deepcopy(all_dataset)
dataset1.update_cells(dataset1.batch_indices.ravel()==0)
dataset2 = deepcopy(all_dataset)
dataset2.update_cells(dataset2.batch_indices.ravel()==1)
# all_dataset.X = csr_matrix(all_dataset.X)
# CompareModels(all_dataset, dataset1, dataset2, plotname, 'writedata')

f = open("../%s/celltypeprop.txt"%plotname, "w+")
f.write("%s\t"*len(all_dataset.cell_types)%tuple(all_dataset.cell_types)+"\n")
freq = [np.mean(all_dataset.labels.ravel()==i) for i in np.unique(all_dataset.labels.ravel())]
f.write("%f\t"*len(all_dataset.cell_types)%tuple(freq)+"\n")
freq1 = [np.mean(all_dataset.labels.ravel()[all_dataset.batch_indices.ravel()==0]==i) for i in np.unique(all_dataset.labels.ravel())]
f.write("%f\t"*len(all_dataset.cell_types)%tuple(freq1)+"\n")
freq2 = [np.mean(all_dataset.labels.ravel()[all_dataset.batch_indices.ravel()==1]==i) for i in np.unique(all_dataset.labels.ravel())]
f.write("%f\t"*len(all_dataset.cell_types)%tuple(freq2)+"\n")
f.close()
CompareModels(all_dataset, dataset1, dataset2, plotname, models)
