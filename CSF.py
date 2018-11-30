use_cuda = True
from scvi.harmonization.utils_chenling import get_matrix_from_dir,assign_label
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import get_matrix_from_dir
from scvi.dataset.dataset import GeneExpressionDataset
from copy import deepcopy
from scvi.harmonization.utils_chenling import CompareModels
import sys
import pandas as pd

models = str(sys.argv[1])
plotname = 'CSF'

class MSDataset(GeneExpressionDataset):
    def __init__(self):
        count, labels, cell_type, gene_names,cellid = self.preprocess()
        super(MSDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                count.tocsr(), labels=labels),
            gene_names=np.char.upper(gene_names), cell_types=cell_type)
    def preprocess(self):
        print("Preprocessing data")
        count, geneid, cellid = get_matrix_from_dir('CSF',storage='/data/')
        print("Finish preprocessing data")
        labels = np.repeat(0,len(cellid))
        cell_type = ['unlabelled']
        return count,labels,cell_type,geneid,cellid



gene_dataset = MSDataset()
celltypes = pd.read_csv('/data/CSF/celltypes.txt')
celltype = np.asarray(celltypes['celltypes'])
donor = np.asarray(celltypes['donor'])

gene_dataset.cell_types,gene_dataset.labels = np.unique(celltype,return_inverse=True)
gene_dataset.batch_names,gene_dataset.batch_indices = np.unique(donor,return_inverse=True)
gene_dataset.batch_indices = gene_dataset.batch_indices.reshape(len(gene_dataset.batch_indices),1)

dataset1 = deepcopy(gene_dataset)
dataset1.update_cells(dataset1.batch_indices.ravel()==3)
dataset2 = deepcopy(gene_dataset)
dataset2.update_cells(dataset2.batch_indices.ravel()==7)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1,dataset2)


f = open("../%s/celltypeprop.txt"%plotname, "w+")
f.write("%s\t"*len(gene_dataset.cell_types)%tuple(gene_dataset.cell_types)+"\n")
freq = [np.mean(gene_dataset.labels.ravel()==i) for i in np.unique(gene_dataset.labels.ravel())]
f.write("%f\t"*len(gene_dataset.cell_types)%tuple(freq)+"\n")
freq1 = [np.mean(gene_dataset.labels.ravel()[gene_dataset.batch_indices.ravel()==0]==i) for i in np.unique(gene_dataset.labels.ravel())]
f.write("%f\t"*len(gene_dataset.cell_types)%tuple(freq1)+"\n")
freq2 = [np.mean(gene_dataset.labels.ravel()[gene_dataset.batch_indices.ravel()==1]==i) for i in np.unique(gene_dataset.labels.ravel())]
f.write("%f\t"*len(gene_dataset.cell_types)%tuple(freq2)+"\n")
f.close()

CompareModels(gene_dataset, dataset1, dataset2, plotname, models)

