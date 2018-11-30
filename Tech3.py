from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import CompareModels
from scvi.dataset.MouseBrain import DentateGyrus10X,DentateGyrusC1
plotname = 'Tech3'
import sys
models = str(sys.argv[1])

dataset1= DentateGyrus10X()
dataset1.subsample_genes(dataset1.nb_genes)
dataset2 = DentateGyrusC1()
dataset2.subsample_genes(dataset2.nb_genes)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1,dataset2)

import numpy as np
f = open("../%s/celltypeprop.txt"%plotname, "w+")
f.write("%s\t"*len(gene_dataset.cell_types)%tuple(gene_dataset.cell_types)+"\n")
freq = [np.mean(gene_dataset.labels.ravel()==i) for i in np.unique(gene_dataset.labels.ravel())]
f.write("%f\t"*len(gene_dataset.cell_types)%tuple(freq)+"\n")
freq1 = [np.mean(gene_dataset.labels.ravel()[gene_dataset.batch_indices.ravel()==0]==i) for i in np.unique(gene_dataset.labels.ravel())]
f.write("%f\t"*len(gene_dataset.cell_types)%tuple(freq1)+"\n")
freq2 = [np.mean(gene_dataset.labels.ravel()[gene_dataset.batch_indices.ravel()==1]==i) for i in np.unique(gene_dataset.labels.ravel())]
f.write("%f\t"*len(gene_dataset.cell_types)%tuple(freq2)+"\n")
f.close()

# from scipy.sparse import csr_matrix
# gene_dataset.X = csr_matrix(gene_dataset.X)
CompareModels(gene_dataset, dataset1, dataset2, plotname, models)
