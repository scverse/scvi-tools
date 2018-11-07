use_cuda = True
from scvi.dataset.BICCN import *
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import trainVAE

import sys
models = str(sys.argv[1])
plotname = 'four'


dataset1 = Zeng10X()
dataset2 = ZengSS2()
dataset3 = Zeng10X(cell_compartment='cell')
dataset4 = ZengSS2(cell_compartment='cell')

dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)
dataset3.subsample_genes(dataset3.nb_genes)
dataset4.subsample_genes(dataset4.nb_genes)

gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2,dataset3,dataset4)
gene_dataset.subsample_genes(10000)
full = trainVAE(gene_dataset,plotname, 0)
