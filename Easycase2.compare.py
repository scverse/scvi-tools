use_cuda = True
from scvi.dataset.BICCN import *
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import CompareModels

import sys
models = str(sys.argv[1])
plotname = 'Macosko_Regev'
print(models)

dataset1 = MacoskoDataset()
dataset2 = RegevDataset()
# dataset1.subsample_genes(20000)
# dataset2.subsample_genes(20000)
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset1.nb_genes)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)

#17846 genes

CompareModels(gene_dataset, dataset1, dataset2, plotname, models)
