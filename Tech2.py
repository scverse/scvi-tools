use_cuda = True

from scvi.dataset.BICCN import *
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import CompareModels
import sys

models = str(sys.argv[1])
plotname = 'Zeng'

dataset1 = Zeng10X()
dataset2 = ZengSS2()
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)

gene_dataset.subsample_genes(2000)
#
# genes = np.genfromtxt('../Seurat_data/'+plotname+'.CCA.genes.txt')
# genes = genes.astype('int')
# gene_dataset.X = gene_dataset.X[:,genes]
# gene_dataset.update_genes(genes)
#
# cells = np.genfromtxt('../Seurat_data/'+plotname+'.CCA.cells.txt')
# print(cells.shape)
# print(gene_dataset.X.shape)

CompareModels(gene_dataset, dataset1, dataset2, plotname, models)
