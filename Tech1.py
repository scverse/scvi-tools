use_cuda = True
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import CompareModels
import sys
import numpy as np

models = str(sys.argv[1])
plotname = 'Tech1'

from scvi.dataset.muris_tabula import TabulaMuris
dataset1 = TabulaMuris('facs')
dataset2 = TabulaMuris('droplet')
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)

gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)
# genes = np.genfromtxt('../Seurat_data/'+plotname+'.CCA.genes.txt')
# genes = genes.astype('int')
# gene_dataset.X = gene_dataset.X[:,genes]
# gene_dataset.update_genes(genes)
# cells = np.genfromtxt('../Seurat_data/'+plotname+'.CCA.cells.txt')

gene_dataset.subsample_genes(500)

CompareModels(gene_dataset, dataset1, dataset2, plotname, models)
# 438.05354394990854

