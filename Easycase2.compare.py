from scvi.harmonization.utils_chenling import run_model, eval_latent
use_cuda = True
from scvi.dataset.BICCN import *
from scvi.dataset.dataset import GeneExpressionDataset

import sys
model_type = str(sys.argv[1])
plotname = 'Macosko_Regev'
print(model_type)

dataset1 = MacoskoDataset()
dataset2 = RegevDataset()
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)

genes = np.genfromtxt('../Seurat_data/'+plotname+'.CCA.genes.txt')
genes = genes.astype('int')
gene_dataset.X = gene_dataset.X[:,genes]
gene_dataset.update_genes(genes)

cells = np.genfromtxt('../Seurat_data/'+plotname+'.CCA.cells.txt')
print(cells.shape)
print(gene_dataset.X.shape)

latent, batch_indices, labels, keys = run_model(model_type, gene_dataset, dataset1, dataset2, plotname)
eval_latent(batch_indices, labels, latent, keys, plotname + '.' + model_type,plotting=True)

for i in [1,2,3]:
    latent, batch_indices, labels, keys = run_model(model_type, gene_dataset, dataset1, dataset2, plotname)
    eval_latent(batch_indices, labels, latent, keys, plotname + '.' + model_type, plotting=False)
