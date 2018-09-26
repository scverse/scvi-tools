use_cuda = True
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import run_model, eval_latent
import sys
import numpy as np
from copy import deepcopy
model_type = str(sys.argv[1])
plotname = 'Tech1'
print(model_type)

from scvi.dataset.muris_tabula import TabulaMuris
dataset1 = TabulaMuris('facs')
dataset2 = TabulaMuris('droplet')
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

latent, batch_indices, labels, keys = run_model(model_type, gene_dataset, dataset1, dataset2, ngenes=500,
                                                filename=plotname)
eval_latent(batch_indices, labels, latent, keys, plotname + '.' + model_type, plotting=True)

for i in [1,2,3]:
    latent, batch_indices, labels, keys = run_model(model_type, gene_dataset, dataset1, dataset2, ngenes=500,filename=plotname)
    eval_latent(batch_indices, labels, latent, keys, plotname + '.' + model_type,plotting=False)

#
# latent, batch_indices, labels,keys = run_model(model_type, gene_dataset, dataset1, dataset2, filename=plotname, ngenes=5000)
# eval_latent(batch_indices, labels, latent, keys, plotname + '.' + model_type)
