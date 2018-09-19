from scvi.harmonization.utils_chenling import run_model, eval_latent
use_cuda = True
from scvi.dataset.BICCN import *
from scvi.dataset.dataset import GeneExpressionDataset

import sys
model_type = str(sys.argv[1])
plotname = 'Macosko_Regev'
print(model_type)

if model_type!='readSeurat':
    dataset1 = MacoskoDataset()
    dataset2 = RegevDataset()
    dataset1.subsample_genes(dataset1.nb_genes)
    dataset2.subsample_genes(dataset2.nb_genes)
    gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)
    gene_dataset.subsample_genes(5000)

    if model_type=='writedata':
        latent, batch_indices, labels, keys = run_model(model_type, gene_dataset, dataset1, dataset2, plotname)
    else:
        latent, batch_indices, labels, keys = run_model(model_type, gene_dataset, dataset1, dataset2)
elif model_type=='readSeurat':
    latent, batch_indices, labels, keys = run_model(model_type, 0, 0, 0, plotname)

if model_type.startswith('scanvi'):
    eval_latent(batch_indices, labels, latent, keys, plotname + '.' + model_type)
else:
    eval_latent(batch_indices, labels, latent, keys, plotname+'.'+model_type)
