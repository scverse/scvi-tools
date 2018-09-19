use_cuda = True
from scvi.harmonization.utils_chenling import eval_latent, run_model

import numpy as np
from scipy.sparse import csr_matrix

from scvi.dataset.dataset import GeneExpressionDataset
import sys

model_type = str(sys.argv[1])
plotname = 'Sim3'
print(model_type)

count = np.load('../sim_data/Sim_EVFbatch.UMI.npy')
count = count.T
meta = np.load('../sim_data/Sim_EVFbatch.meta.npy')

dataset1 = GeneExpressionDataset(
    *GeneExpressionDataset.get_attributes_from_matrix(
        csr_matrix(count[meta[:, 2] == 0, :]), labels=meta[meta[:, 2] == 0, 1]),
    gene_names=['gene' + str(i) for i in range(2000)], cell_types=['type' + str(i + 1) for i in range(5)])

dataset2 = GeneExpressionDataset(
    *GeneExpressionDataset.get_attributes_from_matrix(
        csr_matrix(count[meta[:, 2] == 1, :]), labels=meta[meta[:, 2] == 1, 1]),
    gene_names=['gene' + str(i) for i in range(2000)], cell_types=['type' + str(i + 1) for i in range(5)])

gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)

latent, batch_indices, labels,keys = run_model(model_type, gene_dataset, dataset1, dataset2, filename=plotname, ngenes=5000)

if model_type.startswith('scanvi'):
    eval_latent(batch_indices, labels, latent, keys, plotname + '.' + model_type)
else:
    eval_latent(batch_indices, labels, latent, keys, plotname+'.'+model_type)
