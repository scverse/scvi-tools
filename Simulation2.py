from scvi.harmonization.utils_chenling import eval_latent, run_model

use_cuda = True
import numpy as np
from scipy.sparse import csr_matrix

from scvi.dataset.dataset import GeneExpressionDataset


import sys
model_type = str(sys.argv[1])
plotname = 'Sim2'
print(model_type)

countUMI = np.load('../sim_data/count.UMI.npy')
countnonUMI = np.load('../sim_data/count.nonUMI.npy')
labelUMI = np.load('../sim_data/label.UMI.npy')
labelnonUMI = np.load('../sim_data/label.nonUMI.npy')

UMI = GeneExpressionDataset(
            *GeneExpressionDataset.get_attributes_from_matrix(
                csr_matrix(countUMI.T), labels=labelUMI),
            gene_names=['gene'+str(i) for i in range(2000)], cell_types=['type'+str(i+1) for i in range(5)])

nonUMI = GeneExpressionDataset(
            *GeneExpressionDataset.get_attributes_from_matrix(
                csr_matrix(countnonUMI.T), labels=labelnonUMI),
            gene_names=['gene'+str(i) for i in range(2000)], cell_types=['type'+str(i+1) for i in range(5)])

gene_dataset = GeneExpressionDataset.concat_datasets(UMI, nonUMI)

latent, batch_indices, labels,keys = run_model(model_type, gene_dataset, UMI, nonUMI,filename=plotname, ngenes=5000)

if model_type.startswith('scanvi'):
    eval_latent(batch_indices, labels, latent, keys, plotname + '.' + model_type)
else:
    eval_latent(batch_indices, labels, latent, keys, plotname+'.'+model_type)
