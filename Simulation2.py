from scvi.harmonization.utils_chenling import CompareModels

use_cuda = True
import numpy as np
from scipy.sparse import csr_matrix

from scvi.dataset.dataset import GeneExpressionDataset


import sys
models = str(sys.argv[1])
plotname = 'Sim2'

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

UMI.subsample_genes(UMI.nb_genes)
nonUMI.subsample_genes(nonUMI.nb_genes)

gene_dataset = GeneExpressionDataset.concat_datasets(UMI, nonUMI)

CompareModels(gene_dataset, UMI, nonUMI, plotname, models)
