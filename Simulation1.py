from scvi.harmonization.utils_chenling import CompareModels
use_cuda = True
import numpy as np
from scipy.sparse import csr_matrix
from scvi.dataset.dataset import GeneExpressionDataset

import sys
models = str(sys.argv[1])
plotname = 'Sim1'

count1 = np.load('../sim_data/count1.npy')
count2 = np.load('../sim_data/count2.npy')
label1 = np.load('../sim_data/label1.npy').astype('int')
label2 = np.load('../sim_data/label2.npy').astype('int')

dataset1 = GeneExpressionDataset(
            *GeneExpressionDataset.get_attributes_from_matrix(
                csr_matrix(count1), labels=label1),
            gene_names=['gene_'+str(i) for i in range(2000)], cell_types=['type'+str(i+1) for i in range(5)])

dataset2 = GeneExpressionDataset(
            *GeneExpressionDataset.get_attributes_from_matrix(
                csr_matrix(count2), labels=label2),
            gene_names=['gene_'+str(i) for i in range(2000)], cell_types=['type'+str(i+1) for i in range(5)])

dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)

gene_dataset = GeneExpressionDataset.concat_datasets(dataset1,dataset2)

# CompareModels(gene_dataset, dataset1, dataset2, plotname, 'writedata')

CompareModels(gene_dataset, dataset1, dataset2, plotname, models)

