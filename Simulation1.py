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
label1 = np.load('../sim_data/label1.npy')
label2 = np.load('../sim_data/label2.npy')

dataset1 = GeneExpressionDataset(
            *GeneExpressionDataset.get_attributes_from_matrix(
                csr_matrix(count1), labels=label1),
            gene_names=['gene'+str(i) for i in range(2000)], cell_types=['type'+str(i+1) for i in range(5)])

dataset2 = GeneExpressionDataset(
            *GeneExpressionDataset.get_attributes_from_matrix(
                csr_matrix(count2), labels=label2),
            gene_names=['gene'+str(i) for i in range(2000)], cell_types=['type'+str(i+1) for i in range(5)])

gene_dataset = GeneExpressionDataset.concat_datasets(dataset1,dataset2)
gene_dataset.subsample_genes(1000)
CompareModels(gene_dataset, dataset1, dataset2, plotname, models)
