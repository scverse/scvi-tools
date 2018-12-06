use_cuda = True
from scvi.dataset.MouseBrain import ZeiselMoleArchData
from scvi.harmonization.utils_chenling import CompareModels

import sys
models = str(sys.argv[1])

gene_dataset = ZeiselMoleArchData(coarse=False)
from copy import deepcopy
dataset1 = deepcopy(gene_dataset)
dataset1.update_cells(dataset1.batch_indices.ravel()==0)
dataset2 = deepcopy(gene_dataset)
dataset2.update_cells(dataset2.batch_indices.ravel()==1)


CompareModels(gene_dataset, dataset1, dataset2, 'Zeisel', models)
