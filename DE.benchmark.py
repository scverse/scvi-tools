from scvi.dataset.dataset10X import Dataset10X
from scvi.dataset.pbmc import PbmcDataset

use_cuda = True
from scvi.harmonization.utils_chenling import get_matrix_from_dir,assign_label
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import CompareModels
import sys
models = str(sys.argv[1])
plotname = 'DE'

pbmc = PbmcDataset()
pbmc.update_cells(pbmc.batch_indices.ravel()==0)
# pbmc.labels = pbmc.labels.reshape(len(pbmc),1)

donor = Dataset10X('fresh_68k_pbmc_donor_a')
donor.gene_names = donor.gene_symbols
donor.labels = np.repeat(0,len(donor)).reshape(len(donor),1)
donor.cell_types = ['unlabelled']

pbmc.subsample_genes(pbmc.nb_genes)
donor.subsample_genes(donor.nb_genes)
all_dataset = GeneExpressionDataset.concat_datasets(pbmc, donor)

CompareModels(all_dataset, pbmc, donor, plotname, models)

from scvi.dataset.dataset import SubsetGenes
dataset1, dataset2, gene_dataset = SubsetGenes(pbmc, donor, all_dataset, "DE")
