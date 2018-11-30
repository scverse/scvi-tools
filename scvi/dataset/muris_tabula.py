import pandas as pd
use_cuda = True

from scvi.dataset.dataset import GeneExpressionDataset
import numpy as np
plotname = 'Tech1'
from anndata import read_h5ad
from scipy.sparse import csr_matrix

class TabulaMuris(GeneExpressionDataset):
    def __init__(self, dataname, save_path='/data/muris_tabula/', tissue='Marrow'):
        self.save_path = save_path
        self.dataname = dataname
        self.tissue = tissue
        self.urls =  ['https://github.com/czbiohub/tabula-muris-vignettes/raw/master/data/TM_droplet_metadata.csv',
                'https://github.com/czbiohub/tabula-muris-vignettes/raw/master/data/TM_facs_metadata.csv',
                'https://s3.amazonaws.com/czbiohub-tabula-muris/TM_droplet_mat.h5ad',
                'https://s3.amazonaws.com/czbiohub-tabula-muris/TM_facs_mat.h5ad']
        self.download_names = ['TM_droplet_metadata.csv','TM_facs_metadata.csv','TM_droplet_mat.h5ad','TM_facs_mat.h5ad']
        self.download()
        count, labels, cell_type, gene_names = self.preprocess()
        count = csr_matrix(count.astype('int'))
        super(TabulaMuris, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                count, labels=labels),
            gene_names=np.char.upper(gene_names), cell_types=cell_type)
    def preprocess(self):
        gene_len = pd.read_csv("../muris_tabula/gene_len.txt", sep=" ", low_memory=False, header=None)
        gene_len = np.asarray(gene_len)
        if self.dataname =='facs':
            data = read_h5ad('../muris_tabula/TM_facs_mat.h5ad').T
            meta = pd.read_csv("../muris_tabula/TM_facs_metadata.csv", sep=",", low_memory=False)
            check_meta = [data.obs.index[i] == meta.cell[i] for i in range(len(data.obs.index))]
            print(np.sum(check_meta) == len(data.obs.index))
            data = data[meta.tissue == self.tissue]
            scale = np.mean(gene_len[:, 1])
            data_dict = dict(zip(np.asarray(data._var.index), data.X.T))
            count = []
            genenames = []
            for name, length in gene_len:
                try:
                    count.append(data_dict[name] / length * scale)
                    genenames.append(name)
                except KeyError:
                    continue
            count = np.vstack(count).astype('int64')
            nonzero_cells = np.asarray(np.sum(count.T, axis=1)).ravel() > 100
            count = count[:,nonzero_cells]
            genenames = np.asarray(genenames)
            labels = meta.cell_ontology_class[meta.tissue==self.tissue]
            cell_type, labels = np.unique(np.asarray(labels[nonzero_cells]).astype('str'), return_inverse=True)
        elif self.dataname == 'droplet':
            data = read_h5ad('../muris_tabula/TM_droplet_mat.h5ad').T
            meta = pd.read_csv("../muris_tabula/TM_droplet_metadata.csv", sep=",", low_memory=False)
            data = data[meta.tissue == self.tissue]
            data_dict = dict(zip(np.asarray(data._var.index), data.X.T))
            count = []
            genenames = []
            for name, length in gene_len:
                try:
                    count.append(data_dict[name])
                    genenames.append(name)
                except KeyError:
                    continue
            count = np.vstack(count).astype('int64')
            genenames = np.asarray(genenames)
            labels = meta.cell_ontology_class[meta.tissue==self.tissue]
            cell_type, labels = np.unique(np.asarray(labels).astype('str'), return_inverse=True)
        return(count.T, labels, cell_type, genenames)

