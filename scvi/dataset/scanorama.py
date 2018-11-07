from scvi.dataset import GeneExpressionDataset
import scipy.sparse as sparse
import numpy as np
import os
import sys

class DatasetSCANORAMA(GeneExpressionDataset):

    def __init__(self, filename, save_path='/data/scanorama/'):
        self.save_path = save_path + '%s' % filename
        count, gene_names = self.preprocess()
        super(DatasetSCANORAMA, self).__init__(*GeneExpressionDataset.get_attributes_from_matrix(
                count),
            gene_names=np.char.upper(gene_names))

    def preprocess(self):
        print("Preprocessing dataset: " + self.save_path)
        if os.path.isfile(self.save_path + '.h5.npz'):
            X = sparse.load_npz(self.save_path + '.h5.npz')
            with open(self.save_path + '.h5.genes.txt') as f:
                genes = np.array(f.read().rstrip().split())
        elif os.path.isfile(self.save_path + '.npz'):
            data = np.load(self.save_path + '.npz')
            X = data['X']
            genes = data['genes']
            data.close()
        elif os.path.isfile(self.save_path + '/tab.npz'):
            X = sparse.load_npz(self.save_path + '/tab.npz')
            with open(self.save_path + '/tab.genes.txt') as f:
                genes = np.array(f.read().rstrip().split())
        else:
            sys.stderr.write('Could not find: {}\n'.format(self.save_path))
            exit(1)
        genes = np.array([gene.upper() for gene in genes])
        return X, genes


