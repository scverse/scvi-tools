from scvi.dataset.dataset import GeneExpressionDataset
from sklearn.datasets import dump_svmlight_file, load_svmlight_file
import numpy as np
import loompy
import os


class ZeiselMoleArchData(GeneExpressionDataset):
    def __init__(self, save_path='/data/mouse_brain/',coarse=True):
        self.save_path = save_path
        count, labels, cell_type, gene_names,labels_groups,groups = self.preprocess()
        labels = labels.astype('int')
        if coarse==True:
            labels = labels_groups[labels]
            cell_type = groups
        super(ZeiselMoleArchData, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                count, labels=labels),
            gene_names=np.char.upper(gene_names), cell_types=cell_type)
        if coarse==False:
            self.labels_groups = labels_groups
            self.groups = groups
    def preprocess(self):
        if os.path.isfile(self.save_path + 'zeisel.svmlight'):
            count, labels = load_svmlight_file(self.save_path + 'zeisel.svmlight')
            cell_type = np.load(self.save_path + 'zeisel.celltypes.npy')
            labels_groups = np.load(self.save_path + 'zeisel.labels_groups.npy')
            gene_names = np.load(self.save_path + 'zeisel.gene_names.npy')
            groups = np.load(self.save_path + 'zeisel.groups.npy')
            return(count, labels, cell_type, gene_names, labels_groups,groups)
        else:
            cells = loompy.connect(self.save_path + 'l5_all.loom')
            select = cells[:, :].sum(axis=0) > 0  # Take out cells that doesn't express any gene
            count = cells[:, select].T  # change matrix to cells by genes
            cluster = loompy.connect(self.save_path + 'l5_all.agg.loom')
            fine = cells.ca['TaxonomyRank4']
            cell_types, labels = np.unique(fine, return_inverse=True)
            assert np.sum(cell_types[labels] == fine) == len(fine)
            coarse1 = cluster.ca['TaxonomyRank3']
            coarse2 = cluster.ca['TaxonomyRank2']
            coarse1[coarse2 == 'CNS neurons'] = 'CNS neurons'
            coarse = coarse1
            mapping = dict(zip(cluster.ca['TaxonomyRank4'], coarse))
            group_names = [mapping[x] for x in cell_types]
            groups, labels_groups = np.unique(group_names, return_inverse=True)
            coarse1 = cells.ca['TaxonomyRank3']
            coarse2 = cells.ca['TaxonomyRank2']
            coarse1[coarse2 == 'CNS neurons'] = 'CNS neurons'
            assert np.sum(groups[labels_groups[labels]] == coarse1) == len(coarse1)
            dump_svmlight_file(count, labels, self.save_path + 'zeisel.svmlight')
            np.save(self.save_path + 'zeisel.celltypes.npy', cell_types)
            np.save(self.save_path + 'zeisel.gene_names.npy', cells.ra['Gene'])
            np.save(self.save_path + 'zeisel.labels_groups.npy', labels_groups)
            np.save(self.save_path + 'zeisel.groups.npy', groups)
            return (count, labels, cell_types, cells.ra['Gene'], labels_groups,groups)




class ZeiselCortexOnly(GeneExpressionDataset):
    def __init__(self, save_path='/data/mouse_brain/cortex1/',coarse=True):
        self.save_path = save_path
        count, labels, cell_type, gene_names,labels_groups,groups,batch = self.preprocess()
        labels = labels.astype('int')
        batch_names,batch_indices = np.unique(batch,return_inverse=True)
        if coarse==True:
            labels = labels_groups[labels]
            cell_type = groups
        super(ZeiselCortexOnly, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                count, labels=labels),
            gene_names=np.char.upper(gene_names), cell_types=cell_type)
        if coarse==False:
            self.labels_groups = labels_groups
            self.groups = groups
        self.batch_indices = batch_indices.reshape(len(batch_indices),1)
        self.batch_names = batch_names
    def preprocess(self):
        if os.path.isfile(self.save_path + 'zeiselCortex.svmlight'):
            count, labels = load_svmlight_file(self.save_path + 'zeiselCortex.svmlight')
            cell_type = np.load(self.save_path + 'zeiselCortex.celltypes.npy')
            labels_groups = np.load(self.save_path + 'zeiselCortex.labels_groups.npy')
            gene_names = np.load(self.save_path + 'zeiselCortex.gene_names.npy')
            groups = np.load(self.save_path + 'zeiselCortex.groups.npy')
            batch = np.load(self.save_path + 'zeiselCortex.batch.npy')
            return(count, labels, cell_type, gene_names, labels_groups,groups,batch)
        else:
            cells = loompy.connect(self.save_path + 'l1_cortex1.loom')
            select = cells[:, :].sum(axis=0) > 0  # Take out cells that doesn't express any gene
            count = cells[:, select].T  # change matrix to cells by genes
            batch = cells.ca['DonorID']
            fine = cells.ca['Subclass']
            cell_types, labels = np.unique(fine, return_inverse=True)
            assert np.sum(cell_types[labels] == fine) == len(fine)
            mapping = dict(zip(cells.ca['Subclass'], cells.ca['Class']))
            group_names = [mapping[x] for x in cell_types]
            groups, labels_groups = np.unique(group_names, return_inverse=True)
            dump_svmlight_file(count, labels, self.save_path + 'zeiselCortex.svmlight')
            np.save(self.save_path + 'zeiselCortex.celltypes.npy', cell_types)
            np.save(self.save_path + 'zeiselCortex.gene_names.npy', cells.ra['Gene'])
            np.save(self.save_path + 'zeiselCortex.labels_groups.npy', labels_groups)
            np.save(self.save_path + 'zeiselCortex.groups.npy', groups)
            np.save(self.save_path + 'zeiselCortex.batch.npy', batch)
            return (count, labels, cell_types, cells.ra['Gene'], labels_groups,groups,batch)


class DentateGyrus10X(GeneExpressionDataset):
    def __init__(self, save_path='/data/mouse_brain/'):
        self.save_path = save_path
        count, labels, cell_type, gene_names= self.preprocess()
        labels = labels.astype('int')
        super(DentateGyrus10X, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                count, labels=labels),
            gene_names=np.char.upper(gene_names), cell_types=cell_type)
    def preprocess(self):
        cells = loompy.connect('/data/mouse_brain/' + 'dentate_gyrus_10X_V1.loom')
        select = cells[:, :].sum(axis=0) > 0  # Take out cells that doesn't express any gene
        count = cells[:, select].T  # change matrix to cells by genes
        cell_types, labels = np.unique(cells.ca['Cluster'], return_inverse=True)
        return (count, labels, cell_types, cells.ra['Gene'])


class DentateGyrusC1(GeneExpressionDataset):
    def __init__(self, save_path='/data/mouse_brain/'):
        self.save_path = save_path
        count, labels, cell_type, gene_names= self.preprocess()
        count = np.floor(count)
        labels = labels.astype('int')
        super(DentateGyrusC1, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                count, labels=labels),
            gene_names=np.char.upper(gene_names), cell_types=cell_type)
    def preprocess(self):
        cells = loompy.connect('/data/mouse_brain/' + 'dentate_gyrus_C1.loom')
        select = cells[:, :].sum(axis=0) > 0  # Take out cells that doesn't express any gene
        count = cells[:, select].T  # change matrix to cells by genes
        cell_types, labels = np.unique(cells.ca['Cluster'], return_inverse=True)
        return (count, labels, cell_types, cells.ra['Gene'])

