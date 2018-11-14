from scvi.harmonization.utils_chenling import get_matrix_from_h5, TryFindCells
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset
from sklearn.datasets import dump_svmlight_file, load_svmlight_file
import os
from scipy.sparse import csr_matrix

class MacoskoDataset(GeneExpressionDataset):
    def __init__(self, save_path='../AIBS/',coarse=True):
        self.save_path = save_path
        count, gene_names, labels, cell_type, label_groups, groups = self.preprocess()
        labels = labels.astype('int')
        print('finished preprocessing')
        if coarse==True:
            new_labels_dict = dict(zip(*[np.unique(labels), groups[label_groups]]))
            labels = np.asarray([new_labels_dict[x] for x in labels])
            groups, labels = np.unique(labels, return_inverse=True)
            cell_type = groups
            print('finished coarse labeling')
        super(MacoskoDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                count, labels=labels),
            gene_names=np.char.upper(gene_names), cell_types=cell_type)
        self.labels_groups = label_groups
        self.groups = groups
    def preprocess(self):
        groups = ['NP', 'OLIGO', 'SST', 'SST_CHODL', 'IT', 'IT_CAR', 'PT', 'PVALB', 'L6b', 'VIP', 'LAMP5', 'OPC',
                  'ASTRO', 'Endo']
        groups = np.asarray([x.upper() for x in groups])
        if os.path.isfile(self.save_path + 'macosko_data.svmlight'):
            count, labels = load_svmlight_file(self.save_path + 'macosko_data.svmlight')
            cell_type = np.load(self.save_path + 'macosko_data.celltypes.npy')
            label_groups = np.load(self.save_path + 'macosko_data.labels_groups.npy')
            gene_names = np.load(self.save_path + 'macosko_data.gene_names.npy')
            return(count, gene_names, labels, cell_type, label_groups, groups)
        else:
            macosko_batches = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
            label = np.genfromtxt(self.save_path + '10X_nuclei_Macosko/cluster.membership.csv', dtype='str',
                                  delimiter=',')
            label_batch = np.asarray([str(int(int(x.split('-')[1].split('"')[0]) / 11)) for x in label[1:, 0]])
            label_barcode = np.asarray([x.split('-')[0].split('"')[1] for x in label[1:, 0]])
            label_cluster = np.asarray([x.split('"')[1] for x in label[1:, 1]])
            label_map = np.genfromtxt(self.save_path + '10X_nuclei_Macosko/cluster.annotation.csv', dtype='str',
                                      delimiter=',')
            label_map = dict(
                zip([x.split('"')[1] for x in label_map[:, 0]], [x.split('"')[1] for x in label_map[:, 1]]))
            macosko_data = []
            for batch_i, batch in enumerate(macosko_batches):
                geneid, cellid, count = get_matrix_from_h5(
                    self.save_path + '10X_nuclei_Macosko/' + '171218_p56m1' + batch + '/outs/filtered_gene_bc_matrices_h5.h5',
                    'mm10_premrna')
                count = count.T.tocsr()
                print(count.shape, len(geneid), len(cellid))
                cellid = [id.split('-')[0] for id in cellid]
                label_dict = dict(
                    zip(label_barcode[label_batch == str(batch_i + 1)], label_cluster[label_batch == str(batch_i + 1)]))
                new_count, matched_label = TryFindCells(label_dict, cellid, count)
                cell_type = [label_map[x] for x in np.unique(matched_label)]
                matched_label = np.asarray(matched_label.astype('int')) - 1
                dataset = GeneExpressionDataset(
                    *GeneExpressionDataset.get_attributes_from_matrix(new_count, labels=matched_label),
                    gene_names=geneid, cell_types=cell_type)
                print(dataset.X.shape, len(dataset.labels))
                macosko_data.append(dataset)
            dataset = GeneExpressionDataset.concat_datasets(*macosko_data)
            cell_type = dataset.cell_types
            cell_type = np.asarray([x.upper() for x in cell_type])
            label_groups = [0,1,2,4,7,4,8,1,3,4,6,9,7,5,10,9,6,4,7,13,4,4,2,8,8,1,11,7,12,7,7,10,4,9,4,2,4,6]
            # for i in np.unique(label_groups):
            #     print(np.asarray(cell_type)[np.asarray(label_groups)==i])
            #
            dump_svmlight_file(dataset.X, dataset.labels.ravel(), self.save_path + 'macosko_data.svmlight')
            np.save(self.save_path + 'macosko_data.celltypes.npy', cell_type)
            np.save(self.save_path + 'macosko_data.gene_names.npy', dataset.gene_names)
            np.save(self.save_path + 'macosko_data.labels_groups.npy', label_groups)
            return (dataset.X, dataset.gene_names, dataset.labels.ravel(), cell_type, label_groups, groups)



class RegevDataset(GeneExpressionDataset):
    def __init__(self, save_path='../AIBS/',coarse=True):
        self.save_path = save_path
        count, gene_names, labels, cell_type, label_groups, groups = self.preprocess()
        labels = labels.astype('int')
        print('finished preprocessing')
        if coarse==True:
            new_labels_dict = dict(zip(*[np.unique(labels), groups[label_groups]]))
            labels = np.asarray([new_labels_dict[x] for x in labels])
            groups, labels = np.unique(labels, return_inverse=True)
            cell_type = groups
            print('finished coarse labeling')
        super(RegevDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                count, labels=labels),
            gene_names=np.char.upper(gene_names), cell_types=cell_type)
        self.labels_groups = label_groups
        self.groups = groups
    def preprocess(self):
        groups = ['L6b', 'IT', 'PVALB', 'IT_CAR', 'LAMP5', 'VLMC', 'PT', 'Endo', 'SST', 'VIP', 'Oligo', 'NP', 'OPC',
                  'Astro']
        groups = np.asarray([x.upper() for x in groups])
        if os.path.isfile(self.save_path + 'regev_data.svmlight'):
            count, labels = load_svmlight_file(self.save_path + 'regev_data.svmlight')
            cell_type = np.load(self.save_path + 'regev_data.celltypes.npy')
            gene_names = np.load(self.save_path + 'regev_data.gene_names.npy')
            label_groups = np.load(self.save_path + 'regev_data.labels_groups.npy')
            return(count, gene_names, labels, cell_type, label_groups, groups)
        else:
            regev_batches = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
            label = np.genfromtxt(self.save_path + '10X_nuclei_Regev/cluster.membership.csv', dtype='str',
                                  delimiter=',')
            label_batch = np.asarray([str(int(int(x.split('-')[1].split('"')[0]))) for x in label[1:, 0]])
            label_barcode = np.asarray([x.split('-')[0].split('"')[1] for x in label[1:, 0]])
            label_cluster = np.asarray([x.split('"')[1] for x in label[1:, 1]])
            label_map = np.genfromtxt(self.save_path + '10X_nuclei_Regev/cluster.annotation.csv', dtype='str',
                                      delimiter=',')
            label_map = dict(
                zip([x.split('"')[1] for x in label_map[:, 0]], [x.split('"')[1] for x in label_map[:, 1]]))
            regev_data = []
            for batch_i, batch in enumerate(regev_batches):
                geneid, cellid, count = get_matrix_from_h5(
                    self.save_path + '10X_nuclei_Regev/' + batch + '1/filtered_gene_bc_matrices_h5.h5',
                    'mm10-1.2.0_premrna')
                count = count.T.tocsr()
                cellid = [id.split('-')[0] for id in cellid]
                label_dict = dict(
                    zip(label_barcode[label_batch == str(batch_i + 1)], label_cluster[label_batch == str(batch_i + 1)]))
                new_count, matched_label = TryFindCells(label_dict, cellid, count)
                cell_type = [label_map[x] for x in np.unique(matched_label)]
                matched_label = np.asarray(matched_label.astype('int'))-1
                dataset = GeneExpressionDataset(
                    *GeneExpressionDataset.get_attributes_from_matrix(new_count, labels=matched_label),
                    gene_names=geneid, cell_types=cell_type)
                regev_data.append(dataset)
            dataset = GeneExpressionDataset.concat_datasets(*regev_data)
            cell_type = dataset.cell_types
            cell_type = [x.upper() for x in cell_type]
            label_groups = [0,1,2,1,1,1,0,1,1,1,1,0,3,1,1,4,5,6,7,1,0,2,7,0,8,4,9,1,10,11,12,13,1,1,7,2,13]
            # for i in np.unique(label_groups):
            #     print(np.asarray(cell_type)[np.asarray(label_groups)==i])
            #
            dump_svmlight_file(dataset.X, dataset.labels.ravel(), self.save_path+'regev_data.svmlight')
            np.save(self.save_path + 'regev_data.celltypes.npy', cell_type)
            np.save(self.save_path + 'regev_data.gene_names.npy', dataset.gene_names)
            np.save(self.save_path + 'regev_data.labels_groups.npy', label_groups)
            return (dataset.X, dataset.gene_names, dataset.labels.ravel(), cell_type, label_groups, groups)


class Zeng10X(GeneExpressionDataset):
    def __init__(self, save_path='../AIBS/', coarse=True):
        self.save_path = save_path
        count, gene_names, labels, cell_type, label_groups, groups = self.preprocess()
        labels = labels.astype('int')
        if coarse==True:
            new_labels_dict = dict(zip(*[np.unique(labels), groups[label_groups]]))
            labels = np.asarray([new_labels_dict[x] for x in labels])
            groups, labels = np.unique(labels, return_inverse=True)
            cell_type = groups
        super(Zeng10X, self).__init__(
                    *GeneExpressionDataset.get_attributes_from_matrix(
                        count, labels=labels),
                    gene_names=np.char.upper(gene_names), cell_types=cell_type)
        self.labels_groups = label_groups
        self.groups = groups
    def preprocess(self):
        groups = ['Lamp5', 'VIP', 'SST_CHODL', 'SST','PVALB','IT', 'NP', 'L6B','Astro', 'OPC', 'OLIGO', 'ENDO']
        groups = [x.upper() for x in groups]
        groups = np.asarray(groups)
        if os.path.isfile(self.save_path + 'Zeng10X.svmlight'):
            count, labels = load_svmlight_file(self.save_path + 'Zeng10X.svmlight')
            cell_type = np.load(self.save_path + 'Zeng10X.celltypes.npy')
            label_groups = np.load(self.save_path + 'Zeng10X.labels_groups.npy')
            gene_names = np.load(self.save_path + 'Zeng10X.gene_names.npy')
            return(count, gene_names, labels, cell_type, label_groups, groups)
        else:
            geneid, cellid, count = get_matrix_from_h5(self.save_path + '10X_nuclei_AIBS/umi_counts.h5',
                                               'mm10-1.2.0_premrna')
            label = np.genfromtxt(self.save_path + '10X_nuclei_AIBS/cluster.membership.csv', dtype='str', delimiter=',')
            label_map = np.genfromtxt(self.save_path + '10X_nuclei_AIBS/cluster.annotation.csv', dtype='str',
                              delimiter=',')
            count = count.T.tocsr()
            label_cluster = np.asarray([x.split('"')[1] for x in label[1:, 1]])
            label_barcode = np.asarray([x.split('"')[1] for x in label[1:, 0]])
            label_dict = dict(zip(label_barcode, label_cluster))
            new_count, matched_label = TryFindCells(label_dict, cellid, count)
            matched_label = matched_label.astype('int')
            matched_label = matched_label - 1
            cell_type = np.asarray([x.split('"')[1] for x in label_map[1:, 1]])
            cell_type = [x.upper() for x in cell_type]
            label_groups = [0]*4+[1]*4 + [2] + [3]*3 + [4]*2 + [5]*10 + [6]*2 + [7]*2 + [8,9,10,11]
            # for i in np.unique(label_groups):
            #     print(np.asarray(cell_type)[np.asarray(label_groups)==i])
            # freq = np.unique(matched_label,return_counts=True)[1]
            # for i in np.unique(label_groups):
            #     print(np.sum(np.asarray(freq)[np.asarray(label_groups) == i]))
            dump_svmlight_file(new_count, matched_label, self.save_path + 'Zeng10X.svmlight')
            np.save(self.save_path + 'Zeng10X.celltypes.npy', cell_type)
            np.save(self.save_path + 'Zeng10X.gene_names.npy', geneid)
            np.save(self.save_path + 'Zeng10X.labels_groups.npy', label_groups)
            return(new_count, geneid, matched_label, cell_type, label_groups, groups)



class ZengSS2(GeneExpressionDataset):
    def __init__(self, save_path='../AIBS/', coarse=True):
        self.save_path = save_path
        count, gene_names, labels, cell_type, label_groups, groups = self.preprocess()
        labels = labels.astype('int')
        print('finished preprocessing')
        if coarse==True:
            new_labels_dict = dict(zip(*[np.unique(labels), groups[label_groups]]))
            labels = np.asarray([new_labels_dict[x] for x in labels])
            groups, labels = np.unique(labels, return_inverse=True)
            cell_type = groups
            print('finished coarse labeling')
        super(ZengSS2, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                count, labels=labels),
            gene_names=np.char.upper(gene_names), cell_types=cell_type)
        self.labels_groups = label_groups
        self.groups = groups
    def preprocess(self):
        groups = ['Lamp5', 'VIP', 'SST', 'PVALB', 'IT', 'NP', 'L6B', 'Astro', 'OPC', 'OLIGO', 'ENDO',
                  'Micro']
        groups = [x.upper() for x in groups]
        groups = np.asarray(groups)
        if os.path.isfile(self.save_path + 'ZengSS2.svmlight'):
            count, labels = load_svmlight_file(self.save_path + 'ZengSS2.svmlight')
            cell_type = np.load(self.save_path + 'ZengSS2.celltypes.npy')
            label_groups = np.load(self.save_path + 'ZengSS2.labels_groups.npy')
            gene_names = np.load(self.save_path + 'ZengSS2.gene_names.npy')
            return(count, gene_names, labels, cell_type, label_groups, groups)
        else:
            mat = np.genfromtxt(self.save_path + 'SmartSeq_nuclei_AIBS/exon.counts.csv', dtype='str', delimiter=',')
            cellid = np.asarray([x.split('"')[1] for x in mat[0, 1:]])
            geneid = np.asarray([x.split('"')[1] for x in mat[1:, 0]])
            count = mat[1:, 1:].astype('int')
            count = csr_matrix(count.T)
            label = np.genfromtxt(self.save_path + 'SmartSeq_nuclei_AIBS/cluster.membership.csv', dtype='str', delimiter=',')
            label_map = np.genfromtxt(self.save_path + 'SmartSeq_nuclei_AIBS/cluster.annotation.csv', dtype='str',
                                      delimiter=',')
            label_cluster = np.asarray([x.split('"')[1] for x in label[1:, 1]])
            label_barcode = np.asarray([x.split('"')[1] for x in label[1:, 0]])
            label_dict = dict(zip(label_barcode, label_cluster))
            new_count, matched_label = TryFindCells(label_dict, cellid, count)
            temp = np.sum(new_count, axis=0)
            new_count = new_count[:,np.asarray(temp).squeeze()!=0]
            geneid = geneid[np.asarray(temp).squeeze()!=0]
            np.save(self.save_path + 'ZengSS2.gene_names.npy', geneid)
            matched_label = matched_label.astype('int')
            matched_label = matched_label - 1
            dump_svmlight_file(new_count, matched_label, self.save_path + 'ZengSS2.svmlight')
            cell_type = np.asarray([x.split('"')[1] for x in label_map[1:, 1]])
            cell_type = [x.upper() for x in cell_type]
            np.save(self.save_path + 'ZengSS2.celltypes.npy', cell_type)
            label_groups = [0] * 5 + [1] * 6 + [2] * 7 + [3] * 5 + [4] * 9 + [5] * 2 + [6] * 5 + [7] + [8] + [
                9, 10, 11] + [4] * 3 + [6]
            np.save(self.save_path + 'ZengSS2.labels_groups.npy', label_groups)
            return(new_count, geneid, matched_label, cell_type, label_groups, groups)


