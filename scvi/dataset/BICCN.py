from scvi.harmonization.utils_chenling import get_matrix_from_h5, TryFindCells
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset
from sklearn.datasets import dump_svmlight_file, load_svmlight_file
import os
from scipy.sparse import csr_matrix
import seaborn as sns

color_dictionary = dict()

macosko_regev_labels = [['Pvalb low', 'Pvalb', 'Pvalb 1', 'Pvalb 2'],
                        ['Pvalb Ex_1', 'Pvalb Ex_2', 'Pvalb Ex'],
                        ['Pvalb Astro_1', 'Pvalb Astro_2'],
                        ['L2/3 IT Astro', 'L2/3 IT Macc1', 'L2/3 IT Sla_Astro', 'L2/3 IT', 'L2/3 IT Sla',
                         'L2/3 IT Sla_Inh'],
                        ['Sst Tac2', 'Sst Myh8', 'Sst Etv1', 'Sst Chodl', 'Sst'],
                        ['L5 PT_2', 'L5 PT IT', 'L5 PT_1'],
                        ['L5 IT Tcap_1_3', 'L5 IT Tcap_2', 'L5 IT Tcap_Astro', 'L5 IT Tcap_1', 'L5 IT Tcap_L2/3',
                         'L5 IT Tcap_Foxp2', 'L5 IT Tcap_3'],
                        ['L5 IT Aldh1a7_2', 'L5 IT Aldh1a7', 'L5 IT Aldh1a7_1'],
                        ['L5 NP', 'L5 NP Slc17a8'],
                        ['L6 IT Car3', 'L6 CT Olig', 'L6 IT Maf', 'L6 IT Ntn5 Mgp', 'L6 IT Ntn5 Inpp4b'],
                        ['L6 CT Nxph2', 'L6 CT Astro', 'L6 CT', 'L6 CT Grp'],
                        ['L6b', 'L6b F2r'],
                        ['Lamp5 Sncg', 'Lamp5 Egln3', 'Lamp5 Slc35d3'],
                        ['Vip Rspo4', 'Vip Serpinf1', 'Vip'],
                        ['Astro Ex', 'Astro Aqp4'],
                        ['OPC Pdgfra'],
                        ['VLMC Osr1'],
                        ['Oligo Enpp6_1', 'Oligo Enpp6_2', 'Oligo Opalin'],
                        ['Sncg Ptprk'],
                        ['Endo Slc38a5', 'Endo Slc38a5_Peri_2', 'Endo Slc38a5_Peri_1']]

macosko_regev_colors = [sns.color_palette("Greens")[2:6],  # Pvalb
                        sns.light_palette("green", 5)[0:3],  # Pvalb Ex
                        sns.light_palette("green", 5)[3:5],  # Pvalb Astro
                        sns.light_palette("orange", 6),  # L2/3
                        sns.light_palette('red')[1:6],  # Sst
                        sns.light_palette("cyan", 3),  # L5 PT
                        sns.light_palette('purple', 8)[1:8],  # L5 IT Tcap
                        sns.light_palette('purple', 7)[4:7],  # L5 IT Aldh1a7
                        sns.light_palette("navy", 7)[3:5],  # L5 NP
                        sns.light_palette("brown", 7)[2:7],  # L6 IT
                        sns.dark_palette("brown", 8)[1:5],  # L6 CT
                        sns.dark_palette("green", 8)[5:7],  # L6
                        sns.dark_palette("yellow", 7)[1:4],  # Lamp5
                        sns.dark_palette("yellow", 7)[4:7],  # Vip
                        sns.color_palette("Paired", 4),  # Astro OPC VLMC
                        sns.color_palette('Greys', 3),  # Oligo
                        sns.dark_palette('tan'),  # sncg
                        sns.light_palette('hotpink', 3)]  # endo]

macosko_regev_labels_flattened = [l1 for l2 in macosko_regev_labels for l1 in l2]
macosko_regev_colors = [c1 for c2 in macosko_regev_colors for c1 in c2]

macosko_regev_dictionary = dict(zip(macosko_regev_labels_flattened, macosko_regev_colors))

class MacoskoDataset(GeneExpressionDataset):
    def __init__(self, save_path='../AIBS/',coarse=True):
        self.save_path = save_path
        count, labels, cell_type, gene_names,labels_groups = self.preprocess()
        labels = labels.astype('int')
        if coarse==True:
            groups = ['Pvalb', 'L2/3', 'Sst', 'L5 PT', 'L5 IT Tcap', 'L5 IT Aldh1a7', 'L5 IT Foxp2', 'L5 NP',
                      'L6 IT', 'L6 CT', 'L6 NP', 'L6b', 'Lamp5', 'Vip', 'Astro', 'OPC', 'VLMC', 'Oligo', 'Sncg', 'Endo',
                      'SMC', 'MICRO']
            groups = np.asarray([x.upper() for x in groups])
            cell_type_bygroup = np.concatenate([[x for x in cell_type if x.startswith(y)] for y in groups])
            new_labels_dict = dict(zip(cell_type_bygroup, np.arange(len(cell_type_bygroup))))
            labels = np.asarray([cell_type[x] for x in labels])
            new_labels = np.asarray([new_labels_dict[x] for x in labels])
            labels_groups = [[i for i, x in enumerate(groups) if y.startswith(x)][0] for y in cell_type_bygroup]
            coarse_labels_dict = dict(zip(np.arange(len(labels_groups)), labels_groups))
            coarse_labels = np.asarray([coarse_labels_dict[x] for x in new_labels]).astype('int')
            groups = groups[np.unique(coarse_labels)]
            mapping = dict(zip(np.unique(coarse_labels), np.arange(len(np.unique(coarse_labels)))))
            coarse_labels = np.asarray([mapping[x] for x in coarse_labels])
            cell_type = groups
            labels = coarse_labels
        super(MacoskoDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                count, labels=labels),
            gene_names=np.char.upper(gene_names), cell_types=cell_type)
        self.labels_groups = labels_groups
    def preprocess(self):
        if os.path.isfile(self.save_path + 'macosko_data.svmlight'):
            count, labels = load_svmlight_file(self.save_path + 'macosko_data.svmlight')
            cell_type = np.load(self.save_path + 'macosko_data.celltypes.npy')
            labels_groups = np.load(self.save_path + 'macosko_data.labels_groups.npy')
            gene_names = np.load(self.save_path + 'macosko_data.gene_names.npy')
            return(count, labels, cell_type, gene_names, labels_groups)
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
                new_label = np.repeat(0, len(matched_label))
                for i, x in enumerate(np.unique(matched_label)):
                    new_label[matched_label == x] = i
                cell_type = [label_map[x] for x in np.unique(matched_label)]
                dataset = GeneExpressionDataset(
                    *GeneExpressionDataset.get_attributes_from_matrix(new_count, labels=new_label),
                    gene_names=geneid, cell_types=cell_type)
                print(dataset.X.shape, len(dataset.labels))
                if len(macosko_data) > 0:
                    macosko_data = GeneExpressionDataset.concat_datasets(macosko_data, dataset)
                else:
                    macosko_data = dataset
            dataset = macosko_data
            cell_type = dataset.cell_types
            groups = ['Pvalb', 'L2/3', 'Sst', 'L5 PT', 'L5 IT Tcap', 'L5 IT Aldh1a7', 'L5 IT Foxp2', 'L5 NP',
                      'L6 IT', 'L6 CT', 'L6 NP', 'L6b', 'Lamp5', 'Vip', 'Astro', 'OPC', 'VLMC', 'Oligo', 'Sncg', 'Endo',
                      'SMC', 'MICRO']
            cell_type = [x.upper() for x in cell_type]
            groups = [x.upper() for x in groups]
            labels = np.asarray([cell_type[x] for x in np.concatenate(dataset.labels)])
            cell_type_bygroup = np.concatenate([[x for x in cell_type if x.startswith(y)] for y in groups])
            new_labels_dict = dict(zip(cell_type_bygroup, np.arange(len(cell_type_bygroup))))
            new_labels = np.asarray([new_labels_dict[x] for x in labels])
            labels_groups = [[i for i, x in enumerate(groups) if y.startswith(x)][0] for y in cell_type_bygroup]
            dump_svmlight_file(dataset.X, new_labels, self.save_path + 'macosko_data.svmlight')
            np.save(self.save_path + 'macosko_data.celltypes.npy', cell_type_bygroup)
            np.save(self.save_path + 'macosko_data.gene_names.npy', dataset.gene_names)
            np.save(self.save_path + 'macosko_data.labels_groups.npy', labels_groups)
            return (dataset.X, new_labels, cell_type_bygroup, dataset.gene_names, labels_groups)


class RegevDataset(GeneExpressionDataset):
    def __init__(self, save_path='../AIBS/',coarse=True):
        self.save_path = save_path
        count, labels, cell_type, gene_names, labels_groups = self.preprocess()
        labels = labels.astype('int')
        if coarse==True:
            groups = ['Pvalb', 'L2/3', 'Sst', 'L5 PT', 'L5 IT Tcap', 'L5 IT Aldh1a7', 'L5 IT Foxp2', 'L5 NP',
                      'L6 IT', 'L6 CT', 'L6 NP', 'L6b', 'Lamp5', 'Vip', 'Astro', 'OPC', 'VLMC', 'Oligo', 'Sncg', 'Endo',
                      'SMC', 'MICRO']
            groups = np.asarray([x.upper() for x in groups])
            cell_type_bygroup = np.concatenate([[x for x in cell_type if x.startswith(y)] for y in groups])
            new_labels_dict = dict(zip(cell_type_bygroup, np.arange(len(cell_type_bygroup))))
            labels = np.asarray([cell_type[x] for x in labels])
            new_labels = np.asarray([new_labels_dict[x] for x in labels])
            labels_groups = [[i for i, x in enumerate(groups) if y.startswith(x)][0] for y in cell_type_bygroup]
            coarse_labels_dict = dict(zip(np.arange(len(labels_groups)), labels_groups))
            coarse_labels = np.asarray([coarse_labels_dict[x] for x in new_labels]).astype('int')
            groups = groups[np.unique(coarse_labels)]
            mapping = dict(zip(np.unique(coarse_labels), np.arange(len(np.unique(coarse_labels)))))
            coarse_labels = np.asarray([mapping[x] for x in coarse_labels])
            cell_type = groups
            labels = coarse_labels
        super(RegevDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                count, labels=labels),
            gene_names=np.char.upper(gene_names), cell_types=cell_type)
        self.labels_groups = labels_groups

    def preprocess(self):
        if os.path.isfile(self.save_path + 'regev_data.svmlight'):
            count, labels = load_svmlight_file(self.save_path + 'regev_data.svmlight')
            cell_type = np.load(self.save_path + 'regev_data.celltypes.npy')
            gene_names = np.load(self.save_path + 'regev_data.gene_names.npy')
            labels_groups = np.load(self.save_path + 'regev_data.labels_groups.npy')
            return(count, labels, cell_type, gene_names, labels_groups)
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
                new_label = np.repeat(0, len(matched_label))
                for i, x in enumerate(np.unique(matched_label)):
                    new_label[matched_label == x] = i
                cell_type = [label_map[x] for x in np.unique(matched_label)]
                dataset = GeneExpressionDataset(
                    *GeneExpressionDataset.get_attributes_from_matrix(new_count, labels=new_label),
                    gene_names=geneid, cell_types=cell_type)
                print(dataset.X.shape, len(dataset.labels))
                if len(regev_data) > 0:
                    regev_data = GeneExpressionDataset.concat_datasets(regev_data, dataset)
                else:
                    regev_data = dataset
            dataset = regev_data
            cell_type = dataset.cell_types
            groups = ['Pvalb', 'L2/3', 'Sst', 'L5 PT', 'L5 IT Tcap', 'L5 IT Aldh1a7', 'L5 IT Foxp2', 'L5 NP',
                      'L6 IT', 'L6 CT', 'L6 NP', 'L6b', 'Lamp5', 'Vip', 'Astro', 'OPC', 'VLMC', 'Oligo', 'Sncg', 'Endo',
                      'SMC', 'MICRO']
            cell_type = [x.upper() for x in cell_type]
            groups = [x.upper() for x in groups]
            labels = np.asarray([cell_type[x] for x in np.concatenate(dataset.labels)])
            cell_type_bygroup = np.concatenate([[x for x in cell_type if x.startswith(y)] for y in groups])
            new_labels_dict = dict(zip(cell_type_bygroup, np.arange(len(cell_type_bygroup))))
            new_labels = np.asarray([new_labels_dict[x] for x in labels])
            labels_groups = [[i for i, x in enumerate(groups) if y.startswith(x)][0] for y in cell_type_bygroup]
            dump_svmlight_file(dataset.X, new_labels, self.save_path+'regev_data.svmlight')
            np.save(self.save_path + 'regev_data.celltypes.npy', cell_type_bygroup)
            np.save(self.save_path + 'regev_data.gene_names.npy', dataset.gene_names)
            np.save(self.save_path + 'regev_data.labels_groups.npy', labels_groups)
            return (dataset.X, new_labels, cell_type_bygroup, dataset.gene_names, labels_groups)


def FindGroup(celltype,groups):
    res = [i for i,x in enumerate(groups)  if celltype.startswith(x)]
    return(res)


class Zeng10X(GeneExpressionDataset):
    def __init__(self, save_path='../AIBS/', coarse=True):
        self.save_path = save_path
        count, gene_names, labels, cell_type, labels_groups, groups = self.preprocess()
        assert len(labels_groups) == len(cell_type)
        labels = labels.astype('int')
        if coarse==True:
            new_labels_dict = dict(zip(*[np.unique(labels), groups[labels_groups]]))
            labels = np.asarray([new_labels_dict[x] for x in labels])
            groups, labels = np.unique(labels, return_inverse=True)
            cell_type = groups
        super(Zeng10X, self).__init__(
                    *GeneExpressionDataset.get_attributes_from_matrix(
                        count, labels=labels),
                    gene_names=np.char.upper(gene_names), cell_types=cell_type)
        self.labels_groups = labels_groups
        self.groups = groups
    def preprocess(self):
        groups = ['Pvalb', 'L2/3', 'Sst', 'L5 PT', 'L5 IT Tcap', 'L5 IT Aldh1a7', 'L5 IT Foxp2', 'L5 NP',
                  'L6 IT', 'L6 CT', 'L6 NP', 'L6b', 'Lamp5', 'Vip', 'Astro', 'OPC', 'VLMC', 'Oligo', 'Sncg', 'Endo',
                  'SMC', 'MICRO']
        groups = [x.upper() for x in groups]
        groups = np.asarray(groups)
        if os.path.isfile(self.save_path + 'Zeng10X.svmlight'):
            count, labels = load_svmlight_file(self.save_path + 'Zeng10X.svmlight')
            cell_type = np.load(self.save_path + 'Zeng10X.celltypes.npy')
            gene_names = np.load(self.save_path + 'Zeng10X.gene_names.npy')
            labels_groups = np.load(self.save_path + 'Zeng10X.labels_groups.npy')
            return(count, gene_names, labels, cell_type, labels_groups, groups)
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
            cell_type[cell_type=='IT L6b'] = 'L6B IT'
            cell_type = np.asarray([x.upper() for x in cell_type])
            labels_groups = [FindGroup(y,groups) for y in cell_type]
            labels_groups = np.concatenate(labels_groups)
            dump_svmlight_file(new_count, matched_label, self.save_path + 'Zeng10X.svmlight')
            np.save(self.save_path + 'Zeng10X.celltypes.npy', cell_type)
            np.save(self.save_path + 'Zeng10X.gene_names.npy', geneid)
            np.save(self.save_path + 'Zeng10X.labels_groups.npy', labels_groups)
            return(new_count, geneid, matched_label, cell_type, labels_groups, groups)


class ZengSS2(GeneExpressionDataset):
    def __init__(self, save_path='../AIBS/', coarse=True):
        self.save_path = save_path
        count, gene_names, labels, cell_type, labels_groups, groups = self.preprocess()
        assert len(labels_groups) == len(cell_type)
        labels = labels.astype('int')
        print('finished preprocessing')
        if coarse==True:
            new_labels_dict = dict(zip(*[np.unique(labels), groups[labels_groups]]))
            labels = np.asarray([new_labels_dict[x] for x in labels])
            groups, labels = np.unique(labels, return_inverse=True)
            cell_type = groups
            print('finished coarse labeling')
        super(ZengSS2, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                count, labels=labels),
            gene_names=np.char.upper(gene_names), cell_types=cell_type)
        self.labels_groups = labels_groups
        self.groups = groups
    def preprocess(self):
        groups = ['Pvalb', 'L2/3', 'Sst', 'L5 PT', 'L5 IT Tcap', 'L5 IT Aldh1a7', 'L5 IT Foxp2', 'L5 NP',
          'L6 IT', 'L6 CT', 'L6 NP', 'L6b', 'Lamp5', 'Vip', 'Astro', 'OPC', 'VLMC', 'Oligo', 'Sncg', 'Endo',
          'SMC', 'MICRO']
        groups = [x.upper() for x in groups]
        groups = np.asarray(groups)
        if os.path.isfile(self.save_path + 'ZengSS2.cells.svmlight'):
            count, labels = load_svmlight_file(self.save_path + 'ZengSS2.cells.svmlight')
            cell_type = np.load(self.save_path + 'ZengSS2.cells.celltypes.npy')
            gene_names = np.load(self.save_path + 'ZengSS2.cells.gene_names.npy')
            labels_groups = np.load(self.save_path + 'ZengSS2.cells.labels_groups.npy')
            return(count, gene_names, labels, cell_type, labels_groups, groups)
        else:
            mat = np.genfromtxt(self.save_path + 'SmartSeq_cells_AIBS/exon.counts.csv', dtype='str', delimiter=',')
            cellid = np.asarray([x.split('"')[1] for x in mat[0, 1:]])
            geneid = np.asarray([x.split('"')[1] for x in mat[1:, 0]])
            count = mat[1:, 1:].astype('int')
            count = csr_matrix(count.T)
            label = np.genfromtxt(self.save_path + 'SmartSeq_cells_AIBS/cluster.membership.csv', dtype='str', delimiter=',')
            label_map = np.genfromtxt(self.save_path + 'SmartSeq_cells_AIBS/cluster.annotation.csv', dtype='str',
                                      delimiter=',')
            label_cluster = np.asarray([x.split('"')[1] for x in label[1:, 1]])
            label_barcode = np.asarray([x.split('"')[1] for x in label[1:, 0]])
            label_dict = dict(zip(label_barcode, label_cluster))
            new_count, matched_label = TryFindCells(label_dict, cellid, count)
            temp = np.sum(new_count, axis=0)
            new_count = new_count[:,np.asarray(temp).squeeze()!=0]
            geneid = geneid[np.asarray(temp).squeeze()!=0]
            matched_label = matched_label.astype('int')
            matched_label = matched_label - 1
            cell_type = np.asarray([x.split('"')[1] for x in label_map[1:, 1]])
            cell_type[cell_type=='IT L6b'] = 'L6B IT'
            cell_type = np.asarray([x.upper() for x in cell_type])
            labels_groups = [FindGroup(y,groups) for y in cell_type]
            labels_groups = np.concatenate(labels_groups)
            dump_svmlight_file(new_count, matched_label, self.save_path + 'ZengSS2.cells.svmlight')
            np.save(self.save_path + 'ZengSS2.cells.gene_names.npy', geneid)
            np.save(self.save_path + 'ZengSS2.cells.celltypes.npy', cell_type)
            np.save(self.save_path + 'ZengSS2.cells.labels_groups.npy', labels_groups)
            return(new_count, geneid, matched_label, cell_type, labels_groups, groups)

