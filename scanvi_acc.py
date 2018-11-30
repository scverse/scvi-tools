use_cuda = True
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.dataset.dataset import SubsetGenes
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from scvi.harmonization.utils_chenling import trainSCANVI
import numpy as np


def SCANVI_acc(gene_dataset:GeneExpressionDataset,dataset1,dataset2, plotname: str,rep='1'):
    fname = '../%s/%s.scanvi_acc.txt'%(plotname,plotname)
    methods = ['scanvi','scanvi1','scanvi2']
    f = open(fname, "w+")
    f.write('method\t' +  "%s\t" * len(gene_dataset.cell_types) % tuple(gene_dataset.cell_types) + "\n")
    dataset1, dataset2, gene_dataset = SubsetGenes(dataset1, dataset2, gene_dataset, plotname)
    for i,method in enumerate(methods):
        trainer_scanvi = trainSCANVI(gene_dataset,method,plotname,rep)
        unlabelled_idx = trainer_scanvi.unlabelled_set.indices
        full = trainer_scanvi.create_posterior(trainer_scanvi.model, gene_dataset, indices=np.arange(len(gene_dataset)))
        labels, labels_pred = full.sequential().compute_predictions()
        acc = [np.mean(labels_pred[unlabelled_idx][labels[unlabelled_idx] == i] == i) for i in np.unique(labels)]
        f.write(method + "\t" + "%.4f\t" * len(acc) % tuple(acc) + "\n")
    f.close()

from scvi.dataset.muris_tabula import TabulaMuris
dataset1 = TabulaMuris('facs')
dataset2 = TabulaMuris('droplet')
plotname = 'Tech1'
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)
SCANVI_acc(gene_dataset, dataset1, dataset2, plotname)


from scvi.harmonization.utils_chenling import get_matrix_from_dir,assign_label
from scvi.dataset.pbmc import PbmcDataset
plotname = 'Easy1'
dataset1 = PbmcDataset(filter_out_de_genes=False)
dataset1.update_cells(dataset1.batch_indices.ravel()==0)
dataset1.subsample_genes(dataset1.nb_genes)
count, geneid, cellid = get_matrix_from_dir('cite')
count = count.T.tocsr()
seurat = np.genfromtxt('../cite/cite.seurat.labels', dtype='str', delimiter=',')
cellid = np.asarray([x.split('-')[0] for x in cellid])
labels_map = [0, 0, 1, 2, 3, 4, 5, 6]
labels = seurat[1:, 4]
cell_type = ['CD4 T cells', 'NK cells', 'CD14+ Monocytes', 'B cells','CD8 T cells', 'FCGR3A+ Monocytes', 'Other']
dataset2 = assign_label(cellid, geneid, labels_map, count, cell_type, seurat)
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1,dataset2)
SCANVI_acc(gene_dataset, dataset1, dataset2, plotname)

from copy import deepcopy
from scvi.dataset.scanorama import DatasetSCANORAMA
plotname='Tech4'
dirs = ['data/pancreas/pancreas_inDrop', 'data/pancreas/pancreas_multi_celseq2_expression_matrix', 'data/pancreas/pancreas_multi_celseq_expression_matrix', 'data/pancreas/pancreas_multi_fluidigmc1_expression_matrix', 'data/pancreas/pancreas_multi_smartseq2_expression_matrix']

datasets = [DatasetSCANORAMA(d) for d in dirs]

labels = (open('/data/scanorama/data/cell_labels/pancreas_cluster.txt').read().rstrip().split())

all_dataset = GeneExpressionDataset.concat_datasets(*datasets)
batch_id = all_dataset.batch_indices.ravel()

all_dataset = GeneExpressionDataset.concat_datasets(datasets[0],datasets[1])
all_dataset.cell_types,all_dataset.labels = np.unique(np.asarray(labels)[(batch_id==0)+(batch_id==1)],return_inverse=True)
all_dataset.labels = all_dataset.labels.reshape(len(all_dataset.labels),1)
all_dataset.n_labels = len(np.unique(all_dataset.labels))
all_dataset.subsample_genes(all_dataset.nb_genes)
dataset1 = deepcopy(all_dataset)
dataset1.update_cells(dataset1.batch_indices.ravel()==0)
dataset2 = deepcopy(all_dataset)
dataset2.update_cells(dataset2.batch_indices.ravel()==1)
SCANVI_acc(all_dataset, dataset1, dataset2, plotname)


from scvi.dataset.dataset import GeneExpressionDataset
from scvi.dataset.MouseBrain import DentateGyrus10X,DentateGyrusC1
plotname = 'Tech3'
dataset1 = DentateGyrus10X()
dataset1.subsample_genes(dataset1.nb_genes)
dataset2 = DentateGyrusC1()
dataset2.subsample_genes(dataset2.nb_genes)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1,dataset2)
SCANVI_acc(gene_dataset, dataset1, dataset2, plotname)

import pandas as pd
from copy import deepcopy
from scvi.harmonization.utils_chenling import get_matrix_from_dir
class MSDataset(GeneExpressionDataset):
    def __init__(self):
        count, labels, cell_type, gene_names,cellid = self.preprocess()
        super(MSDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                count.tocsr(), labels=labels),
            gene_names=np.char.upper(gene_names), cell_types=cell_type)
    def preprocess(self):
        print("Preprocessing data")
        count, geneid, cellid = get_matrix_from_dir('CSF',storage='/data/')
        print("Finish preprocessing data")
        labels = np.repeat(0,len(cellid))
        cell_type = ['unlabelled']
        return count,labels,cell_type,geneid,cellid


gene_dataset = MSDataset()
celltypes = pd.read_csv('/data/CSF/celltypes.txt')
celltype = np.asarray(celltypes['celltypes'])
donor = np.asarray(celltypes['donor'])

gene_dataset.cell_types,gene_dataset.labels = np.unique(celltype,return_inverse=True)
gene_dataset.batch_names,gene_dataset.batch_indices = np.unique(donor,return_inverse=True)
gene_dataset.batch_indices = gene_dataset.batch_indices.reshape(len(gene_dataset.batch_indices),1)

dataset1 = deepcopy(gene_dataset)
dataset1.update_cells(dataset1.batch_indices.ravel()==3)
dataset2 = deepcopy(gene_dataset)
dataset2.update_cells(dataset2.batch_indices.ravel()==7)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1,dataset2)
SCANVI_acc(gene_dataset, dataset1, dataset2, 'CSF')
#
