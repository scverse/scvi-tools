use_cuda = True
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import run_model
from scvi.harmonization.utils_chenling import entropy_batch_mixing
from scvi.metrics.clustering import select_indices_evenly
from scvi.dataset.dataset import SubsetGenes
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt

def BE_curve(gene_dataset:GeneExpressionDataset, dataset1, dataset2, plotname: str,recompute=False):
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    fname = '../%s/%s.BE.txt'%(plotname,plotname)
    colors = ('r', 'g', 'b', 'y', 'm', 'c')
    methods = ['vae', 'scanvi', 'readSeurat', 'Combat', 'MNN', 'PCA']
    model_names = ['scVI', 'SCAN-VI', 'CCA', 'Combat', 'MNN', 'PCA']
    if (not os.path.isfile(fname)) or recompute==True:
        f = open(fname, "w+")
        _, batch_indices, labels, keys, stats = run_model('readSeurat', gene_dataset, dataset1, dataset2,
                                                               filename=plotname)
        plt.figure(figsize=(5, 5))
        dataset1, dataset2, gene_dataset = SubsetGenes(dataset1, dataset2, gene_dataset, plotname)
        for i,method in enumerate(methods):
            latent,  _, _, _, _ = run_model(method, gene_dataset, dataset1, dataset2,filename=plotname, rep='0')
            KNeighbors = np.concatenate([np.arange(10, 100, 10), np.arange(100, 500, 50)])
            sample = select_indices_evenly(2000, batch_indices)
            BE = [entropy_batch_mixing(latent[sample, :], batch_indices[sample], n_neighbors=k,n_samples_per_pool=500) for k in KNeighbors]
            plt.plot(KNeighbors, BE, colors[i], label=model_names[i])
            f.write(method + "\t" + "%.4f\t"*len(BE)%tuple(BE) + "\n")
        plt.legend(loc='lower right', shadow=False)
        plt.savefig("../%s/%s.BE.pdf" % (plotname,plotname))
        f.close()
    else:
        import pandas as pd
        stats = pd.read_table(fname, delim_whitespace=True)
        res = []
        for x in np.unique(methods):
            stat = np.mean(np.asarray(stats[methods == x])[:, 1:], axis=0)
            res.append(stat)
        res = np.asarray(res)
        sorted_res = []
        model_order = ['vae', 'scanvi', 'readSeurat', 'Combat', 'MNN', 'PCA']
        model_names = ['scVI', 'SCAN-VI', 'CCA', 'Combat', 'MNN', 'PCA']
        for x in model_order:
            sorted_res.append(res[methods == x, :])
        sorted_res = np.asarray(sorted_res).squeeze()
        import matplotlib
        KNeighbors = np.concatenate([np.arange(10, 100, 10), np.arange(100, 500, 50)])
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42
        import matplotlib.pyplot as plt
        plt.figure(figsize=(5, 5))
        colors = ('r', 'g', 'b', 'y', 'm', 'c')
        for i, x in enumerate(model_names):
            plt.plot(KNeighbors, sorted_res[:, i], colors[i], label=x)
        plt.legend(loc='lower right', shadow=False)
        plt.savefig("../%s/%s.BE.pdf" % (plotname, plotname))


# from scvi.dataset.muris_tabula import TabulaMuris
# dataset1 = TabulaMuris('facs')
# dataset2 = TabulaMuris('droplet')
# plotname = 'Tech1'
# dataset1.subsample_genes(dataset1.nb_genes)
# dataset2.subsample_genes(dataset2.nb_genes)
# gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)
# BE_curve(dataset1, dataset2, gene_dataset, plotname,recompute=True)
#
# from scvi.dataset.BICCN import *
# plotname = 'Zeng'
# dataset1 = Zeng10X()
# dataset2 = ZengSS2()
# dataset1.subsample_genes(dataset1.nb_genes)
# dataset2.subsample_genes(dataset2.nb_genes)
# BE_curve(dataset1, dataset2,plotname)
#
# plotname = 'Macosko_Regev'
# dataset1 = MacoskoDataset()
# dataset2 = RegevDataset()
# dataset1.subsample_genes(dataset1.nb_genes)
# dataset2.subsample_genes(dataset1.nb_genes)
# BE_curve(dataset1, dataset2,plotname)

# from scvi.harmonization.utils_chenling import get_matrix_from_dir,assign_label
# from scvi.dataset.pbmc import PbmcDataset
# plotname = 'Easy1'
# dataset1 = PbmcDataset(filter_out_de_genes=False)
# dataset1.update_cells(dataset1.batch_indices.ravel()==0)
# dataset1.subsample_genes(dataset1.nb_genes)
# count, geneid, cellid = get_matrix_from_dir('cite')
# count = count.T.tocsr()
# seurat = np.genfromtxt('../cite/cite.seurat.labels', dtype='str', delimiter=',')
# cellid = np.asarray([x.split('-')[0] for x in cellid])
# labels_map = [0, 0, 1, 2, 3, 4, 5, 6]
# labels = seurat[1:, 4]
# cell_type = ['CD4 T cells', 'NK cells', 'CD14+ Monocytes', 'B cells','CD8 T cells', 'FCGR3A+ Monocytes', 'Other']
# dataset2 = assign_label(cellid, geneid, labels_map, count, cell_type, seurat)
# dataset1.subsample_genes(dataset1.nb_genes)
# dataset2.subsample_genes(dataset2.nb_genes)
# BE_curve(dataset1, dataset2,plotname)
#
# from scvi.dataset.scanorama import DatasetSCANORAMA
# plotname='Tech4'
# dirs = ['data/pancreas/pancreas_inDrop', 'data/pancreas/pancreas_multi_celseq2_expression_matrix', 'data/pancreas/pancreas_multi_celseq_expression_matrix', 'data/pancreas/pancreas_multi_fluidigmc1_expression_matrix', 'data/pancreas/pancreas_multi_smartseq2_expression_matrix']
#
# datasets = [DatasetSCANORAMA(d) for d in dirs]
#
# labels = (open('/data/scanorama/data/cell_labels/pancreas_cluster.txt').read().rstrip().split())
#
# all_dataset = GeneExpressionDataset.concat_datasets(*datasets)
# batch_id = all_dataset.batch_indices.ravel()
#
# all_dataset = GeneExpressionDataset.concat_datasets(datasets[0],datasets[1])
# all_dataset.cell_types,all_dataset.labels = np.unique(np.asarray(labels)[(batch_id==0)+(batch_id==1)],return_inverse=True)
# all_dataset.labels = all_dataset.labels.reshape(len(all_dataset.labels),1)
# all_dataset.n_labels = len(np.unique(all_dataset.labels))
# all_dataset.subsample_genes(all_dataset.nb_genes)
# dataset1 = deepcopy(all_dataset)
# dataset1.update_cells(dataset1.batch_indices.ravel()==0)
# dataset2 = deepcopy(all_dataset)
# dataset2.update_cells(dataset2.batch_indices.ravel()==1)
# BE_curve(all_dataset, plotname,recompute=True)
#
# from scvi.dataset.dataset import GeneExpressionDataset
# from scvi.dataset.MouseBrain import DentateGyrus10X,DentateGyrusC1
# plotname = 'Tech3'
# dataset1 = DentateGyrus10X()
# dataset1.subsample_genes(dataset1.nb_genes)
# dataset2 = DentateGyrusC1()
# dataset2.subsample_genes(dataset2.nb_genes)
# gene_dataset = GeneExpressionDataset.concat_datasets(dataset1,dataset2)
# BE_curve(gene_dataset,dataset1,dataset2, plotname,recompute=True)

import pandas as pd
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

BE_curve(gene_dataset,dataset1,dataset2, 'CSF',recompute=True)
