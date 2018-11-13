use_cuda = True
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import run_model
from scvi.harmonization.utils_chenling import entropy_batch_mixing
from scvi.metrics.clustering import select_indices_evenly
from scvi.dataset.dataset import SubsetGenes
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def BE_curve(dataset1: GeneExpressionDataset, dataset2: GeneExpressionDataset, plotname: str,recompute=False):
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    fname = '../%s/%s.BE.txt'%(plotname,plotname)
    colors = ('r', 'g', 'b', 'y', 'm', 'c')
    methods = ['vae', 'scanvi', 'readSeurat', 'Combat', 'MNN', 'PCA']
    model_names = ['scVI', 'SCAN-VI', 'CCA', 'Combat', 'MNN', 'PCA']
    if (not os.path.isfile(fname)) or recompute==True:
        f = open(fname, "w+")
        gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)
        _, batch_indices, labels, keys, stats = run_model('readSeurat', gene_dataset, dataset1, dataset2,
                                                               filename=plotname)
        plt.figure(figsize=(5, 5))
        dataset1, dataset2, gene_dataset = SubsetGenes(dataset1, dataset2, gene_dataset, plotname)
        for i,method in enumerate(methods):
            latent,  _, _, _, _ = run_model(method, gene_dataset, dataset1, dataset2,filename=plotname, rep='0')
            KNeighbors = np.concatenate([np.arange(10, 100, 10), np.arange(100, 500, 50)])
            sample = select_indices_evenly(np.min(np.unique(batch_indices,return_counts=True)[1]), batch_indices)
            BE = [entropy_batch_mixing(latent[sample, :], batch_indices[sample], n_neighbors=k) for k in KNeighbors]
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


from scvi.dataset.muris_tabula import TabulaMuris
dataset1 = TabulaMuris('facs')
dataset2 = TabulaMuris('droplet')
plotname = 'Tech1'
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)
BE_curve(dataset1, dataset2,plotname,recompute=True)

from scvi.dataset.BICCN import *
plotname = 'Zeng'
dataset1 = Zeng10X()
dataset2 = ZengSS2()
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)
BE_curve(dataset1, dataset2,plotname)

plotname = 'Macosko_Regev'
dataset1 = MacoskoDataset()
dataset2 = RegevDataset()
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset1.nb_genes)
BE_curve(dataset1, dataset2,plotname)

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
BE_curve(dataset1, dataset2,plotname)

