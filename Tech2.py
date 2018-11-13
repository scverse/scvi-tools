use_cuda = True
from scvi.dataset.BICCN import *
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import CompareModels

import sys
models = str(sys.argv[1])
plotname = 'Zeng'


dataset1 = Zeng10X(coarse=False)
dataset2 = ZengSS2(coarse=False)
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)


# labels, props = np.unique(gene_dataset.labels,return_counts=True)
# props = props/np.sum(props)
# prop1 = np.asarray([
#     np.mean(gene_dataset.labels.ravel()[gene_dataset.batch_indices.ravel()==0] == i)
# for i in labels])
# prop2 = np.asarray([
#     np.mean(gene_dataset.labels.ravel()[gene_dataset.batch_indices.ravel()==1] == i)
# for i in labels])
# f = open('../%s/celltypeprop.txt'%plotname,mode='w')
# f.write("\t".join(gene_dataset.cell_types[labels])+"\n")
# f.write(("%.6f\t"*len(props))%tuple(props)+"\n")
# f.write(("%.6f\t"*len(props))%tuple(prop1)+"\n")
# f.write(("%.6f\t"*len(prop2))%tuple(prop2)+"\n")
# f.close()
CompareModels(gene_dataset, dataset1, dataset2, plotname, models)
#
#
# from scvi.harmonization.utils_chenling import SubsetGenes,run_model,eval_latent
# dataset1, dataset2, gene_dataset = SubsetGenes(dataset1, dataset2, gene_dataset, plotname)
# gene_dataset.subsample_genes(1000)
# latent, batch_indices, labels, keys, stats = run_model('scanvi2', gene_dataset, dataset1, dataset2,
#                                                                    filename=plotname, rep='1')
#
# res_knn, res_knn_partial, res_kmeans, res_kmeans_partial = \
#     eval_latent(batch_indices=batch_indices, labels=labels, latent=latent, keys=keys,
#                 labelled_idx=labelled_idx, unlabelled_idx=unlabelled_idx,
#                 plotname='Zeng.scanvi2.diag', plotting=True, partial_only=False)
