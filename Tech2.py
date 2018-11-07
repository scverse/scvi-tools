use_cuda = True
from scvi.dataset.BICCN import *
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import CompareModels

import sys
models = str(sys.argv[1])
plotname = 'Zeng'


dataset1 = Zeng10X()
dataset2 = ZengSS2()
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
