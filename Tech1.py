use_cuda = True
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import CompareModels
import sys
import numpy as np

models = str(sys.argv[1])
plotname = 'Tech1'

from scvi.dataset.muris_tabula import TabulaMuris
dataset1 = TabulaMuris('facs')
dataset2 = TabulaMuris('droplet')
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)

gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)

CompareModels(gene_dataset, dataset1, dataset2, plotname, models)
# 438.05354394990854
#
#
# from scvi.harmonization.utils_chenling import run_model
# from scvi.harmonization.utils_chenling import scmap_eval
# latent, batch_indices, labels, keys, stats = run_model('scmap', gene_dataset, dataset1, dataset2, filename=plotname)
# pred1 = latent
# pred2 = stats
# res1 = scmap_eval(pred1, labels[batch_indices == 1],labels)
# res2 = scmap_eval(pred2, labels[batch_indices == 0],labels)
# print(res1)
# print(res2)
# temp1, batch_indices, labels, keys, temp2 = run_model('scmap', gene_dataset, dataset1, dataset2, filename=plotname)
# pred1 = temp1
# pred2 = temp2
# res1 = scmap_eval(pred1, labels[batch_indices == 1],labels)
# res2 = scmap_eval(pred2, labels[batch_indices == 0],labels)
