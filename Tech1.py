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
# genes = np.genfromtxt('../Seurat_data/'+plotname+'.CCA.genes.txt')
# genes = genes.astype('int')
# gene_dataset.X = gene_dataset.X[:,genes]
# gene_dataset.update_genes(genes)
# cells = np.genfromtxt('../Seurat_data/'+plotname+'.CCA.cells.txt')

ngenes = 1000

import pandas as pd
genes1 = pd.read_table('../Seurat_data/'+plotname+'.1.hvg_info.csv',delimiter=',')
geneid1 =np.asarray([x.replace('gene_','') for x in genes1[genes1.keys()[0]]]).astype('int')
genenames1 = genes1['genename']
genes2 = pd.read_table('../Seurat_data/'+plotname+'.2.hvg_info.csv',delimiter=',')
geneid2 =np.asarray([x.replace('gene_','') for x in genes2[genes2.keys()[0]]]).astype('int')
genenames2 = genes2['genename']
assert np.sum(np.asarray(genenames1)==gene_dataset.gene_names)==len(gene_dataset.gene_names)
assert np.sum(np.asarray(genenames2)==gene_dataset.gene_names)==len(gene_dataset.gene_names)
geneid = np.union1d(geneid1[:ngenes],geneid2[:ngenes])-1
gene_dataset.X = gene_dataset.X[:,geneid]
gene_dataset.update_genes(geneid)
dataset1.X = dataset1.X[:,geneid]
dataset1.update_genes(geneid)
dataset2.X = dataset2.X[:,geneid]
dataset2.update_genes(geneid)

CompareModels(gene_dataset, dataset1, dataset2, plotname, models)
# 438.05354394990854

