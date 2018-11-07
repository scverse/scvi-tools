use_cuda = True
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import run_model
from scvi.harmonization.utils_chenling import entropy_batch_mixing
from scvi.metrics.clustering import select_indices_evenly
from scvi.dataset.dataset import SubsetGenes

import sys
import numpy as np

plotname = 'Tech1'

from scvi.dataset.muris_tabula import TabulaMuris
dataset1 = TabulaMuris('facs')
dataset2 = TabulaMuris('droplet')
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)

gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)

seurat, batch_indices, labels, keys, stats = run_model('readSeurat', gene_dataset, dataset1, dataset2,
                                                       filename=plotname)

dataset1, dataset2, gene_dataset = SubsetGenes(dataset1, dataset2, gene_dataset, plotname)
vae,  _, _, _, _ = run_model('vae', gene_dataset, dataset1, dataset2,filename=plotname, rep='0')

KNeighbors = np.concatenate([np.arange(10, 100, 10), np.arange(100, 500, 50)])
sample = select_indices_evenly(2000, batch_indices)
seurat_BE = [entropy_batch_mixing(seurat[sample, :], batch_indices[sample], n_neighbors=k) for k in KNeighbors]
vae_BE = [entropy_batch_mixing(vae[sample, :], batch_indices[sample], n_neighbors=k) for k in KNeighbors]

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
plt.figure(figsize=(5, 5))
plt.plot(KNeighbors, seurat_BE,'b',label='Seurat')
plt.plot(KNeighbors, vae_BE,'r',label='VAE')
legend = plt.legend(loc='lower right', shadow=False)
plt.savefig("../%s/%s.BE.pdf" % (plotname,plotname))


from scvi.dataset.BICCN import *
plotname = 'Zeng'

dataset1 = Zeng10X()
dataset2 = ZengSS2()
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)

seurat, batch_indices, labels, keys, stats = run_model('readSeurat', gene_dataset, dataset1, dataset2,
                                                       filename=plotname)
dataset1, dataset2, gene_dataset = SubsetGenes(dataset1, dataset2, gene_dataset, plotname)
vae,  _, _, _, _ = run_model('vae', gene_dataset, dataset1, dataset2,filename=plotname, rep='0')

KNeighbors = np.concatenate([np.arange(10, 100, 10), np.arange(100, 500, 50)])
sample = select_indices_evenly(2000, batch_indices)
seurat_BE = [entropy_batch_mixing(seurat[sample, :], batch_indices[sample], n_neighbors=k) for k in KNeighbors]
vae_BE = [entropy_batch_mixing(vae[sample, :], batch_indices[sample], n_neighbors=k) for k in KNeighbors]

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
plt.figure(figsize=(5, 5))
plt.plot(KNeighbors, seurat_BE,'b',label='Seurat')
plt.plot(KNeighbors, vae_BE,'r',label='VAE')
legend = plt.legend(loc='lower right', shadow=False)
plt.savefig("../%s/%s.BE.pdf" % (plotname,plotname))

