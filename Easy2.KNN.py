use_cuda = True
from scvi.dataset.BICCN import *
from scvi.harmonization.utils_chenling import run_model
from scvi.harmonization.utils_chenling import KNNJaccardIndex
from scvi.dataset.dataset import SubsetGenes

import numpy as np

plotname = 'Macosko_Regev'

dataset1 = MacoskoDataset()
dataset2 = RegevDataset()
dataset1.subsample_genes(dataset1.nb_genes)
dataset2.subsample_genes(dataset2.nb_genes)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)

seurat1 = np.genfromtxt('../Seurat_data/' + plotname + '.1.CCA.txt')
seurat2 = np.genfromtxt('../Seurat_data/' + plotname + '.2.CCA.txt')
seurat, batch_indices, labels, keys, stats = run_model('readSeurat', gene_dataset, dataset1, dataset2,
                                                       filename=plotname)

dataset1, dataset2, gene_dataset = SubsetGenes(dataset1, dataset2, gene_dataset, plotname)
vae1, _, _, _, _ = run_model('vae', dataset1, 0, 0, filename=plotname, rep='vae1')
vae2, _, _, _, _ = run_model('vae', dataset2, 0, 0, filename=plotname, rep='vae2')
vae,  _, _, _, _ = run_model('vae', gene_dataset, dataset1, dataset2,filename=plotname, rep='0')

KNeighbors = np.concatenate([np.arange(10, 100, 10), np.arange(100, 500, 50)])
seurat_seurat = [KNNJaccardIndex(seurat1, seurat2, seurat, batch_indices, k, True)[0] for k in KNeighbors]
vae_seurat = [KNNJaccardIndex(vae1, vae2, seurat, batch_indices, k, True)[0] for k in KNeighbors]
vae_vae = [KNNJaccardIndex(vae1, vae2, vae, batch_indices, k, True)[0] for k in KNeighbors]
seurat_vae = [KNNJaccardIndex(seurat1, seurat2, vae, batch_indices, k, True)[0] for k in KNeighbors]

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 10))
plt.plot(KNeighbors, seurat_seurat,'r',label='Seurat_Seurat')
plt.plot(KNeighbors, vae_vae,'b',label='VAE_VAE')
plt.plot(KNeighbors, vae_seurat,'g',label='VAE_Seurat')
plt.plot(KNeighbors, seurat_vae,'y',label='Seurat_VAE')
legend = plt.legend(loc='lower right', shadow=False)
plt.savefig("../%s/%s.KNN.pdf" % (plotname,plotname))
