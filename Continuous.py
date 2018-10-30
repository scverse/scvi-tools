from scvi.dataset.paul_tusi import Paul,Tusi
from scvi.dataset.dataset import GeneExpressionDataset
import numpy as np
use_cuda = True
from scvi.harmonization.utils_chenling import run_model
from scvi.dataset.dataset import SubsetGenes
from scvi.harmonization.utils_chenling import CompareModels
from scvi.harmonization.utils_chenling import KNNJaccardIndex
from scipy.stats import entropy
from sklearn.neighbors import NearestNeighbors

plotname = 'Continuous'

dataset1 = Paul()
dataset1.subsample_genes(dataset1.nb_genes)
dataset2 = Tusi()
dataset2.time_traj = dataset2.time_traj[dataset2.batch_indices.ravel()==1]
dataset2.update_cells(dataset2.batch_indices.ravel()==1)
dataset2.subsample_genes(dataset2.nb_genes)
gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)

CompareModels(gene_dataset, dataset1, dataset2, plotname, 'writedata')

seurat1 = np.genfromtxt('../Seurat_data/' + plotname + '.1.CCA.txt')
seurat2 = np.genfromtxt('../Seurat_data/' + plotname + '.2.CCA.txt')
seurat, batch_indices, labels, keys, stats = run_model('readSeurat', gene_dataset, dataset1, dataset2,
                                                       filename=plotname)

dataset1, dataset2, gene_dataset = SubsetGenes(dataset1, dataset2, gene_dataset, plotname)
vae1, _, _, _, _ = run_model('vae', dataset1, 0, 0, filename=plotname, rep='vae1')
vae2, _, _, _, _ = run_model('vae', dataset2, 0, 0, filename=plotname, rep='vae2')
vae,  _, _, _, _ = run_model('vae', gene_dataset, dataset1, dataset2,filename=plotname, rep='0')
scanvi1,  _, _, _, _ = run_model('scanvi1', gene_dataset, dataset1, dataset2,filename=plotname, rep='0')
scanvi2,  _, _, _, _ = run_model('scanvi2', gene_dataset, dataset1, dataset2,filename=plotname, rep='0')

KNeighbors = np.concatenate([np.arange(10, 100, 10), np.arange(100, 500, 50)])
seurat_seurat = [KNNJaccardIndex(seurat1, seurat2, seurat, batch_indices, k)[0] for k in KNeighbors]
vae_seurat = [KNNJaccardIndex(vae1, vae2, seurat, batch_indices, k)[0] for k in KNeighbors]
vae_vae = [KNNJaccardIndex(vae1, vae2, vae, batch_indices, k)[0] for k in KNeighbors]
seurat_vae = [KNNJaccardIndex(seurat1, seurat2, vae, batch_indices, k)[0] for k in KNeighbors]
vae_scanvi1 = [KNNJaccardIndex(vae1, vae2, scanvi1, batch_indices, k)[0] for k in KNeighbors]
vae_scanvi2 = [KNNJaccardIndex(vae1, vae2, scanvi2, batch_indices, k)[0] for k in KNeighbors]

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 10))
plt.plot(KNeighbors, seurat_seurat,'r',label='Seurat_Seurat')
plt.plot(KNeighbors, vae_vae,'b',label='VAE_VAE')
plt.plot(KNeighbors, vae_scanvi1,'g',label='VAE_SCANVI1')
plt.plot(KNeighbors, vae_scanvi2,'y',label='VAE_SCANVI2')
# plt.plot(KNeighbors, vae_seurat,'g',label='VAE_Seurat')
# plt.plot(KNeighbors, seurat_vae,'y',label='Seurat_VAE')
legend = plt.legend(loc='lower right', shadow=False)
plt.savefig('../%s/%s.KNN.pdf' % (plotname,plotname))



def entropy_from_indices(indices):
    return entropy(np.array(np.unique(indices, return_counts=True)[1].astype(np.int32)))


def entropy_batch_mixing_subsampled(latent, batches, sample, n_neighbors=50, n_pools=50, n_samples_per_pool=100):
    X = latent[sample,:]
    nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(latent)
    indices = nbrs.kneighbors(X, return_distance=False)[:, 1:]
    batch_indices = np.vectorize(lambda i: batches[i])(indices)
    entropies = np.apply_along_axis(entropy_from_indices, axis=1, arr=batch_indices)
    if n_pools == 1:
        res = np.mean(entropies)
    else:
        res = np.mean([
            np.mean(entropies[np.random.choice(len(entropies), size=n_samples_per_pool)])
            for _ in range(n_pools)
        ])
    return res

import pandas as pd

def logit(p):
    p[p == 0] = np.min(p[p > 0])
    p[p == 1] = np.max(p[p < 1])
    return np.log(p / (1-p))

probs = np.exp(dataset2.time_traj)
fac, _ = pd.qcut(probs,20,retbins=True)
intervals,bin = np.unique(np.asarray(fac),return_inverse=True)
midint = [(x.left+x.right)/2 for x in intervals]

BE_vae=[]
for i in np.unique(bin):
    sample = np.arange(len(batch_indices))[batch_indices==1][bin==i]
    BE_vae.append(entropy_batch_mixing_subsampled(vae,batch_indices,sample))

BE_seurat=[]
for i in np.unique(bin):
    sample = np.arange(len(batch_indices))[batch_indices==1][bin==i]
    BE_seurat.append(entropy_batch_mixing_subsampled(seurat,batch_indices,sample))


BE_scanvi1=[]
for i in np.unique(bin):
    sample = np.arange(len(batch_indices))[batch_indices==1][bin==i]
    BE_scanvi1.append(entropy_batch_mixing_subsampled(scanvi1,batch_indices,sample))


BE_scanvi2=[]
for i in np.unique(bin):
    sample = np.arange(len(batch_indices))[batch_indices==1][bin==i]
    BE_scanvi2.append(entropy_batch_mixing_subsampled(scanvi2,batch_indices,sample))

import matplotlib.pyplot as plt
plt.figure(figsize=(10, 10))
plt.plot(np.arange(20), BE_vae,'r',label='BE VAE')
plt.plot(np.arange(20), BE_seurat,'b',label='BE Seurat')
plt.plot(np.arange(20), BE_scanvi1,'g',label='BE SCANVI1')
plt.plot(np.arange(20), BE_scanvi2,'y',label='BE SCANVI2')
legend = plt.legend(loc='lower right', shadow=False)
plt.savefig('../%s/%s.BE_potential.pdf' % (plotname,plotname))
