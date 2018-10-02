from scvi.dataset.paul_tusi import Paul,Tusi
from scvi.dataset.dataset import GeneExpressionDataset
import numpy as np
use_cuda = True
from scvi.harmonization.utils_chenling import run_model
from copy import deepcopy

model_type = 'vae'
plotname = 'Continuous'

dataset1 = Paul()
dataset1.subsample_genes(dataset1.nb_genes)
dataset2 = Tusi()
dataset2.update_cells(dataset2.batch_indices.ravel()==1)
dataset2.subsample_genes(dataset2.nb_genes)

gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)
gene_dataset.subsample_genes(500)

dataset1 = deepcopy(gene_dataset)
dataset1.update_cells(gene_dataset.batch_indices.ravel() == 0)
dataset1.subsample_genes(dataset1.nb_genes)
latent1, _, _, _, _ = run_model('vae', dataset1, 0, 0, filename=plotname,rep='Paul')
np.save(file='../Continuous/Paul.vae.latent.npy',arr=latent1.astype('float64'))
dataset2 = deepcopy(gene_dataset)
dataset2.update_cells(gene_dataset.batch_indices.ravel() == 1)
dataset2.subsample_genes(dataset2.nb_genes)
latent2, _, _, _, _ = run_model('vae', dataset2, 0, 0, filename=plotname,rep='Tusi')
np.save(file='../Continuous/Tusi.vae.latent.npy',arr=latent2.astype('float64'))

_, _, _, _, _ = run_model('writedata', gene_dataset, 0, 0, filename=plotname,rep='')

latent, batch_indices, labels, keys, stats = run_model('vae', gene_dataset, 0, 0, filename=plotname,rep='merged')
np.save(file='../Continuous/merged.vae.latent.npy',arr=latent.astype('float64'))


latent1, _, _, _, _ = run_model('scanvi', dataset1, 0, 0, filename=plotname,rep='Paul')
latent2, _, _, _, _ = run_model('scanvi', dataset2, 0, 0, filename=plotname,rep='Tusi')
np.save(file='../Continuous/Paul.scanvi.latent.npy',arr=latent1.astype('float64'))
np.save(file='../Continuous/Tusi.scanvi.latent.npy',arr=latent2.astype('float64'))

latent, batch_indices, labels, keys, stats = run_model('scanvi1', gene_dataset, 0, 0, filename=plotname,rep='merged')
np.save(file='../Continuous/merged.scanvi1.latent.npy',arr=latent.astype('float64'))

