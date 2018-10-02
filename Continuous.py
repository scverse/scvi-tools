from scvi.dataset.paul_tusi import Paul,Tusi
from scvi.dataset.dataset import GeneExpressionDataset
import numpy as np
use_cuda = True
from scvi.harmonization.utils_chenling import run_model,eval_latent
from copy import deepcopy

ngenes=1000
model_type = 'vae'
plotname = 'Continuous'

dataset1 = Paul()
dataset1.subsample_genes(dataset1.nb_genes)
dataset2 = Tusi()
dataset2.subsample_genes(dataset2.nb_genes)
dataset2.update_cells(dataset2.batch_indices.ravel()==1)

gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)
gene_dataset.subsample_genes(500)
latent, batch_indices, labels, keys, stats = run_model(model_type, gene_dataset, dataset1, dataset2,
                                                       filename=plotname)

dataset1 = deepcopy(gene_dataset)
dataset1.update_cells(gene_dataset.batch_indices.ravel() == 0)
dataset1.subsample_genes(dataset1.nb_genes)
latent1, _, _, _, _ = run_model('vae', dataset1, 0, 0, filename=plotname)
dataset2 = deepcopy(gene_dataset)
dataset2.update_cells(gene_dataset.batch_indices.ravel() == 1)
dataset2.subsample_genes(dataset2.nb_genes)
latent2, _, _, _, _ = run_model('vae', dataset2, 0, 0, filename=plotname)

np.save(file='../Paul.latent.npy',arr=latent1.astype('float64'))
np.save(file='../Tusi.latent.npy',arr=latent2.astype('float64'))
np.save(file='../merged.latent.npy',arr=latent.astype('float64'))

from scvi.models.vae import VAE
from scvi.models.scanvi import SCANVI
from scvi.inference import UnsupervisedTrainer, SemiSupervisedTrainer
import torch

vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches, n_labels=gene_dataset.n_labels,
          n_hidden=128, n_latent=10, n_layers=2, dispersion='gene')
trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=1.0)
trainer.train(n_epochs=500)

scanvi = SCANVI(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels, n_layers=2)
scanvi.load_state_dict(vae.state_dict(), strict=False)
trainer_scanvi = SemiSupervisedTrainer(scanvi, gene_dataset, classification_ratio=1,
                                       n_epochs_classifier=1, lr_classification=5 * 1e-3)
trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices.ravel() ==0))
trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(
    indices=(gene_dataset.batch_indices.ravel() == 1)
)
trainer_scanvi.train(n_epochs=50)

torch.save(trainer.model,'cont.vae.model.pkl')
torch.save(trainer_scanvi.model,'cont.scanvi.model.pkl')


