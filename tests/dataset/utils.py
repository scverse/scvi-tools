import torch

from scvi.core.trainers.inference import UnsupervisedTrainer
from scvi.core.models.vae import VAE
from anndata import AnnData

use_cuda = torch.cuda.is_available()


def unsupervised_training_one_epoch(adata: AnnData):
    stats = adata.uns["_scvi"]["summary_stats"]
    vae = VAE(stats["n_genes"], stats["n_batch"], stats["n_labels"])
    trainer = UnsupervisedTrainer(vae, adata, train_size=0.5, use_cuda=use_cuda)
    trainer.train(n_epochs=1)
