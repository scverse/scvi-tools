import numpy as np

from scvi.models import VAE
from scvi.inference import UnsupervisedTrainer
from scvi.inference import Posterior
from .base import AbstractModelClass


class SCVI(AbstractModelClass):
    def __init__(
        self,
        adata,
        n_batch=0,
        n_labels=0,
        n_hidden=128,
        n_latent=10,
        n_layers=1,
        dropout_rate=0.1,
        dispersion="gene",
        log_variational=True,
        reconstruction_loss="zinb",
        latent_distribution="normal",
    ):
        assert (
            "scvi_data_registry" in adata.uns.keys()
        ), "Please setup your AnnData with setup_anndata() first"

        self.adata = adata
        summary_stats = adata.uns["scvi_summary_stats"]
        self.model = VAE(
            n_input=summary_stats["n_genes"],
            n_batch=n_batch,
            n_labels=n_labels,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            log_variational=log_variational,
            reconstruction_loss=reconstruction_loss,
            latent_distribution=latent_distribution,
        )
        self.is_trained = False

    # what args do we actually want here?
    def train(
        self,
        n_epochs=400,
        train_size=1.0,
        test_size=None,
        n_iter_kl_warmup=None,
        n_epochs_kl_warmup=400,
        **train_kwargs
    ):
        self.trainer = UnsupervisedTrainer(
            self.model,
            self.adata,
            train_size=train_size,
            test_size=test_size,
            n_iter_kl_warmup=n_iter_kl_warmup,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
        )
        self.trainer.train(n_epochs=n_epochs)
        self.is_trained = True

    def get_z(self, adata=None, indices=None):
        if self.is_trained is False:
            raise "Please train the model first."

        post = self._make_posterior()
        for tensors in post:
            self.vae.guide(tensors)["z"]
            for tensors in post:
                self.vae.sample_from_posterior_z()

    def _make_posterior(self, adata=None, indices=None):
        if adata is None:
            adata = self.adata
        if indices is None:
            indices = np.arange(adata.n_obs)
        return Posterior(adata)

    # def save(self, file_name):
    #   # save the model state dict and the trainer state dict only
    #   pass

    # def load(self):
    #   # load state dicts, maybe a class method?
    #   pass


# def _make_posterior(self, adata = None):
# 	if adata = None:
#     adata = self.adata
#   return post = Posterior(adata)


# def differential_expression(self, adata = None, group1 = None, group2 = "rest", within_key=None):
#   # group 1 and group 2 are valid obs keys in the anndata
#   # runs 1vsall or 1v1 based on group1 and group2 choices
#   # runs within cluster
#   # new differential expression class
#   pass


# def posterior_predictive_sample(self):
#   # insert posterior predictive code generate function
#   pass

# def get_sample_scale(self, transform_batch: List[int]):
#   post = self._make_posterior()

#   for tensors in post:
#     # get sample scale code from current posterior
