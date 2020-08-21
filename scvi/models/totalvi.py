import logging
from anndata import AnnData

from scvi._compat import Literal
from scvi.core.models import TOTALVAE

from scvi.models._base import BaseModelClass, CiteSeqMixin
from scvi.core.posteriors import TotalPosterior
from scvi.core.trainers import TotalTrainer

logger = logging.getLogger(__name__)


class TOTALVI(CiteSeqMixin, BaseModelClass):
    """total Variational Inference [GayosoSteier20]_

    Parameters
    ----------
    adata
        AnnData object that has been registered with scvi
    n_latent
        Dimensionality of the latent space
    gene_dispersion
        One of the following

        * ``'gene'`` - genes_dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - genes_dispersion can differ between different batches
        * ``'gene-label'`` - genes_dispersion can differ between different labels
    protein_dispersion
        One of the following

        * ``'protein'`` - protein_dispersion parameter is constant per protein across cells
        * ``'protein-batch'`` - protein_dispersion can differ between different batches NOT TESTED
        * ``'protein-label'`` - protein_dispersion can differ between different labels NOT TESTED
    gene_likelihood
        One of

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
    latent_distribution
        One of

        * ``'normal'`` - Normal distribution
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)

    Examples
    --------

    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.dataset.setup_anndata(adata, batch_key="batch", protein_expression_obsm_key="protein_expression")
    >>> vae = scvi.models.TOTALVI(adata)
    >>> vae.train(n_epochs=400)
    >>> adata.obsm["X_totalVI"] = vae.get_latent_representation()
    """

    def __init__(
        self,
        adata: AnnData,
        n_latent: int = 20,
        gene_dispersion: Literal[
            "gene", "gene-batch", "gene-label", "gene-cell"
        ] = "gene",
        protein_dispersion: Literal[
            "protein", "protein-batch", "protein-label"
        ] = "protein",
        gene_likelihood: Literal["zinb", "nb"] = "nb",
        latent_distribution: Literal["normal", "ln"] = "ln",
        use_cuda: bool = True,
        **model_kwargs,
    ):
        super(TOTALVI, self).__init__(adata, use_cuda=use_cuda)
        if "totalvi_batch_mask" in adata.uns.keys():
            batch_mask = adata.uns["totalvi_batch_mask"]
        else:
            batch_mask = None
        self.model = TOTALVAE(
            n_input_genes=self.summary_stats["n_genes"],
            n_input_proteins=self.summary_stats["n_proteins"],
            n_batch=self.summary_stats["n_batch"],
            gene_dispersion=gene_dispersion,
            protein_dispersion=protein_dispersion,
            reconstruction_loss_gene=gene_likelihood,
            latent_distribution=latent_distribution,
            protein_batch_mask=batch_mask,
            **model_kwargs,
        )
        self._model_summary_string = (
            "TotalVI Model with following params: \nn_latent: {}, "
            "gene_dispersion: {}, protein_dispersion: {}, gene_likelihood: {}, latent_distribution: {}"
        ).format(
            n_latent,
            gene_dispersion,
            protein_dispersion,
            gene_likelihood,
            latent_distribution,
        )
        self._posterior_class = TotalPosterior
        self._trainer_class = TotalTrainer

    def train(
        self,
        n_epochs=400,
        train_size=0.9,
        test_size=None,
        lr=1e-3,
        n_iter_kl_warmup=None,
        n_epochs_kl_warmup=400,
        batch_size=256,
        metric_frequency=1,
        trainer_kwargs={},
        train_kwargs={},
    ):

        if "totalvi_batch_mask" in self.adata.uns.keys():
            imputation = True
        else:
            imputation = False
        self.trainer = TotalTrainer(
            self.model,
            self.adata,
            train_size=train_size,
            test_size=test_size,
            n_iter_kl_warmup=n_iter_kl_warmup,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
            frequency=metric_frequency,
            batch_size=batch_size,
            use_adversarial_loss=imputation,
            use_cuda=self.use_cuda,
            **trainer_kwargs,
        )
        self.trainer.train(n_epochs=n_epochs, lr=lr, **train_kwargs)
        self.is_trained = True
        self.train_indices = self.trainer.train_set.indices
        self.test_indices = self.trainer.test_set.indices
        self.validation_indices = self.trainer.validation_set.indices
