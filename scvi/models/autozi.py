import numpy as np
import logging
import torch
from anndata import AnnData

from typing import Optional, Dict, Union
from scvi._compat import Literal
from scvi.models._modules.autozivae import AutoZIVAE
from scvi.models import SCVI

from scvi.inference.posterior import Posterior

logger = logging.getLogger(__name__)


class AUTOZI(SCVI):
    """Automatic identification of ZI [Clivio19]_

    Parameters
    ----------
    adata
        AnnData object that has been registered with scvi
    n_hidden
        Number of nodes per hidden layer
    n_latent
        Dimensionality of the latent space
    n_layers
        Number of hidden layers used for encoder NN
    dropout_rate
        Dropout rate for neural networks
    dispersion
        One of the following

        * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - dispersion can differ between different batches
        * ``'gene-label'`` - dispersion can differ between different labels
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
    latent_distribution
        One of

        * ``'normal'`` - Normal distribution
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)
    alpha_prior
        Float denoting the alpha parameter of the prior Beta distribution of
        the zero-inflation Bernoulli parameter. Should be between 0 and 1, not included.
        When set to ``None'', will be set to 1 - beta_prior if beta_prior is not ``None'',
        otherwise the prior Beta distribution will be learned on an Empirical Bayes fashion.
    beta_prior
        Float denoting the beta parameter of the prior Beta distribution of
        the zero-inflation Bernoulli parameter. Should be between 0 and 1, not included.
        When set to ``None'', will be set to 1 - alpha_prior if alpha_prior is not ``None'',
        otherwise the prior Beta distribution will be learned on an Empirical Bayes fashion.
    minimal_dropout
        Float denoting the lower bound of the cell-gene ZI rate in the ZINB component.
        Must be non-negative. Can be set to 0 but not recommended as this may make
        the mixture problem ill-defined.
    zero_inflation: One of the following

        * ``'gene'`` - zero-inflation Bernoulli parameter of AutoZI is constant per gene across cells
        * ``'gene-batch'`` - zero-inflation Bernoulli parameter can differ between different batches
        * ``'gene-label'`` - zero-inflation Bernoulli parameter can differ between different labels
        * ``'gene-cell'`` - zero-inflation Bernoulli parameter can differ for every gene in every cell


    Examples
    --------

    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.dataset.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.models.AutoZIVAE(adata)
    >>> vae.train(n_epochs=400)
    """

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        latent_distribution: Literal["normal", "ln"] = "normal",
        use_cuda: bool = True,
        alpha_prior: Optional[float] = 0.5,
        beta_prior: Optional[float] = 0.5,
        minimal_dropout: float = 0.01,
        zero_inflation: str = "gene",
        **model_kwargs,
    ):
        assert (
            "scvi_data_registry" in adata.uns.keys()
        ), "Please setup your AnnData with scvi.dataset.setup_anndata(adata) first"

        self.adata = adata
        summary_stats = adata.uns["scvi_summary_stats"]
        self.model = AutoZIVAE(
            n_input=summary_stats["n_genes"],
            n_batch=summary_stats["n_batch"],
            n_labels=summary_stats["n_labels"],
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers_encoder=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            latent_distribution=latent_distribution,
            **model_kwargs,
        )
        self.is_trained = False
        self.use_cuda = use_cuda and torch.cuda.is_available()
        self.batch_size = 128
        self._posterior_class = Posterior

    def get_alphas_betas(
        self, as_numpy: bool = True
    ) -> Dict[str, Union[torch.Tensor, np.ndarray]]:
        """Return parameters of Bernoulli Beta distributions in a dictionary"""

        return self.model.get_alphas_beta(as_numpy=as_numpy)
