import numpy as np
import os
import logging
import torch
import pandas as pd

from typing import Optional, Union, List
from scvi._compat import Literal
from scvi.models._modules.vae import VAE
from scvi.models._base import AbstractModelClass

# from scvi.models._differential import DifferentialExpression
from scvi.models._modules.distributions import (
    NegativeBinomial,
    ZeroInflatedNegativeBinomial,
)
from scvi import _CONSTANTS
from scvi.inference.inference import UnsupervisedTrainer
from scvi.inference.posterior import Posterior

logger = logging.getLogger(__name__)


class SCVI(AbstractModelClass):
    def __init__(
        self,
        adata,
        n_hidden=128,
        n_latent=10,
        n_layers=1,
        dropout_rate=0.1,
        dispersion="gene",
        reconstruction_loss="zinb",
        latent_distribution="normal",
        use_cuda=True,
        **model_kwargs,
    ):
        assert (
            "scvi_data_registry" in adata.uns.keys()
        ), "Please setup your AnnData with setup_anndata() first"

        self.adata = adata
        summary_stats = adata.uns["scvi_summary_stats"]
        self.model = VAE(
            n_input=summary_stats["n_genes"],
            n_batch=summary_stats["n_batch"],
            n_labels=summary_stats["n_labels"],
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            reconstruction_loss=reconstruction_loss,
            latent_distribution=latent_distribution,
            **model_kwargs,
        )
        self.is_trained = False
        self.use_cuda = use_cuda and torch.cuda.is_available()

    # assume posterior is a data loader + elbo + marginal ll
    def _make_posterior(self, adata=None, indices=None, batch_size=128):
        if adata is None:
            adata = self.adata
        if indices is None:
            indices = np.arange(adata.n_obs)
        post = Posterior(
            self.model,
            adata,
            shuffle=False,
            indices=indices,
            use_cuda=self.use_cuda,
            batch_size=batch_size,
        ).sequential()
        return post

    def save(self, dir_path):
        # save the model state dict and the trainer state dict only
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        else:
            raise ValueError(
                "{} already exists. Please provide an unexisting directory for saving.".format(
                    dir_path
                )
            )
        torch.save(self.model.state_dict(), os.path.join(dir_path, "model_params.pt"))
        torch.save(
            self.trainer.optimizer.state_dict(),
            os.path.join(dir_path, "optimizer_params.pt"),
        )

    def load(self, dir_path):
        # load state dicts, maybe a class method?
        # Loading scVI model
        model_path = os.path.join(dir_path, "model_params.pt")
        optimizer_path = os.path.join(dir_path, "optimizer_params.pt")
        if self.use_cuda:
            self.model.load_state_dict(torch.load(model_path))
            self.trainer.optimizer.load_state_dict(torch.load(optimizer_path))
            self.model.cuda()
        else:
            device = torch.device("cpu")
            self.model.load_state_dict(torch.load(model_path, map_location=device))
            self.trainer.optimizer.load_state_dict(
                torch.load(optimizer_path, map_location=device)
            )
        self.model.eval()

    # what args do we actually want here?
    def train(
        self,
        n_epochs=400,
        train_size=0.9,
        test_size=None,
        learning_rate=1e-3,
        n_iter_kl_warmup=None,
        n_epochs_kl_warmup=400,
        trainer_kwargs={},
        train_kwargs={},
    ):

        self.trainer = UnsupervisedTrainer(
            self.model,
            self.adata,
            train_size=train_size,
            test_size=test_size,
            n_iter_kl_warmup=n_iter_kl_warmup,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
            **trainer_kwargs,
        )
        self.trainer.train(n_epochs=n_epochs, lr=learning_rate, **train_kwargs)
        self.is_trained = True

    @torch.no_grad()
    def get_latent_representation(
        self, adata=None, indices=None, give_mean=True, mc_samples=5000
    ):
        if self.is_trained is False:
            raise RuntimeError("Please train the model first.")

        post = self._make_posterior(adata=adata, indices=indices, shuffle=False)
        latent = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            z = self.model.sample_from_posterior_z(
                x, give_mean=give_mean, n_samples=mc_samples
            )
            latent += [z.cpu()]
        return np.array(torch.cat(latent))

    @torch.no_grad()
    def get_latent_library_size(self, adata=None, indices=None, give_mean=True):
        if self.is_trained is False:
            raise RuntimeError("Please train the model first.")

        post = self._make_posterior(adata=adata, indices=indices, shuffle=False)
        libraries = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            out = self.model.get_latents(x)
            if give_mean is False:
                library = out["library"]
            else:
                library = torch.distributions.LogNormal(
                    out["ql_m"], out["ql_v"].sqrt()
                ).mean
            libraries += [library.cpu()]
        return np.array(torch.cat(libraries))

    @torch.no_grad()
    def get_normalized_expression(
        self,
        adata=None,
        indices=None,
        transform_batch: Optional[int] = None,
        gene_list: Optional[Union[np.ndarray, List[int]]] = None,
        library_size: Optional[Union[float, Literal["observed"]]] = 1,
        n_samples: int = 1,
        return_mean: bool = True,
        return_numpy: Optional[bool] = None,
    ) -> Union[np.ndarray, pd.DataFrame]:
        r"""Returns the normalized (decoded) gene expression.

        This is denoted as :math:`\rho_n` in the scVI paper.

        Parameters
        ----------
        transform_batch
            Batch to condition on.
            If transform_batch is:

            - None, then real observed batch is used
            - int, then batch transform_batch is used
        gene_list
            Return frequencies of expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
        library_size
            Scale the expression frequencies to a common library size.
            This allows gene expression levels to be interpreted on a common scale of relevant
            magnitude.
        n_samples
            Get sample scale from multiple samples.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a `np.ndarray` instead of a `pd.DataFrame`. Includes gene
            names as columns. If either n_samples=1 or return_mean=True, defaults to False.
            Otherwise, it defaults to True.

        Returns
        -------
        - **normalized_expression** - array of normalized expression

        If ``n_samples`` > 1 and ``return_mean`` is False, then the shape is ``(samples, cells, genes)``.
        Otherwise, shape is ``(cells, genes)``. Return type is ``pd.DataFrame`` unless ``return_numpy`` is True.

        """
        post = self._make_posterior(adata=adata, indices=indices, shuffle=False)

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = adata.var_names
            gene_mask = [True if gene in gene_list else False for gene in all_genes]

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                logger.warning(
                    "return_numpy must be True if n_samples > 1 and return_mean is False, returning np.ndarray"
                )
            return_numpy = True

        if library_size == "observed":
            model_fn = self.model.get_sample_rate
            scaling = 1
        else:
            model_fn = self.model.get_sample_scale
            scaling = library_size

        exprs = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            batch_idx = tensors[_CONSTANTS.BATCH_KEY]
            labels = tensors[_CONSTANTS.LABELS_KEY]
            exprs += [
                np.array(
                    (
                        model_fn(
                            x,
                            batch_index=batch_idx,
                            y=labels,
                            n_samples=n_samples,
                            transform_batch=transform_batch,
                        )[..., gene_mask]
                        * scaling
                    ).cpu()
                )
            ]

        if n_samples > 1:
            # The -2 axis correspond to cells.
            exprs = np.concatenate(exprs, axis=-2)
        else:
            exprs = np.concatenate(exprs, axis=0)

        if n_samples > 1 and return_mean:
            exprs = exprs.mean(0)

        if return_numpy is None or return_numpy is False:
            return pd.DataFrame(exprs, columns=adata.gene_names[gene_mask])
        else:
            return exprs

    def differential_expression(
        self, groupby, group1=None, group2="rest", adata=None, within_key=None
    ):
        # group 1 and group 2 are valid obs keys in the anndata
        # runs 1vsall or 1v1 based on group1 and group2 choices
        # runs within cluster
        # new differential expression class

        raise NotImplementedError

        # if group2 == "rest":
        #     idx2 = ~group1
        # else:
        #     idx2 = 0

        # DE = DifferentialExpression(self.get_sample_scale)
        # pass

    @torch.no_grad()
    def posterior_predictive_sample(
        self,
        adata=None,
        indices=None,
        n_samples: int = 1,
        gene_list: Union[list, np.ndarray] = None,
    ) -> np.ndarray:
        r"""Generate observation samples from the posterior predictive distribution

        The posterior predictive distribution is written as :math:`p(\hat{x} \mid x)`.

        Parameters
        ----------
        n_samples
            Number of required samples for each cell
        gene_list
            Indices or names of genes of interest

        Returns
        -------
        x_new : :py:class:`torch.Tensor`
            tensor with shape (n_cells, n_genes, n_samples)
        """
        assert self.model.reconstruction_loss in ["zinb", "nb", "poisson"]

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = adata.var_names
            gene_mask = [True if gene in gene_list else False for gene in all_genes]

        post = self._make_posterior(adata=adata, indices=indices)

        x_new = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            batch_idx = tensors[_CONSTANTS.BATCH_KEY]
            labels = tensors[_CONSTANTS.LABELS_KEY]
            outputs = self.model.inference(
                x, batch_index=batch_idx, y=labels, n_samples=n_samples
            )
            px_r = outputs["px_r"]
            px_rate = outputs["px_rate"]
            px_dropout = outputs["px_dropout"]

            if self.model.reconstruction_loss == "poisson":
                l_train = px_rate
                l_train = torch.clamp(l_train, max=1e8)
                dist = torch.distributions.Poisson(
                    l_train
                )  # Shape : (n_samples, n_cells_batch, n_genes)
            elif self.model.reconstruction_loss == "nb":
                dist = NegativeBinomial(mu=px_rate, theta=px_r)
            elif self.model.reconstruction_loss == "zinb":
                dist = ZeroInflatedNegativeBinomial(
                    mu=px_rate, theta=px_r, zi_logits=px_dropout
                )
            else:
                raise ValueError(
                    "{} reconstruction error not handled right now".format(
                        self.model.reconstruction_loss
                    )
                )
            exprs = dist.sample().permute(
                [1, 2, 0]
            )  # Shape : (n_cells_batch, n_genes, n_samples)

            if gene_list is not None:
                x_new = x_new[:, gene_mask, :]

            x_new.append(exprs.cpu())

        x_new = torch.cat(x_new)  # Shape (n_cells, n_genes, n_samples)

        return x_new.numpy()
