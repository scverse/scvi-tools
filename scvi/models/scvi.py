import numpy as np
import os
import torch
import pandas as pd

from typing import Optional, Union, List, Callable
from scvi.models.vae import VAE
from ._base import AbstractModelClass
from .de import DifferentialExpression
from scvi import _CONSTANTS
from scvi.inference import UnsupervisedTrainer
from scvi.inference import Posterior


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
        use_cuda=True,
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
        self.use_cuda = use_cuda and torch.cuda.is_available()

    # what args do we actually want here?
    def train(
        self,
        n_epochs=400,
        train_size=1.0,
        test_size=None,
        n_iter_kl_warmup=None,
        n_epochs_kl_warmup=400,
        **train_kwargs,
    ):
        self.model.train()

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

    @torch.no_grad()
    def get_z(self, adata=None, indices=None):
        self.model.eval()

        if self.is_trained is False:
            raise "Please train the model first."
        post = self._make_posterior(adata=adata, indices=indices, shuffle=False)
        latent = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            z = self.model.sample_from_posterior_z(x)
            latent += [z.cpu()]
        return np.array(torch.cat(latent))

        # for tensors in post:
        #     self.vae.guide(tensors)["z"]

    # assume posterior is a data loader + elbo + marginal ll
    def _make_posterior(self, adata=None, indices=None, shuffle=False):
        if adata is None:
            adata = self.adata
        if indices is None:
            indices = np.arange(adata.n_obs)
        return Posterior(
            self.model, adata, shuffle=shuffle, indices=indices, use_cuda=self.use_cuda
        )

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

    @torch.no_grad()
    def get_sample_scale(
        self,
        adata=None,
        indices=None,
        transform_batch: Optional[int] = None,
        gene_list: Optional[Union[np.ndarray, List[int]]] = None,
        library_size: float = 1,
        return_df: Optional[bool] = None,
        n_samples: int = 1,
        return_mean: bool = True,
    ) -> Union[np.ndarray, pd.DataFrame]:
        r"""Returns the frequencies of expression for the data.

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
        return_df
            Return a DataFrame instead of an `np.ndarray`. Includes gene
            names as columns. Requires either n_samples=1 or return_mean=True.
            When `gene_list` is not None and contains more than one gene, this is option is True.
            Otherwise, it defaults to False.
        n_samples
            Get sample scale from multiple samples.
        return_mean
            Whether to return the mean of the samples.

        Returns
        -------
        - **denoised_expression** - array of decoded expression adjusted for library size

        If ``n_samples`` > 1 and ``return_mean`` is False, then the shape is ``(samples, cells, genes)``.
        Otherwise, shape is ``(cells, genes)``. Return type is ``np.ndarray`` unless ``return_df`` is True.

        """
        post = self._make_posterior(adata=adata, indices=indices, shuffle=False)

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = self.gene_dataset.adata.var_names
            gene_mask = [True if gene in gene_list else False for gene in all_genes]
            if return_df is None and sum(gene_mask) > 1:
                return_df = True

        if n_samples > 1 and return_mean is False and return_df is True:
            raise ValueError(
                "return_df must be False if n_samples > 1 and return_mean is True"
            )

        px_scales = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            batch_idx = tensors[_CONSTANTS.BATCH_KEY]
            labels = tensors[_CONSTANTS.LABELS_KEY]
            px_scales += [
                np.array(
                    (
                        self.model.get_sample_scale(
                            x,
                            batch_index=batch_idx,
                            y=labels,
                            n_samples=n_samples,
                            transform_batch=transform_batch,
                        )[..., gene_mask]
                        * library_size
                    ).cpu()
                )
            ]

        if n_samples > 1:
            # The -2 axis correspond to cells.
            px_scales = np.concatenate(px_scales, axis=-2)
        else:
            px_scales = np.concatenate(px_scales, axis=0)

        if n_samples > 1 and return_mean:
            px_scales = px_scales.mean(0)

        if return_df is True:
            return pd.DataFrame(
                px_scales, columns=self.gene_dataset.gene_names[gene_mask]
            )
        else:
            return px_scales

    @torch.no_grad()
    def imputation(
        self,
        adata=None,
        indices=None,
        n_samples: Optional[int] = 1,
        transform_batch: Optional[Union[int, List[int]]] = None,
    ) -> np.ndarray:
        """Imputes px_rate over self cells

        Parameters
        ----------
        n_samples
            number of posterior samples
        transform_batch
            Batches to condition on.
            If transform_batch is:

            - None, then real observed batch is used
            - int, then batch transform_batch is used
            - list of int, then px_rates are averaged over provided batches.
        Returns
        -------
        type
            n_samples, n_cells, n_genes) px_rates squeezed array

        """
        post = self._make_posterior(adata=adata, indices=indices, shuffle=False)

        if (transform_batch is None) or (isinstance(transform_batch, int)):
            transform_batch = [transform_batch]
        imputed_arr = []
        for batch in transform_batch:
            imputed_list_batch = []
            for tensors in post:
                x = tensors[_CONSTANTS.X_KEY]
                batch_idx = tensors[_CONSTANTS.BATCH_KEY]
                labels = tensors[_CONSTANTS.LABELS_KEY]
                px_rate = self.model.get_sample_rate(
                    x,
                    batch_index=batch_idx,
                    y=labels,
                    n_samples=n_samples,
                    transform_batch=batch,
                )
                imputed_list_batch += [np.array(px_rate.cpu())]
            # axis 1 is cells if n_samples > 1
            imputed_arr.append(
                np.concatenate(imputed_list_batch, axis=1 if n_samples > 1 else 0)
            )
        imputed_arr = np.array(imputed_arr)
        # shape: (len(transformed_batch), n_samples, n_cells, n_genes) if n_samples > 1
        # else shape: (len(transformed_batch), n_cells, n_genes)
        return imputed_arr.mean(0).squeeze()

    def differential_expression(
        self, groupby, group1=None, group2="rest", adata=None, within_key=None
    ):
        # group 1 and group 2 are valid obs keys in the anndata
        # runs 1vsall or 1v1 based on group1 and group2 choices
        # runs within cluster
        # new differential expression class

        if group2 == "rest":
            idx2 = ~group1
        else:
            idx2 = 0

        DE = DifferentialExpression(self.get_sample_scale)
        pass

    # def create_posterior(
    #     self,
    #     model=None,
    #     adata: anndata.AnnData = None,
    #     shuffle=False,
    #     indices=None,
    #     type_class=Posterior,
    # ):

    #     model = self.model if model is None and hasattr(self, "model") else model
    #     adata = self.adata if adata is None and hasattr(self, "model") else adata
    #     return type_class(
    #         model,
    #         adata,
    #         shuffle=shuffle,
    #         indices=indices,
    #         use_cuda=self.use_cuda,
    #         batch_size=self.batch_size,
    #         data_loader_kwargs=self.data_loader_kwargs,
    #     )


# def posterior_predictive_sample(self):
#   # insert posterior predictive code generate function
#   pass

# def get_sample_scale(self, transform_batch: List[int]):
#   post = self._make_posterior()

#   for tensors in post:
#     # get sample scale code from current posterior
