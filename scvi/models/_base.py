import torch
import os
import logging
import pandas as pd
import numpy as np

from anndata import AnnData
from scvi.dataset._anndata import get_from_registry
from collections.abc import Iterable
from functools import partial
from scvi.core._distributions import NegativeBinomial, ZeroInflatedNegativeBinomial
from scvi.models._differential import DifferentialComputation
from scvi.models._utils import scrna_raw_counts_properties
from scvi import _CONSTANTS
from typing import Optional, Union, List, Dict, Tuple
from scvi._compat import Literal
from scvi.core.trainers import UnsupervisedTrainer
from abc import ABC, abstractmethod

logger = logging.getLogger(__name__)


class RNASeqMixin:
    def train(
        self,
        n_epochs=400,
        train_size=0.9,
        test_size=None,
        lr=1e-3,
        n_iter_kl_warmup=None,
        n_epochs_kl_warmup=400,
        metric_frequency=None,
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
            frequency=metric_frequency,
            use_cuda=self.use_cuda,
            **trainer_kwargs,
        )
        self.trainer.train(n_epochs=n_epochs, lr=lr, **train_kwargs)
        self.is_trained = True
        self.train_indices = self.trainer.train_set.indices
        self.test_indices = self.trainer.test_set.indices
        self.validation_indices = self.trainer.validation_set.indices

    @torch.no_grad()
    def get_elbo(self, adata=None, indices=None):

        post = self._make_posterior(adata=adata, indices=indices)

        return -post.elbo()

    @torch.no_grad()
    def get_marginal_ll(self, adata=None, indices=None, n_mc_samples=1000):

        post = self._make_posterior(adata=adata, indices=indices)

        return -post.marginal_ll(n_mc_samples=n_mc_samples)

    @torch.no_grad()
    def get_reconstruction_error(self, adata=None, indices=None):

        post = self._make_posterior(adata=adata, indices=indices)

        return -post.reconstruction_error()

    @torch.no_grad()
    def get_normalized_expression(
        self,
        adata=None,
        indices=None,
        transform_batch: Optional[
            int
        ] = None,  # should take categories and ints (de requires ints )
        gene_list: Optional[Union[np.ndarray, List[int]]] = None,
        library_size: Optional[Union[float, Literal["latent"]]] = 1,
        n_samples: int = 1,
        return_mean: bool = True,
        return_numpy: Optional[bool] = None,
    ) -> Union[np.ndarray, pd.DataFrame]:
        # whatever is here will pass it into cat and batch
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

        adata = adata if adata is not None else self.adata
        post = self._make_posterior(adata=adata, indices=indices)

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

        if library_size == "latent":
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
            return pd.DataFrame(
                exprs, columns=adata.var_names[gene_mask], index=adata.obs_names
            )
        else:
            return exprs

    def differential_expression(
        self,
        groupby,
        group1=None,
        group2=None,
        adata=None,
        mode="vanilla",
        within_key=None,
        all_stats=True,
    ):
        if adata is None:
            adata = self.adata
        cell_idx1 = adata.obs[groupby] == group1
        if group2 is None:
            cell_idx2 = ~cell_idx1
        else:
            cell_idx2 = adata.obs[groupby] == group2
        dc = DifferentialComputation(self.get_normalized_expression, adata)
        all_info = dc.get_bayes_factors(cell_idx1, cell_idx2)

        gene_names = adata.var_names
        if all_stats is True:
            (
                mean1,
                mean2,
                nonz1,
                nonz2,
                norm_mean1,
                norm_mean2,
            ) = scrna_raw_counts_properties(cell_idx1, cell_idx2)
            genes_properties_dict = dict(
                raw_mean1=mean1,
                raw_mean2=mean2,
                non_zeros_proportion1=nonz1,
                non_zeros_proportion2=nonz2,
                raw_normalized_mean1=norm_mean1,
                raw_normalized_mean2=norm_mean2,
            )
            all_info = {**all_info, **genes_properties_dict}

        res = pd.DataFrame(all_info, index=gene_names)
        sort_key = "proba_de" if mode == "change" else "bayes_factor"
        res = res.sort_values(by=sort_key, ascending=False)
        return res

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

        adata = adata if adata is not None else self.adata
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
                exprs = exprs[:, gene_mask, :]

            x_new.append(exprs.cpu())

        x_new = torch.cat(x_new)  # Shape (n_cells, n_genes, n_samples)

        return x_new.numpy()

    @torch.no_grad()
    def _get_denoised_samples(
        self,
        adata=None,
        indices=None,
        n_samples: int = 25,
        batch_size: int = 64,
        rna_size_factor: int = 1000,
        transform_batch: Optional[int] = None,
    ) -> np.ndarray:
        """Return samples from an adjusted posterior predictive.

        Parameters
        ----------
        n_samples
            How may samples per cell
        batch_size
            Mini-batch size for sampling. Lower means less GPU memory footprint
        rna_size_factor
            size factor for RNA prior to sampling gamma distribution
        transform_batch
            int of which batch to condition on for all cells

        Returns
        -------

        """

        post = self._make_posterior(adata=adata, indices=indices)

        posterior_list = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            batch_idx = tensors[_CONSTANTS.BATCH_KEY]
            labels = tensors[_CONSTANTS.LABELS_KEY]

            outputs = self.model.inference(
                x, batch_index=batch_idx, y=labels, n_samples=n_samples
            )
            px_scale = outputs["px_scale"]
            px_r = outputs["px_r"]

            rate = rna_size_factor * px_scale
            if len(px_r.size()) == 2:
                px_dispersion = px_r
            else:
                px_dispersion = torch.ones_like(x) * px_r

            # This gamma is really l*w using scVI manuscript notation
            p = rate / (rate + px_dispersion)
            r = px_dispersion
            l_train = torch.distributions.Gamma(r, (1 - p) / p).sample()
            data = l_train.cpu().numpy()
            # """
            # In numpy (shape, scale) => (concentration, rate), with scale = p /(1 - p)
            # rate = (1 - p) / p  # = 1/scale # used in pytorch
            # """
            posterior_list += [data]

            posterior_list[-1] = np.transpose(posterior_list[-1], (1, 2, 0))

        return np.concatenate(posterior_list, axis=0)

    @torch.no_grad()
    def get_feature_correlation_matrix(
        self,
        adata=None,
        indices=None,
        n_samples: int = 10,
        batch_size: int = 64,
        rna_size_factor: int = 1000,
        transform_batch: Optional[Union[int, List[int]]] = None,
        correlation_type: Literal["spearman", "pearson"] = "spearman",
    ) -> pd.DataFrame:
        """Generate gene-gene correlation matrix using scvi uncertainty and expression

        Parameters
        ----------
        n_samples
            How may samples per cell
        batch_size
            Mini-batch size for sampling. Lower means less GPU memory footprint
        rna_size_factor
            size factor for RNA prior to sampling gamma distribution
        transform_batch
            Batches to condition on.
            If transform_batch is:

            - None, then real observed batch is used
            - int, then batch transform_batch is used
            - list of int, then values are averaged over provided batches.
        correlation_type
            One of "pearson", "spearman"

        Returns
        -------
        Gene-gene correlation matrix
        """

        from scipy.stats import spearmanr

        if (transform_batch is None) or (isinstance(transform_batch, int)):
            transform_batch = [transform_batch]
        corr_mats = []
        for b in transform_batch:
            denoised_data = self._get_denoised_samples(
                adata=adata,
                indicies=indices,
                n_samples=n_samples,
                batch_size=batch_size,
                rna_size_factor=rna_size_factor,
                transform_batch=b,
            )
            flattened = np.zeros(
                (denoised_data.shape[0] * n_samples, denoised_data.shape[1])
            )
            for i in range(n_samples):
                flattened[
                    denoised_data.shape[0] * (i) : denoised_data.shape[0] * (i + 1)
                ] = denoised_data[:, :, i]
            if correlation_type == "pearson":
                corr_matrix = np.corrcoef(flattened, rowvar=False)
            elif correlation_type == "spearman":
                corr_matrix, _ = spearmanr(flattened)
            else:
                raise ValueError(
                    "Unknown correlation type. Choose one of 'spearman', 'pearson'."
                )
            corr_mats.append(corr_matrix)
        corr_matrix = np.mean(np.stack(corr_mats), axis=0)
        return pd.DataFrame(
            corr_matrix, index=self.adata.var_names, columns=self.adata.var_names
        )

    @torch.no_grad()
    def get_likelihood_parameters(
        self,
        adata=None,
        indices=None,
        n_samples: Optional[int] = 1,
        give_mean: Optional[bool] = False,
    ) -> Dict[str, np.ndarray]:

        r"""Estimates for the parameters of the likelihood :math:`p(x \mid z)`.


        """

        post = self._make_posterior(adata=adata, indices=indices)

        dropout_list = []
        mean_list = []
        dispersion_list = []
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

            n_batch = px_rate.size(0) if n_samples == 1 else px_rate.size(1)
            dispersion_list += [
                np.repeat(np.array(px_r.cpu())[np.newaxis, :], n_batch, axis=0)
            ]
            mean_list += [np.array(px_rate.cpu())]
            dropout_list += [np.array(px_dropout.cpu())]

        dropout = np.concatenate(dropout_list)
        means = np.concatenate(mean_list)
        dispersions = np.concatenate(dispersion_list)
        if give_mean and n_samples > 1:
            dropout = dropout.mean(0)
            means = means.mean(0)

        return_dict = {}
        return_dict["mean"] = means

        if self.model.reconstruction_loss == "zinb":
            return_dict["dropout"] = dropout
            return_dict["dispersions"] = dispersions
        if self.model.reconstruction_loss == "nb":
            return_dict["dispersions"] = dispersions

        return return_dict

    @torch.no_grad()
    def get_latent_representation(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Union[np.ndarray, List[int]]] = None,
        give_mean: bool = True,
        mc_samples: int = 5000,
    ) -> np.ndarray:
        """Return the latent representation for each cell

        Parameters
        ----------
        adata
            AnnData object that has been registered with scvi. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        give_mean
            Give mean of distribution or sample from it
        mc_samples
            For distributions with no closed-form mean (e.g., `logistic normal`), how many Monte Carlo
            samples to take for computing mean.

        Returns
        -------
        latent_representation : np.ndarray
            Low-dimensional representation for each cell

        Examples
        --------

        >>> vae = scvi.model.SCVI(adata)
        >>> vae.train(n_epochs=400)
        >>> adata.obsm["X_scVI"] = vae.get_latent_representation()

        We can also get the latent representation for a subset of cells

        >>> adata_subset = adata[adata.obs.cell_type == "really cool cell type"]
        >>> latent_subset = vae.get_latent_representation(adata_subset)

        """

        if self.is_trained is False:
            raise RuntimeError("Please train the model first.")

        post = self._make_posterior(adata=adata, indices=indices)
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

        post = self._make_posterior(adata=adata, indices=indices)
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


class CiteSeqMixin(RNASeqMixin):
    @torch.no_grad()
    def get_reconstruction_error(self, adata=None, indices=None, mode="total"):

        post = self._make_posterior(adata=adata, indices=indices)

        return -post.reconstruction_error(mode=mode)

    @torch.no_grad()
    def get_latent_representation(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Union[np.ndarray, List[int]]] = None,
        give_mean: bool = True,
        mc_samples: int = 5000,
    ) -> np.ndarray:
        """Return the latent representation for each cell

        Parameters
        ----------
        adata
            AnnData object that has been registered with scvi. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        give_mean
            Give mean of distribution or sample from it
        mc_samples
            For distributions with no closed-form mean (e.g., `logistic normal`), how many Monte Carlo
            samples to take for computing mean.

        Returns
        -------
        latent_representation : np.ndarray
            Low-dimensional representation for each cell

        Examples
        --------

        >>> vae = scvi.model.TOTALVI(adata)
        >>> vae.train(n_epochs=400)
        >>> adata.obsm["X_totalVI"] = vae.get_latent_representation()

        We can also get the latent representation for a subset of cells

        >>> adata_subset = adata[adata.obs.cell_type == "really cool cell type"]
        >>> latent_subset = vae.get_latent_representation(adata_subset)

        """

        if self.is_trained is False:
            raise RuntimeError("Please train the model first.")

        post = self._make_posterior(adata=adata, indices=indices)
        latent = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            y = tensors[_CONSTANTS.PROTEIN_EXP_KEY]
            batch = tensors[_CONSTANTS.BATCH_KEY]
            z = self.model.sample_from_posterior_z(
                x, y, batch, give_mean=give_mean, n_samples=mc_samples
            )
            latent += [z.cpu()]
        return np.array(torch.cat(latent))

    @torch.no_grad()
    def get_latent_library_size(self, adata=None, indices=None, give_mean=True):
        if self.is_trained is False:
            raise RuntimeError("Please train the model first.")

        post = self._make_posterior(adata=adata, indices=indices)
        libraries = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            y = tensors[_CONSTANTS.PROTEIN_EXP_KEY]
            batch = tensors[_CONSTANTS.BATCH_KEY]
            library = self.model.sample_from_posterior_l(
                x, y, batch, give_mean=give_mean
            )
            libraries += [library.cpu()]
        return np.array(torch.cat(libraries))

    @torch.no_grad()
    def get_normalized_expression(
        self,
        adata=None,
        indices=None,
        transform_batch: Optional[int] = None,
        gene_list: Optional[Union[np.ndarray, List[int]]] = None,
        protein_list: Optional[Union[np.ndarray, List[int]]] = None,
        library_size: Optional[Union[float, Literal["latent"]]] = 1,
        n_samples: int = 1,
        sample_protein_mixing: bool = True,
        scale_protein: bool = False,
        include_protein_background: bool = False,
        return_mean: bool = True,
        return_numpy: Optional[bool] = None,
    ) -> Tuple[Union[np.ndarray, pd.DataFrame], Union[np.ndarray, pd.DataFrame]]:
        r"""Returns the normalized gene expression and protein expression

        This is denoted as :math:`\rho_n` in the totalVI paper for genes, and TODO
        for proteins, :math:`(1-\pi_{nt})\alpha_{nt}\beta_{nt}`.

        Parameters
        ----------
        transform_batch
            Batch to condition on.
            If transform_batch is:

            - None, then real observed batch is used
            - int, then batch transform_batch is used
            - List[int], then average over batches in list
        gene_list
            Return frequencies of expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
        protein_list
            Return protein expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
        library_size
            Scale the expression frequencies to a common library size.
            This allows gene expression levels to be interpreted on a common scale of relevant
            magnitude.
        n_samples
            Get sample scale from multiple samples.
        sample_protein_mixing
            Sample mixing bernoulli, setting background to zero
        scale_protein
            Make protein expression sum to 1
        include_protein_background
            Include background component for protein expression
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a `np.ndarray` instead of a `pd.DataFrame`. Includes gene
            names as columns. If either n_samples=1 or return_mean=True, defaults to False.
            Otherwise, it defaults to True.

        Returns
        -------
        - **gene_normalized_expression** - normalized expression for RNA
        - **protein_normalized_expression** - normalized expression for proteins

        If ``n_samples`` > 1 and ``return_mean`` is False, then the shape is ``(samples, cells, genes)``.
        Otherwise, shape is ``(cells, genes)``. Return type is ``pd.DataFrame`` unless ``return_numpy`` is True.

        """

        adata = adata if adata is not None else self.adata
        post = self._make_posterior(adata=adata, indices=indices)

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = adata.var_names
            gene_mask = [True if gene in gene_list else False for gene in all_genes]
        if protein_list is None:
            protein_mask = slice(None)
        else:
            all_proteins = adata.uns["scvi_protein_names"]
            protein_mask = [True if p in protein_list else False for p in all_proteins]

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                logger.warning(
                    "return_numpy must be True if n_samples > 1 and return_mean is False, returning np.ndarray"
                )
            return_numpy = True

        scale_list_gene = []
        scale_list_pro = []
        if not isinstance(transform_batch, Iterable):
            transform_batch = [transform_batch]
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            y = tensors[_CONSTANTS.PROTEIN_EXP_KEY]
            batch_index = tensors[_CONSTANTS.BATCH_KEY]
            label = tensors[_CONSTANTS.LABELS_KEY]
            px_scale = torch.zeros_like(x)
            py_scale = torch.zeros_like(y)
            if n_samples > 1:
                px_scale = torch.stack(n_samples * [px_scale])
                py_scale = torch.stack(n_samples * [py_scale])
            for b in transform_batch:
                outputs = self.model.inference(
                    x,
                    y,
                    batch_index=batch_index,
                    label=label,
                    n_samples=n_samples,
                    transform_batch=b,
                )
                if library_size == "observed":
                    px_scale += outputs["px_"]["rate"]
                else:
                    px_scale += outputs["px_"]["scale"]
                px_scale = px_scale[..., gene_mask]

                py_ = outputs["py_"]
                # probability of background
                protein_mixing = 1 / (1 + torch.exp(-py_["mixing"]))
                if sample_protein_mixing is True:
                    protein_mixing = torch.distributions.Bernoulli(
                        protein_mixing
                    ).sample()
                protein_val = py_["rate_fore"] * (1 - protein_mixing)
                if include_protein_background is True:
                    protein_val += py_["rate_back"] * protein_mixing

                if scale_protein is True:
                    protein_val = torch.nn.functional.normalize(
                        protein_val, p=1, dim=-1
                    )
                protein_val = protein_val[..., protein_mask]
                py_scale += protein_val
            px_scale /= len(transform_batch)
            py_scale /= len(transform_batch)
            scale_list_gene.append(px_scale.cpu())
            scale_list_pro.append(py_scale.cpu())

        if n_samples > 1:
            # concatenate along batch dimension -> result shape = (samples, cells, features)
            scale_list_gene = torch.cat(scale_list_gene, dim=1)
            scale_list_pro = torch.cat(scale_list_pro, dim=1)
            # (cells, features, samples)
            scale_list_gene = scale_list_gene.permute(1, 2, 0)
            scale_list_pro = scale_list_pro.permute(1, 2, 0)
        else:
            scale_list_gene = torch.cat(scale_list_gene, dim=0)
            scale_list_pro = torch.cat(scale_list_pro, dim=0)

        if return_mean is True and n_samples > 1:
            scale_list_gene = torch.mean(scale_list_gene, dim=-1)
            scale_list_pro = torch.mean(scale_list_pro, dim=-1)

        scale_list_gene = scale_list_gene.cpu().numpy()
        scale_list_pro = scale_list_pro.cpu().numpy()

        if return_numpy is None or return_numpy is False:
            gene_df = pd.DataFrame(
                scale_list_gene,
                columns=adata.var_names[gene_mask],
                index=adata.obs_names,
            )
            pro_df = pd.DataFrame(
                scale_list_pro,
                columns=adata.uns["scvi_protein_names"][protein_mask],
                index=adata.obs_names,
            )

            return gene_df, pro_df
        else:
            return scale_list_gene, scale_list_pro

    @torch.no_grad()
    def get_protein_foreground_probability(
        self,
        adata=None,
        indices=None,
        transform_batch: Optional[int] = None,
        protein_list: Optional[Union[np.ndarray, List[int]]] = None,
        n_samples: int = 1,
        return_mean: bool = True,
        return_numpy: Optional[bool] = None,
    ):
        r"""Returns the foreground probability for proteins

        This is denoted as :math:`\pi_{nt}` in the totalVI paper.

        Parameters
        ----------
        transform_batch
            Batch to condition on.
            If transform_batch is:

            - None, then real observed batch is used
            - int, then batch transform_batch is used
            - List[int], then average over batches in list
        protein_list
            Return protein expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
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
        - **foreground_probability** - probability foreground for each protein

        If ``n_samples`` > 1 and ``return_mean`` is False, then the shape is ``(samples, cells, proteins)``.
        Otherwise, shape is ``(cells, proteins)``. Return type is ``pd.DataFrame`` unless ``return_numpy`` is True.

        """

        adata = adata if adata is not None else self.adata
        post = self._make_posterior(adata=adata, indices=indices)

        if protein_list is None:
            protein_mask = slice(None)
        else:
            all_proteins = adata.uns["scvi_protein_names"]
            protein_mask = [True if p in protein_list else False for p in all_proteins]

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                logger.warning(
                    "return_numpy must be True if n_samples > 1 and return_mean is False, returning np.ndarray"
                )
            return_numpy = True

        py_mixings = []
        if (transform_batch is None) or (isinstance(transform_batch, int)):
            transform_batch = [transform_batch]
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            y = tensors[_CONSTANTS.PROTEIN_EXP_KEY]
            batch_index = tensors[_CONSTANTS.BATCH_KEY]
            label = tensors[_CONSTANTS.LABELS_KEY]
            py_mixing = torch.zeros_like(y)
            if n_samples > 1:
                py_mixing = torch.stack(n_samples * [py_mixing])
            for b in transform_batch:
                outputs = self.model.inference(
                    x,
                    y,
                    batch_index=batch_index,
                    label=label,
                    n_samples=n_samples,
                    transform_batch=b,
                )
                py_mixing += torch.sigmoid(outputs["py_"]["mixing"])[..., protein_mask]
            py_mixing /= len(transform_batch)
            py_mixings += [py_mixing.cpu()]
        if n_samples > 1:
            # concatenate along batch dimension -> result shape = (samples, cells, features)
            py_mixings = torch.cat(py_mixings, dim=1)
            # (cells, features, samples)
            py_mixings = py_mixings.permute(1, 2, 0)
        else:
            py_mixings = torch.cat(py_mixings, dim=0)

        if return_mean is True and n_samples > 1:
            py_mixings = torch.mean(py_mixings, dim=-1)

        py_mixings = py_mixings.cpu().numpy()

        if return_numpy is True:
            return 1 - py_mixings
        else:
            pro_names = self.adata.uns["scvi_protein_names"]
            foreground_prob = pd.DataFrame(
                1 - py_mixings, columns=pro_names, index=adata.obs_names
            )
            return foreground_prob

    def differential_expression(
        self,
        groupby,
        adata=None,
        group1=None,
        group2=None,
        mode="change",
        delta=0.5,
        all_stats=True,
        protein_prior_count=0.5,
        scale_protein=False,
        sample_protein_mixing=False,
        include_protein_background=False,
    ):
        if adata is None:
            adata = self.adata
        cell_idx1 = np.asarray(adata.obs[groupby] == group1)
        if group2 is None:
            cell_idx2 = ~cell_idx1
        else:
            cell_idx2 = np.asarray(adata.obs[groupby] == group2)

        def _expression_for_de(
            adata=None,
            indices=None,
            transform_batch: Optional[int] = None,
            scale_protein=False,
            sample_protein_mixing=False,
            include_protein_background=False,
            protein_prior_count=0.5,
        ):
            rna, protein = self.get_normalized_expression(
                adata=adata,
                indices=indices,
                transform_batch=transform_batch,
                return_numpy=True,
                scale_protein=scale_protein,
                sample_protein_mixing=sample_protein_mixing,
                include_protein_background=include_protein_background,
            )
            protein += protein_prior_count

            joint = np.concatenate([rna, protein], axis=1)
            return joint

        model_fn = partial(
            _expression_for_de,
            scale_protein=scale_protein,
            sample_protein_mixing=sample_protein_mixing,
            include_protein_background=include_protein_background,
            protein_prior_count=protein_prior_count,
        )

        dc = DifferentialComputation(model_fn, adata)
        all_info = dc.get_bayes_factors(cell_idx1, cell_idx2, mode=mode, delta=delta)

        col_names = np.concatenate(
            [np.asarray(adata.var_names), adata.uns["scvi_protein_names"]]
        )
        if all_stats is True:
            nan = np.array([np.nan] * len(adata.uns["scvi_protein_names"]))
            (
                mean1,
                mean2,
                nonz1,
                nonz2,
                norm_mean1,
                norm_mean2,
            ) = scrna_raw_counts_properties(adata, cell_idx1, cell_idx2)
            protein_exp = get_from_registry(adata, _CONSTANTS.PROTEIN_EXP_KEY)
            mean2_pro = np.asarray(protein_exp[cell_idx2].mean(0))
            mean1_pro = np.asarray(protein_exp[cell_idx1].mean(0))
            nonz1_pro = np.asarray((protein_exp[cell_idx1] > 0).mean(0))
            nonz2_pro = np.asarray((protein_exp[cell_idx2] > 0).mean(0))
            # TODO implement properties for proteins
            genes_properties_dict = dict(
                raw_mean1=np.concatenate([mean1, mean1_pro]),
                raw_mean2=np.concatenate([mean2, mean2_pro]),
                non_zeros_proportion1=np.concatenate([nonz1, nonz1_pro]),
                non_zeros_proportion2=np.concatenate([nonz2, nonz2_pro]),
                raw_normalized_mean1=np.concatenate([norm_mean1, nan]),
                raw_normalized_mean2=np.concatenate([norm_mean2, nan]),
            )
            all_info = {**all_info, **genes_properties_dict}

        res = pd.DataFrame(all_info, index=col_names)
        sort_key = "proba_de" if mode == "change" else "bayes_factor"
        res = res.sort_values(by=sort_key, ascending=False)
        return res

    @torch.no_grad()
    def posterior_predictive_sample(
        self,
        adata=None,
        indices=None,
        n_samples: int = 1,
        gene_list: Union[list, np.ndarray] = None,
        protein_list: Union[list, np.ndarray] = None,
    ) -> np.ndarray:
        r"""Generate observation samples from the posterior predictive distribution

        The posterior predictive distribution is written as :math:`p(\hat{x}, \hat{y} \mid x, y)`.

        Parameters
        ----------
        n_samples
            Number of required samples for each cell
        gene_list
            Indices or names of genes of interest

        Returns
        -------
        x_new : :py:class:`np.ndarray`
            tensor with shape (n_cells, n_genes, n_samples)
        """
        assert self.model.reconstruction_loss_gene in ["nb"]

        adata = adata if adata is not None else self.adata
        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = adata.var_names
            gene_mask = [True if gene in gene_list else False for gene in all_genes]
        if protein_list is None:
            protein_mask = slice(None)
        else:
            all_proteins = adata.uns["scvi_protein_names"]
            protein_mask = [True if p in protein_list else False for p in all_proteins]

        post = self._make_posterior(adata=adata, indices=indices)

        posterior_list = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            batch_idx = tensors[_CONSTANTS.BATCH_KEY]
            labels = tensors[_CONSTANTS.LABELS_KEY]
            y = tensors[_CONSTANTS.PROTEIN_EXP_KEY]
            with torch.no_grad():
                outputs = self.model.inference(
                    x, y, batch_index=batch_idx, label=labels, n_samples=n_samples
                )
            px_ = outputs["px_"]
            py_ = outputs["py_"]

            pi = 1 / (1 + torch.exp(-py_["mixing"]))
            mixing_sample = torch.distributions.Bernoulli(pi).sample()
            protein_rate = (
                py_["rate_fore"] * (1 - mixing_sample)
                + py_["rate_back"] * mixing_sample
            )
            rate = torch.cat(
                (px_["rate"][..., gene_mask], protein_rate[..., protein_mask]), dim=-1
            )
            if len(px_["r"].size()) == 2:
                px_dispersion = px_["r"]
            else:
                px_dispersion = torch.ones_like(x) * px_["r"]
            if len(py_["r"].size()) == 2:
                py_dispersion = py_["r"]
            else:
                py_dispersion = torch.ones_like(y) * py_["r"]

            dispersion = torch.cat(
                (px_dispersion[..., gene_mask], py_dispersion[..., protein_mask]),
                dim=-1,
            )

            # This gamma is really l*w using scVI manuscript notation
            p = rate / (rate + dispersion)
            r = dispersion
            l_train = torch.distributions.Gamma(r, (1 - p) / p).sample()
            data = torch.distributions.Poisson(l_train).sample().cpu().numpy()
            # """
            # In numpy (shape, scale) => (concentration, rate), with scale = p /(1 - p)
            # rate = (1 - p) / p  # = 1/scale # used in pytorch
            # """
            posterior_list += [data]

            posterior_list[-1] = np.transpose(posterior_list[-1], (1, 2, 0))

        posterior_list = np.concatenate(posterior_list, axis=0)

        return posterior_list

    @torch.no_grad()
    def _get_denoised_samples(
        self,
        adata=None,
        indices=None,
        n_samples: int = 25,
        batch_size: int = 64,
        rna_size_factor: int = 1000,
        transform_batch: Optional[int] = None,
    ) -> np.ndarray:
        """Return samples from an adjusted posterior predictive.

        Parameters
        ----------
        n_samples
            How may samples per cell
        batch_size
            Mini-batch size for sampling. Lower means less GPU memory footprint
        rna_size_factor
            size factor for RNA prior to sampling gamma distribution
        transform_batch
            int of which batch to condition on for all cells

        Returns
        -------

        """

        post = self._make_posterior(adata=adata, indices=indices)

        posterior_list = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            batch_idx = tensors[_CONSTANTS.BATCH_KEY]
            label = tensors[_CONSTANTS.LABELS_KEY]
            y = tensors[_CONSTANTS.PROTEIN_EXP_KEY]
            with torch.no_grad():
                outputs = self.model.inference(
                    x,
                    y,
                    batch_index=batch_idx,
                    label=label,
                    n_samples=n_samples,
                    transform_batch=transform_batch,
                )
            px_ = outputs["px_"]
            py_ = outputs["py_"]

            pi = 1 / (1 + torch.exp(-py_["mixing"]))
            mixing_sample = torch.distributions.Bernoulli(pi).sample()
            protein_rate = py_["rate_fore"]
            rate = torch.cat((rna_size_factor * px_["scale"], protein_rate), dim=-1)
            if len(px_["r"].size()) == 2:
                px_dispersion = px_["r"]
            else:
                px_dispersion = torch.ones_like(x) * px_["r"]
            if len(py_["r"].size()) == 2:
                py_dispersion = py_["r"]
            else:
                py_dispersion = torch.ones_like(y) * py_["r"]

            dispersion = torch.cat((px_dispersion, py_dispersion), dim=-1)

            # This gamma is really l*w using scVI manuscript notation
            p = rate / (rate + dispersion)
            r = dispersion
            l_train = torch.distributions.Gamma(r, (1 - p) / p).sample()
            data = l_train.cpu().numpy()
            # make background 0
            data[:, :, self.adata.shape[1] :] = (
                data[:, :, self.adata.shape[1] :] * (1 - mixing_sample).cpu().numpy()
            )
            posterior_list += [data]

            posterior_list[-1] = np.transpose(posterior_list[-1], (1, 2, 0))

        return np.concatenate(posterior_list, axis=0)

    @torch.no_grad()
    def get_feature_correlation_matrix(
        self,
        adata=None,
        indices=None,
        n_samples: int = 10,
        batch_size: int = 64,
        rna_size_factor: int = 1000,
        transform_batch: Optional[Union[int, List[int]]] = None,
        correlation_type: Literal["spearman", "pearson"] = "spearman",
        log_transform: bool = False,
    ) -> pd.DataFrame:
        """Generate gene-gene correlation matrix using scvi uncertainty and expression

        Parameters
        ----------
        n_samples
            How may samples per cell
        batch_size
            Mini-batch size for sampling. Lower means less GPU memory footprint
        rna_size_factor
            size factor for RNA prior to sampling gamma distribution
        transform_batch
            Batches to condition on.
            If transform_batch is:

            - None, then real observed batch is used
            - int, then batch transform_batch is used
            - list of int, then values are averaged over provided batches.
        correlation_type
            One of "pearson", "spearman"
        log_transform
            Whether to log transform denoised values prior to correlation calculation

        Returns
        -------
        Gene-protein-gene-protein correlation matrix
        """

        from scipy.stats import spearmanr

        if (transform_batch is None) or (isinstance(transform_batch, int)):
            transform_batch = [transform_batch]
        corr_mats = []
        for b in transform_batch:
            denoised_data = self.generate_denoised_samples(
                n_samples=n_samples,
                batch_size=batch_size,
                rna_size_factor=rna_size_factor,
                transform_batch=b,
            )
            flattened = np.zeros(
                (denoised_data.shape[0] * n_samples, denoised_data.shape[1])
            )
            for i in range(n_samples):
                flattened[
                    denoised_data.shape[0] * (i) : denoised_data.shape[0] * (i + 1)
                ] = denoised_data[:, :, i]
            if log_transform is True:
                flattened[:, : self.gene_dataset.n_genes] = np.log(
                    flattened[:, : self.gene_dataset.n_genes] + 1e-8
                )
                flattened[:, self.gene_dataset.n_genes :] = np.log1p(
                    flattened[:, self.gene_dataset.n_genes :]
                )
            if correlation_type == "pearson":
                corr_matrix = np.corrcoef(flattened, rowvar=False)
            else:
                corr_matrix = spearmanr(flattened, axis=0)
            corr_mats.append(corr_matrix)
        corr_matrix = np.mean(np.stack(corr_mats), axis=0)

        names = np.concatenate(
            [np.asarray(self.adata.var_names), self.adata.uns["scvi_protein_names"]]
        )
        return pd.DataFrame(corr_matrix, index=names, columns=names)

    # TODO
    @torch.no_grad()
    def get_likelihood_parameters(
        self,
        adata=None,
        indices=None,
        n_samples: Optional[int] = 1,
        give_mean: Optional[bool] = False,
    ) -> Dict[str, np.ndarray]:

        r"""Estimates for the parameters of the likelihood :math:`p(x \mid z)`.


        """

        post = self._make_posterior(adata=adata, indices=indices)

        rna_dropout_list = []
        rna_mean_list = []
        rna_dispersion_list = []
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

            n_batch = px_rate.size(0) if n_samples == 1 else px_rate.size(1)
            rna_dispersion_list += [
                np.repeat(np.array(px_r.cpu())[np.newaxis, :], n_batch, axis=0)
            ]
            rna_mean_list += [np.array(px_rate.cpu())]
            rna_dropout_list += [np.array(px_dropout.cpu())]

        dropout = np.concatenate(rna_dropout_list)
        means = np.concatenate(rna_mean_list)
        dispersions = np.concatenate(rna_dispersion_list)
        if give_mean and n_samples > 1:
            dropout = dropout.mean(0)
            means = means.mean(0)

        return_dict = {}
        return_dict["RNA_mean"] = means

        if self.model.reconstruction_loss_gene == "zinb":
            return_dict["RNA_dropout"] = dropout
            return_dict["RNA_dispersions"] = dispersions
        if self.model.reconstruction_loss_gene == "nb":
            return_dict["RNA_dispersions"] = dispersions

        return return_dict


class BaseModelClass(ABC):
    def __init__(self, adata, use_cuda):
        assert (
            "scvi_data_registry" in adata.uns.keys()
        ), "Please setup your AnnData with scvi.dataset.setup_anndata(adata) first"
        self.adata = adata
        self._check_anndata(adata)

        # TODO make abstract properties
        self.summary_stats = adata.uns["scvi_summary_stats"]
        self.is_trained = False
        # just a regular property
        self.use_cuda = use_cuda and torch.cuda.is_available()
        self._model_summary_string = ""

        self._posterior_class = None
        self._trainer_class = None

        # all methods need a batch_size and it needs to be passed to make posterior

    def _make_posterior(self, adata=None, indices=None, batch_size=128):
        if adata is None:
            adata = self.adata
        if indices is None:
            indices = np.arange(adata.n_obs)
        post = self._posterior_class(
            self.model,
            adata,
            shuffle=False,
            indices=indices,
            use_cuda=self.use_cuda,
            batch_size=batch_size,
        ).sequential()
        return post

    # @abstractmethod
    def _check_anndata(self, adata):
        pass

    @abstractmethod
    def train(self):
        pass

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

    def __repr__(self,):
        summary_string = self._model_summary_string + "\nTraining status: {}".format(
            "Trained" if self.is_trained else "Not Trained"
        )
        return summary_string
