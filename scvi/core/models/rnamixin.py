import logging
from functools import partial
from typing import Dict, Iterable, Optional, Sequence, Union

import numpy as np
import pandas as pd
import torch
from anndata import AnnData

from scvi import _CONSTANTS
from scvi._compat import Literal
from scvi._docs import doc_differential_expression
from scvi._utils import _doc_params
from scvi.core.distributions import NegativeBinomial, ZeroInflatedNegativeBinomial
from scvi.core.models._utils import _de_core
from scvi.model._utils import (
    _get_batch_code_from_category,
    _get_var_names_from_setup_anndata,
    scrna_raw_counts_properties,
)

logger = logging.getLogger(__name__)

Number = Union[int, float]


class RNASeqMixin:
    @torch.no_grad()
    def get_normalized_expression(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        transform_batch: Optional[Sequence[Union[Number, str]]] = None,
        gene_list: Optional[Sequence[str]] = None,
        library_size: Union[float, Literal["latent"]] = 1,
        n_samples: int = 1,
        batch_size: Optional[int] = None,
        return_mean: bool = True,
        return_numpy: Optional[bool] = None,
    ) -> Union[np.ndarray, pd.DataFrame]:
        r"""
        Returns the normalized (decoded) gene expression.

        This is denoted as :math:`\rho_n` in the scVI paper.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        transform_batch
            Batch to condition on.
            If transform_batch is:

            - None, then real observed batch is used.
            - int, then batch transform_batch is used.
        gene_list
            Return frequencies of expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
        library_size
            Scale the expression frequencies to a common library size.
            This allows gene expression levels to be interpreted on a common scale of relevant
            magnitude. If set to `"latent"`, use the latent libary size.
        n_samples
            Number of posterior samples to use for estimation.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`. DataFrame includes
            gene names as columns. If either `n_samples=1` or `return_mean=True`, defaults to `False`.
            Otherwise, it defaults to `True`.

        Returns
        -------
        If `n_samples` > 1 and `return_mean` is False, then the shape is `(samples, cells, genes)`.
        Otherwise, shape is `(cells, genes)`. In this case, return type is :class:`~pandas.DataFrame` unless `return_numpy` is True.
        """
        adata = self._validate_anndata(adata)
        scdl = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)

        transform_batch = _get_batch_code_from_category(adata, transform_batch)

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = _get_var_names_from_setup_anndata(adata)
            gene_mask = [True if gene in gene_list else False for gene in all_genes]

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                logger.warning(
                    "return_numpy must be True if n_samples > 1 and return_mean is False, returning np.ndarray"
                )
            return_numpy = True
        if indices is None:
            indices = np.arange(adata.n_obs)
        if library_size == "latent":
            model_fn = self.model.get_sample_rate
            scaling = 1
        else:
            model_fn = self.model.get_sample_scale
            scaling = library_size

        exprs = []
        for tensors in scdl:
            x = tensors[_CONSTANTS.X_KEY]
            batch_idx = tensors[_CONSTANTS.BATCH_KEY]
            labels = tensors[_CONSTANTS.LABELS_KEY]
            per_batch_exprs = []
            for batch in transform_batch:
                output = model_fn(
                    x,
                    batch_index=batch_idx,
                    y=labels,
                    n_samples=n_samples,
                    transform_batch=batch,
                )[..., gene_mask]
                output *= scaling
                output = output.cpu().numpy()
                per_batch_exprs.append(output)
            per_batch_exprs = np.stack(
                per_batch_exprs
            )  # shape is (len(transform_batch) x batch_size x n_var)
            exprs += [per_batch_exprs.mean(0)]

        if n_samples > 1:
            # The -2 axis correspond to cells.
            exprs = np.concatenate(exprs, axis=-2)
        else:
            exprs = np.concatenate(exprs, axis=0)
        if n_samples > 1 and return_mean:
            exprs = exprs.mean(0)

        if return_numpy is None or return_numpy is False:
            return pd.DataFrame(
                exprs,
                columns=adata.var_names[gene_mask],
                index=adata.obs_names[indices],
            )
        else:
            return exprs

    @_doc_params(
        doc_differential_expression=doc_differential_expression,
    )
    def differential_expression(
        self,
        adata: Optional[AnnData] = None,
        groupby: Optional[str] = None,
        group1: Optional[Iterable[str]] = None,
        group2: Optional[str] = None,
        idx1: Optional[Union[Sequence[int], Sequence[bool]]] = None,
        idx2: Optional[Union[Sequence[int], Sequence[bool]]] = None,
        mode: Literal["vanilla", "change"] = "change",
        delta: float = 0.25,
        batch_size: Optional[int] = None,
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: Optional[Iterable[str]] = None,
        batchid2: Optional[Iterable[str]] = None,
        fdr_target: float = 0.05,
        **kwargs,
    ) -> pd.DataFrame:
        r"""
        A unified method for differential expression analysis.

        Implements `"vanilla"` DE [Lopez18]_ and `"change"` mode DE [Boyeau19]_.

        Parameters
        ----------
        {doc_differential_expression}
        **kwargs
            Keyword args for :func:`scvi.core.utils.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential expression DataFrame.
        """
        adata = self._validate_anndata(adata)

        col_names = _get_var_names_from_setup_anndata(adata)
        model_fn = partial(
            self.get_normalized_expression,
            return_numpy=True,
            n_samples=1,
            batch_size=batch_size,
        )
        result = _de_core(
            adata,
            model_fn,
            groupby,
            group1,
            group2,
            idx1,
            idx2,
            all_stats,
            scrna_raw_counts_properties,
            col_names,
            mode,
            batchid1,
            batchid2,
            delta,
            batch_correction,
            fdr_target,
            **kwargs,
        )

        return result

    @torch.no_grad()
    def posterior_predictive_sample(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        n_samples: int = 1,
        gene_list: Optional[Sequence[str]] = None,
        batch_size: Optional[int] = None,
    ) -> np.ndarray:
        r"""
        Generate observation samples from the posterior predictive distribution.

        The posterior predictive distribution is written as :math:`p(\hat{x} \mid x)`.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        n_samples
            Number of samples for each cell.
        gene_list
            Names of genes of interest.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        x_new : :py:class:`torch.Tensor`
            tensor with shape (n_cells, n_genes, n_samples)
        """
        if self.model.gene_likelihood not in ["zinb", "nb", "poisson"]:
            raise ValueError("Invalid gene_likelihood.")

        adata = self._validate_anndata(adata)
        scdl = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)

        if indices is None:
            indices = np.arange(adata.n_obs)

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = _get_var_names_from_setup_anndata(adata)
            gene_mask = [True if gene in gene_list else False for gene in all_genes]

        x_new = []
        for tensors in scdl:
            x = tensors[_CONSTANTS.X_KEY]
            batch_idx = tensors[_CONSTANTS.BATCH_KEY]
            labels = tensors[_CONSTANTS.LABELS_KEY]
            outputs = self.model.inference(
                x, batch_index=batch_idx, y=labels, n_samples=n_samples
            )
            px_r = outputs["px_r"]
            px_rate = outputs["px_rate"]
            px_dropout = outputs["px_dropout"]

            if self.model.gene_likelihood == "poisson":
                l_train = px_rate
                l_train = torch.clamp(l_train, max=1e8)
                dist = torch.distributions.Poisson(
                    l_train
                )  # Shape : (n_samples, n_cells_batch, n_genes)
            elif self.model.gene_likelihood == "nb":
                dist = NegativeBinomial(mu=px_rate, theta=px_r)
            elif self.model.gene_likelihood == "zinb":
                dist = ZeroInflatedNegativeBinomial(
                    mu=px_rate, theta=px_r, zi_logits=px_dropout
                )
            else:
                raise ValueError(
                    "{} reconstruction error not handled right now".format(
                        self.model.gene_likelihood
                    )
                )
            if n_samples > 1:
                exprs = dist.sample().permute(
                    [1, 2, 0]
                )  # Shape : (n_cells_batch, n_genes, n_samples)
            else:
                exprs = dist.sample()

            if gene_list is not None:
                exprs = exprs[:, gene_mask, ...]

            x_new.append(exprs.cpu())
        x_new = torch.cat(x_new)  # Shape (n_cells, n_genes, n_samples)

        return x_new.numpy()

    @torch.no_grad()
    def _get_denoised_samples(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        n_samples: int = 25,
        batch_size: int = 64,
        rna_size_factor: int = 1000,
        transform_batch: Optional[Sequence[int]] = None,
    ) -> np.ndarray:
        """
        Return samples from an adjusted posterior predictive.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        n_samples
            Number of posterior samples to use for estimation.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        rna_size_factor
            size factor for RNA prior to sampling gamma distribution.
        transform_batch
            int of which batch to condition on for all cells.

        Returns
        -------
        denoised_samples
        """
        adata = self._validate_anndata(adata)
        scdl = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)

        data_loader_list = []
        for tensors in scdl:
            x = tensors[_CONSTANTS.X_KEY]
            batch_idx = tensors[_CONSTANTS.BATCH_KEY]
            labels = tensors[_CONSTANTS.LABELS_KEY]

            outputs = self.model.inference(
                x,
                batch_index=batch_idx,
                y=labels,
                n_samples=n_samples,
                transform_batch=transform_batch,
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
            data_loader_list += [data]

            data_loader_list[-1] = np.transpose(data_loader_list[-1], (1, 2, 0))

        return np.concatenate(data_loader_list, axis=0)

    @torch.no_grad()
    def get_feature_correlation_matrix(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        n_samples: int = 10,
        batch_size: int = 64,
        rna_size_factor: int = 1000,
        transform_batch: Optional[Sequence[Union[Number, str]]] = None,
        correlation_type: Literal["spearman", "pearson"] = "spearman",
    ) -> pd.DataFrame:
        """
        Generate gene-gene correlation matrix using scvi uncertainty and expression.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        n_samples
            Number of posterior samples to use for estimation.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        rna_size_factor
            size factor for RNA prior to sampling gamma distribution.
        transform_batch
            Batches to condition on.
            If transform_batch is:

            - None, then real observed batch is used.
            - int, then batch transform_batch is used.
            - list of int, then values are averaged over provided batches.
        correlation_type
            One of "pearson", "spearman".

        Returns
        -------
        Gene-gene correlation matrix
        """
        from scipy.stats import spearmanr

        adata = self._validate_anndata(adata)

        transform_batch = _get_batch_code_from_category(adata, transform_batch)

        corr_mats = []
        for b in transform_batch:
            denoised_data = self._get_denoised_samples(
                adata=adata,
                indices=indices,
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
        var_names = _get_var_names_from_setup_anndata(adata)
        return pd.DataFrame(corr_matrix, index=var_names, columns=var_names)

    @torch.no_grad()
    def get_likelihood_parameters(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        n_samples: Optional[int] = 1,
        give_mean: Optional[bool] = False,
        batch_size: Optional[int] = None,
    ) -> Dict[str, np.ndarray]:
        r"""
        Estimates for the parameters of the likelihood :math:`p(x \mid z)`

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        n_samples
            Number of posterior samples to use for estimation.
        give_mean
            Return expected value of parameters or a samples
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        """
        adata = self._validate_anndata(adata)
        scdl = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)

        dropout_list = []
        mean_list = []
        dispersion_list = []
        for tensors in scdl:
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

        if self.model.gene_likelihood == "zinb":
            return_dict["dropout"] = dropout
            return_dict["dispersions"] = dispersions
        if self.model.gene_likelihood == "nb":
            return_dict["dispersions"] = dispersions

        return return_dict

    @torch.no_grad()
    def get_latent_library_size(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        give_mean: bool = True,
        batch_size: Optional[int] = None,
    ) -> np.ndarray:
        r"""
        Returns the latent library size for each cell.

        This is denoted as :math:`\ell_n` in the scVI paper.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        give_mean
            Return the mean or a sample from the posterior distribution.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        """
        if self.is_trained_ is False:
            raise RuntimeError("Please train the model first.")
        adata = self._validate_anndata(adata)
        scdl = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)
        libraries = []
        for tensors in scdl:
            x = tensors[_CONSTANTS.X_KEY]
            b = tensors[_CONSTANTS.BATCH_KEY]
            library = self.model.sample_from_posterior_l(x, b, give_mean=give_mean)
            libraries += [library.cpu()]
        return np.array(torch.cat(libraries))
