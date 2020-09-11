import torch
import pickle
import os
import logging
import pandas as pd
import numpy as np
import inspect

from anndata import AnnData
from functools import partial
from scvi.core._distributions import NegativeBinomial, ZeroInflatedNegativeBinomial
from scvi.core.utils import DifferentialComputation
from scvi.models._utils import (
    scrna_raw_counts_properties,
    _get_var_names_from_setup_anndata,
    _get_batch_code_from_category,
)
from scvi.dataset._utils import (
    _check_anndata_setup_equivalence,
    _check_nonnegative_integers,
)
from scvi import _CONSTANTS
from typing import Optional, Union, List, Dict, Sequence, Iterable
from scvi._compat import Literal
from scvi.core.trainers import UnsupervisedTrainer
from abc import ABC, abstractmethod
from scvi.dataset import get_from_registry, transfer_anndata_setup

logger = logging.getLogger(__name__)


class VAEMixin:
    def train(
        self,
        n_epochs=400,
        train_size=0.9,
        test_size=None,
        lr=1e-3,
        n_iter_kl_warmup=None,
        n_epochs_kl_warmup=400,
        frequency=None,
        train_fun_kwargs={},
        **kwargs,
    ):
        train_fun_kwargs = dict(train_fun_kwargs)
        if self.is_trained_ is False:
            self.trainer = UnsupervisedTrainer(
                self.model,
                self.adata,
                train_size=train_size,
                test_size=test_size,
                n_iter_kl_warmup=n_iter_kl_warmup,
                n_epochs_kl_warmup=n_epochs_kl_warmup,
                frequency=frequency,
                use_cuda=self.use_cuda,
                **kwargs,
            )
            self.train_indices_ = self.trainer.train_set.indices
            self.test_indices_ = self.trainer.test_set.indices
            self.validation_indices_ = self.trainer.validation_set.indices
        # for autotune
        if "n_epochs" not in train_fun_kwargs:
            train_fun_kwargs["n_epochs"] = n_epochs
        if "lr" not in train_fun_kwargs:
            train_fun_kwargs["lr"] = lr
        self.trainer.train(**train_fun_kwargs)
        self.is_trained_ = True

    @torch.no_grad()
    def get_elbo(self, adata=None, indices=None, batch_size=128):

        adata = self._validate_anndata(adata)
        post = self._make_posterior(adata=adata, indices=indices, batch_size=batch_size)

        return -post.elbo()

    @torch.no_grad()
    def get_marginal_ll(
        self, adata=None, indices=None, n_mc_samples=1000, batch_size=128
    ):

        adata = self._validate_anndata(adata)
        post = self._make_posterior(adata=adata, indices=indices, batch_size=batch_size)

        return -post.marginal_ll(n_mc_samples=n_mc_samples)

    @torch.no_grad()
    def get_reconstruction_error(self, adata=None, indices=None, batch_size=128):

        adata = self._validate_anndata(adata)
        post = self._make_posterior(adata=adata, indices=indices, batch_size=batch_size)

        return -post.reconstruction_error()

    @torch.no_grad()
    def get_latent_representation(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Union[np.ndarray, List[int]]] = None,
        give_mean: bool = True,
        mc_samples: int = 5000,
        batch_size=128,
    ) -> np.ndarray:
        """
        Return the latent representation for each cell.

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
        batch_size
            Minibatch size for data loading into model

        Returns
        -------
        latent_representation : np.ndarray
            Low-dimensional representation for each cell
        """
        if self.is_trained_ is False:
            raise RuntimeError("Please train the model first.")

        adata = self._validate_anndata(adata)
        post = self._make_posterior(adata=adata, indices=indices, batch_size=batch_size)
        latent = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            z = self.model.sample_from_posterior_z(
                x, give_mean=give_mean, n_samples=mc_samples
            )
            latent += [z.cpu()]
        return np.array(torch.cat(latent))


class RNASeqMixin:
    @torch.no_grad()
    def get_normalized_expression(
        self,
        adata=None,
        indices=None,
        transform_batch: Optional[Union[str, int]] = None,
        gene_list: Optional[Union[np.ndarray, List[int]]] = None,
        library_size: Optional[Union[float, Literal["latent"]]] = 1,
        n_samples: int = 1,
        batch_size=128,
        return_mean: bool = True,
        return_numpy: Optional[bool] = None,
    ) -> Union[np.ndarray, pd.DataFrame]:
        # whatever is here will pass it into cat and batch
        r"""
        Returns the normalized (decoded) gene expression.

        This is denoted as :math:`\rho_n` in the scVI paper.

        Parameters
        ----------
        adata
            AnnData object that has been registered with scvi. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
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
        batch_size
            Minibatch size for data loading into model
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
        adata = self._validate_anndata(adata)
        post = self._make_posterior(adata=adata, indices=indices, batch_size=batch_size)
        if transform_batch is not None:
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
                exprs,
                columns=adata.var_names[gene_mask],
                index=adata.obs_names[indices],
            )
        else:
            return exprs

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
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: Optional[Iterable[str]] = None,
        batchid2: Optional[Iterable[str]] = None,
        **kwargs,
    ) -> pd.DataFrame:
        r"""
        A unified method for differential expression inference.

        Implements `"vanilla"` DE [Lopez18]_ and `"change"` mode DE [Boyeau19]_.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData.
            If None, defaults to the AnnData object used to initialize the model.
        groupby
            The key of the observations grouping to consider.
        group1
            Subset of groups, e.g. [`'g1'`, `'g2'`, `'g3'`], to which comparison
            shall be restricted, or all groups in `groupby` (default).
        group2
            If `None`, compare each group in `group1` to the union of the rest of the groups
            in `groupby`. If a group identifier, compare with respect to this group.
        idx1
            Boolean mask or indices for `group1`. `idx1` and `idx2` can be used as an alternative
            to the AnnData keys. If `idx1` is not `None`, this option overrides `group1`
            and `group2`.
        idx2
            Boolean mask or indices for `group2`. By default, includes all cells not specified in
            `idx1`.
        mode
            Method for differential expression. See user guide for full explanation.
        delta
            specific case of region inducing differential expression.
            In this case, we suppose that :math:`R \setminus [-\delta, \delta]` does not induce differential expression
            (change model default case).
        all_stats
            Concatenate count statistics (e.g., mean expression group 1) to DE results.
        batch_correction
            Whether to correct for batch effects in DE inference.
        batchid1
            Subset of categories from `batch_key` registered in :func:`~scvi.dataset.setup_anndata`,
            e.g. [`'batch1'`, `'batch2'`, `'batch3'`], for group1. By default all categories are used.
            Only used if `batch_correction` is `True`.
        batchid2
            Same as `batchid1` for group2. `batchid2` must either have null intersection with `batchid1`,
            or be exactly equal to `batchid1`.
        **kwargs
            Keyword args for :func:`scvi.core.utils.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential expression DataFrame.
        """
        adata = self._validate_anndata(adata)

        if group1 is None and idx1 is None:
            group1 = adata.obs[groupby].cat.categories.tolist()

        if isinstance(group1, str):
            group1 = [group1]

        # make a temp obs key using indices
        temp_key = None
        if idx1 is not None:
            g1_key = "one"
            obs_col = np.array(["None"] * adata.shape[0], dtype=str)
            obs_col[idx1] = g1_key
            group2 = None if idx2 is None else "two"
            if idx2 is not None:
                obs_col[idx2] = group2
            temp_key = "_scvi_temp_de"
            adata.obs[temp_key] = obs_col
            groupby = temp_key
            group1 = [g1_key]

        df_results = []
        gene_names = _get_var_names_from_setup_anndata(adata)
        model_fn = partial(
            self.get_normalized_expression, return_numpy=True, n_samples=1
        )
        dc = DifferentialComputation(model_fn, adata)
        for g1 in group1:
            cell_idx1 = adata.obs[groupby] == g1
            if group2 is None:
                cell_idx2 = ~cell_idx1
            else:
                cell_idx2 = adata.obs[groupby] == group2

            all_info = dc.get_bayes_factors(
                cell_idx1,
                cell_idx2,
                mode=mode,
                delta=delta,
                batchid1=batchid1,
                batchid2=batchid2,
                use_observed_batches=not batch_correction,
                **kwargs,
            )

            if all_stats is True:
                genes_properties_dict = scrna_raw_counts_properties(
                    adata, cell_idx1, cell_idx2
                )
                all_info = {**all_info, **genes_properties_dict}

            res = pd.DataFrame(all_info, index=gene_names)
            sort_key = "proba_de" if mode == "change" else "bayes_factor"
            res = res.sort_values(by=sort_key, ascending=False)
            if idx1 is not None:
                res["comparison"] = "{} vs {}".format(g1, group2)
            df_results.append(res)

        if temp_key is not None:
            del adata.obs[temp_key]

        result = pd.concat(df_results, axis=0)

        return result

    @torch.no_grad()
    def posterior_predictive_sample(
        self,
        adata=None,
        indices=None,
        n_samples: int = 1,
        gene_list: Union[list, np.ndarray] = None,
        batch_size=128,
    ) -> np.ndarray:
        r"""
        Generate observation samples from the posterior predictive distribution.

        The posterior predictive distribution is written as :math:`p(\hat{x} \mid x)`.

        Parameters
        ----------
        adata
            AnnData object that has been registered with scvi. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        n_samples
            Number of required samples for each cell
        gene_list
            Names of genes of interest
        batch_size
            Minibatch size for data loading into model

        Returns
        -------
        x_new : :py:class:`torch.Tensor`
            tensor with shape (n_cells, n_genes, n_samples)
        """
        if self.model.gene_likelihood not in ["zinb", "nb", "poisson"]:
            raise ValueError("Invalid gene_likelihood.")

        adata = self._validate_anndata(adata)
        post = self._make_posterior(adata=adata, indices=indices, batch_size=batch_size)

        if indices is None:
            indices = np.arange(adata.n_obs)

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = _get_var_names_from_setup_anndata(adata)
            gene_mask = [True if gene in gene_list else False for gene in all_genes]

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
        adata=None,
        indices=None,
        n_samples: int = 25,
        batch_size: int = 64,
        rna_size_factor: int = 1000,
        transform_batch: Optional[int] = None,
    ) -> np.ndarray:
        """
        Return samples from an adjusted posterior predictive.

        Parameters
        ----------
        adata
            AnnData object that has been registered with scvi. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        n_samples
            How may samples per cell
        batch_size
            Minibatch size for data loading into model
        rna_size_factor
            size factor for RNA prior to sampling gamma distribution
        transform_batch
            int of which batch to condition on for all cells

        Returns
        -------
        denoised_samples
        """
        adata = self._validate_anndata(adata)
        post = self._make_posterior(adata=adata, indices=indices, batch_size=batch_size)

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
        """
        Generate gene-gene correlation matrix using scvi uncertainty and expression.

        Parameters
        ----------
        adata
            AnnData object that has been registered with scvi. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        n_samples
            How may samples per cell
        batch_size
            Minibatch size for data loading into model
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

        adata = self._validate_anndata(adata)

        if (transform_batch is None) or (isinstance(transform_batch, int)):
            transform_batch = [transform_batch]
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
        adata=None,
        indices=None,
        n_samples: Optional[int] = 1,
        give_mean: Optional[bool] = False,
        batch_size=128,
    ) -> Dict[str, np.ndarray]:
        r"""Estimates for the parameters of the likelihood :math:`p(x \mid z)`."""
        adata = self._validate_anndata(adata)
        post = self._make_posterior(adata=adata, indices=indices, batch_size=batch_size)

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

        if self.model.gene_likelihood == "zinb":
            return_dict["dropout"] = dropout
            return_dict["dispersions"] = dispersions
        if self.model.gene_likelihood == "nb":
            return_dict["dispersions"] = dispersions

        return return_dict

    @torch.no_grad()
    def get_latent_library_size(
        self, adata=None, indices=None, give_mean=True, batch_size=128
    ):
        if self.is_trained_ is False:
            raise RuntimeError("Please train the model first.")
        adata = self._validate_anndata(adata)
        post = self._make_posterior(adata=adata, indices=indices, batch_size=batch_size)
        libraries = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            library = self.model.sample_from_posterior_l(x, give_mean=give_mean)
            libraries += [library.cpu()]
        return np.array(torch.cat(libraries))


class BaseModelClass(ABC):
    def __init__(self, adata=None, use_cuda=False):
        if adata is not None:
            if "_scvi" not in adata.uns.keys():
                raise ValueError(
                    "Please setup your AnnData with scvi.dataset.setup_anndata(adata) first"
                )
            self.adata = adata
            self.scvi_setup_dict_ = adata.uns["_scvi"]
            self.summary_stats = self.scvi_setup_dict_["summary_stats"]
            self._validate_anndata(adata, copy_if_view=False)

        self.is_trained_ = False
        self.use_cuda = use_cuda and torch.cuda.is_available()
        self._model_summary_string = ""
        self.train_indices_ = None
        self.test_indices_ = None
        self.validation_indices_ = None

    def _make_posterior(
        self, adata: AnnData, indices=None, batch_size=128, **posterior_kwargs
    ):
        if indices is None:
            indices = np.arange(adata.n_obs)
        post = self._posterior_class(
            self.model,
            adata,
            shuffle=False,
            indices=indices,
            use_cuda=self.use_cuda,
            batch_size=batch_size,
            **posterior_kwargs,
        ).sequential()
        return post

    def _validate_anndata(
        self, adata: Optional[AnnData] = None, copy_if_view: bool = True
    ):
        if adata is None:
            adata = self.adata
        if adata.is_view:
            logger.warning("Input anndata is a view.")
            if copy_if_view:
                logger.info("Making copy of anndata.")
                adata = adata.copy()
        if "_scvi" not in adata.uns_keys():
            logger.info(
                "Input adata not setup with scvi. "
                + "attempting to transfer anndata setup"
            )
            transfer_anndata_setup(self.scvi_setup_dict_, adata)
        is_nonneg_int = _check_nonnegative_integers(
            get_from_registry(adata, _CONSTANTS.X_KEY)
        )
        if not is_nonneg_int:
            logger.warning(
                "Make sure the registered X field in anndata contains unnormalized count data."
            )

        _check_anndata_setup_equivalence(self.scvi_setup_dict_, adata)

        return adata

    @property
    @abstractmethod
    def _posterior_class(self):
        pass

    @property
    @abstractmethod
    def _trainer_class(self):
        pass

    @abstractmethod
    def train(self):
        pass

    @property
    def is_trained(self):
        return self.is_trained_

    @property
    def test_indices(self):
        return self.test_indices_

    @property
    def train_indices(self):
        return self.train_indices_

    @property
    def validation_indices(self):
        return self.validation_indices_

    def _get_user_attributes(self):
        # returns all the self attributes defined in a model class, eg, self.is_trained_
        attributes = inspect.getmembers(self, lambda a: not (inspect.isroutine(a)))
        attributes = [
            a for a in attributes if not (a[0].startswith("__") and a[0].endswith("__"))
        ]
        attributes = [a for a in attributes if not a[0].startswith("_abc_")]
        return attributes

    def _get_init_params(self, locals):
        # returns the model init signiture with associated passed in values
        # except the anndata objects passed in
        init = self.__init__
        sig = inspect.signature(init)
        init_params = [p for p in sig.parameters]
        user_params = {p: locals[p] for p in locals if p in init_params}
        user_params = {
            k: v for (k, v) in user_params.items() if not isinstance(v, AnnData)
        }
        return user_params

    def save(self, dir_path, overwrite=False):
        # get all the user attributes
        user_attributes = self._get_user_attributes()
        # only save the public attributes with _ at the very end
        user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}
        # save the model state dict and the trainer state dict only
        if not os.path.exists(dir_path) or overwrite:
            os.makedirs(dir_path, exist_ok=overwrite)
        else:
            raise ValueError(
                "{} already exists. Please provide an unexisting directory for saving.".format(
                    dir_path
                )
            )
        torch.save(self.model.state_dict(), os.path.join(dir_path, "model_params.pt"))
        with open(os.path.join(dir_path, "attr.pkl"), "wb") as f:
            pickle.dump(user_attributes, f)

    @classmethod
    def load(cls, adata: AnnData, dir_path, use_cuda=False):
        model_path = os.path.join(dir_path, "model_params.pt")
        setup_dict_path = os.path.join(dir_path, "attr.pkl")
        with open(setup_dict_path, "rb") as handle:
            attr_dict = pickle.load(handle)
        if "init_params_" not in attr_dict.keys():
            raise ValueError(
                "No init_params_ were saved by the model. Check out the developers guide if creating custom models."
            )
        # get the parameters for the class init signiture
        init_params = attr_dict.pop("init_params_")
        # grab all the parameters execept for kwargs (is a dict)
        non_kwargs = {k: v for k, v in init_params.items() if not isinstance(v, dict)}
        # expand out kwargs
        kwargs = {k: v for k, v in init_params.items() if isinstance(v, dict)}
        kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        model = cls(adata, **non_kwargs, **kwargs)
        for attr, val in attr_dict.items():
            setattr(model, attr, val)
        use_cuda = use_cuda and torch.cuda.is_available()

        if use_cuda:
            model.model.load_state_dict(torch.load(model_path))
            model.model.cuda()
        else:
            device = torch.device("cpu")
            model.model.load_state_dict(torch.load(model_path, map_location=device))
        model.model.eval()
        model._validate_anndata(adata)
        return model

    def __repr__(
        self,
    ):
        summary_string = self._model_summary_string + "\nTraining status: {}".format(
            "Trained" if self.is_trained_ else "Not Trained"
        )
        return summary_string
