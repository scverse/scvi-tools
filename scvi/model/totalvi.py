import logging
from functools import partial
from typing import Dict, Iterable, Optional, Sequence, Tuple, Union, TypeVar
from collections.abc import Iterable as IterableClass

import numpy as np
import pandas as pd
import torch
from anndata import AnnData

from scvi import _CONSTANTS
from scvi._compat import Literal
from scvi._docs import doc_differential_expression
from scvi._utils import _doc_params
from scvi.core.data_loaders import TotalDataLoader
from scvi.core.distributions import NegativeBinomial, NegativeBinomialMixture
from scvi.core.models import ArchesMixin, BaseModelClass, RNASeqMixin, VAEMixin
from scvi.core.models._utils import _de_core
from scvi.core.modules import TOTALVAE
from scvi.core.trainers import TotalTrainer
from scvi.data import get_from_registry
from scvi.data._utils import _check_nonnegative_integers
from scvi.model._utils import (
    _get_var_names_from_setup_anndata,
    _get_batch_code_from_category,
    cite_seq_raw_counts_properties,
)

logger = logging.getLogger(__name__)
Number = TypeVar("Number", int, float)


class TOTALVI(RNASeqMixin, VAEMixin, ArchesMixin, BaseModelClass):
    """
    total Variational Inference [GayosoSteier20]_.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    n_latent
        Dimensionality of the latent space.
    gene_dispersion
        One of the following:

        * ``'gene'`` - genes_dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - genes_dispersion can differ between different batches
        * ``'gene-label'`` - genes_dispersion can differ between different labels
    protein_dispersion
        One of the following:

        * ``'protein'`` - protein_dispersion parameter is constant per protein across cells
        * ``'protein-batch'`` - protein_dispersion can differ between different batches NOT TESTED
        * ``'protein-label'`` - protein_dispersion can differ between different labels NOT TESTED
    gene_likelihood
        One of:

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
    latent_distribution
        One of:

        * ``'normal'`` - Normal distribution
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)
    empirical_protein_background_prior
        Set the initialization of protein background prior empirically. This option fits a GMM for each of
        100 cells per batch and averages the distributions. Note that even with this option set to `True`,
        this only initializes a parameter that is learned during inference. If `False`, randomly initializes.
    use_cuda
        Use the GPU or not.
    **model_kwargs
        Keyword args for :class:`~scvi.core.modules.TOTALVAE`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.data.setup_anndata(adata, batch_key="batch", protein_expression_obsm_key="protein_expression")
    >>> vae = scvi.model.TOTALVI(adata)
    >>> vae.train()
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
        latent_distribution: Literal["normal", "ln"] = "normal",
        empirical_protein_background_prior: bool = True,
        use_cuda: bool = True,
        **model_kwargs,
    ):
        super(TOTALVI, self).__init__(adata, use_cuda=use_cuda)
        if "totalvi_batch_mask" in self.scvi_setup_dict_.keys():
            batch_mask = self.scvi_setup_dict_["totalvi_batch_mask"]
        else:
            batch_mask = None
        if empirical_protein_background_prior:
            prior_mean, prior_scale = _get_totalvi_protein_priors(adata)
        else:
            prior_mean, prior_scale = None, None
        self.model = TOTALVAE(
            n_input_genes=self.summary_stats["n_vars"],
            n_input_proteins=self.summary_stats["n_proteins"],
            n_batch=self.summary_stats["n_batch"],
            n_latent=n_latent,
            gene_dispersion=gene_dispersion,
            protein_dispersion=protein_dispersion,
            gene_likelihood=gene_likelihood,
            latent_distribution=latent_distribution,
            protein_batch_mask=batch_mask,
            protein_background_prior_mean=prior_mean,
            protein_background_prior_scale=prior_scale,
            **model_kwargs,
        )
        self._model_summary_string = (
            "TotalVI Model with the following params: \nn_latent: {}, "
            "gene_dispersion: {}, protein_dispersion: {}, gene_likelihood: {}, latent_distribution: {}"
        ).format(
            n_latent,
            gene_dispersion,
            protein_dispersion,
            gene_likelihood,
            latent_distribution,
        )
        self.init_params_ = self._get_init_params(locals())

    def train(
        self,
        n_epochs: int = 400,
        train_size: float = 0.9,
        test_size: Optional[float] = None,
        lr: float = 4e-3,
        n_epochs_kl_warmup: Optional[int] = None,
        n_iter_kl_warmup: Union[Literal["auto"], int] = "auto",
        batch_size: int = 256,
        frequency: Optional[int] = None,
        train_fun_kwargs: dict = {},
        **kwargs,
    ):
        """
        Train the model.

        Parameters
        ----------
        n_epochs
            Number of passes through the dataset.
        train_size
            Size of training set in the range [0.0, 1.0].
        test_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + test_size < 1`, the remaining cells belong to a validation set.
        lr
            Learning rate for optimization.
        n_epochs_kl_warmup
            Number of passes through dataset for scaling term on KL divergence to go from 0 to 1.
        n_iter_kl_warmup
            Number of minibatches for scaling term on KL divergence to go from 0 to 1.
            To use, set to not `None` and set `n_epochs_kl_warmup` to `None`.
        batch_size
            Minibatch size to use during training.
        frequency
            Frequency with which metrics are computed on the data for train/test/val sets.
        train_fun_kwargs
            Keyword args for the train method of :class:`~scvi.core.trainers.TotalTrainer`.
        **kwargs
            Other keyword args for :class:`~scvi.core.trainers.TotalTrainer`.
        """
        train_fun_kwargs = dict(train_fun_kwargs)
        if "totalvi_batch_mask" in self.scvi_setup_dict_.keys():
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
            frequency=frequency,
            batch_size=batch_size,
            use_adversarial_loss=imputation,
            use_cuda=self.use_cuda,
            **kwargs,
        )
        # for autotune
        if "n_epochs" not in train_fun_kwargs:
            train_fun_kwargs["n_epochs"] = n_epochs
        if "lr" not in train_fun_kwargs:
            train_fun_kwargs["lr"] = lr

        logger.info("Training for {} epochs.".format(n_epochs))

        self.trainer.train(**train_fun_kwargs)
        self.is_trained_ = True
        self.train_indices_ = self.trainer.train_set.indices
        self.test_indices_ = self.trainer.test_set.indices
        self.validation_indices_ = self.trainer.validation_set.indices
        self.history_ = self.trainer.history

    @torch.no_grad()
    def get_reconstruction_error(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        mode: Literal["total", "gene", "protein"] = "total",
        batch_size: Optional[int] = None,
    ):
        r"""
        Return the reconstruction error for the data.

        This is typically written as :math:`p(x, y \mid z)`, the likelihood term given one posterior sample.
        Note, this is not the negative likelihood, higher is better.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        mode
            Compute for genes, proteins, or both.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        """
        adata = self._validate_anndata(adata)
        post = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)

        return -post.reconstruction_error(mode=mode)

    @torch.no_grad()
    def get_latent_representation(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        give_mean: bool = True,
        mc_samples: int = 5000,
        batch_size: Optional[int] = None,
    ) -> np.ndarray:
        """
        Return the latent representation for each cell.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        give_mean
            Give mean of distribution or sample from it
        mc_samples
            For distributions with no closed-form mean (e.g., `logistic normal`), how many Monte Carlo
            samples to take for computing mean.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

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
        if self.is_trained_ is False:
            raise RuntimeError("Please train the model first.")

        adata = self._validate_anndata(adata)
        post = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)
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
    def get_latent_library_size(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        give_mean: bool = True,
        batch_size: Optional[int] = None,
    ) -> np.ndarray:
        r"""
        Returns the latent library size for each cell.

        This is denoted as :math:`\ell_n` in the totalVI paper.

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
        post = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)
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
        transform_batch: Optional[Sequence[Union[Number, str]]] = None,
        gene_list: Optional[Sequence[str]] = None,
        protein_list: Optional[Sequence[str]] = None,
        library_size: Optional[Union[float, Literal["latent"]]] = 1,
        n_samples: int = 1,
        sample_protein_mixing: bool = False,
        scale_protein: bool = False,
        include_protein_background: bool = False,
        batch_size: Optional[int] = None,
        return_mean: bool = True,
        return_numpy: Optional[bool] = None,
    ) -> Tuple[Union[np.ndarray, pd.DataFrame], Union[np.ndarray, pd.DataFrame]]:
        r"""
        Returns the normalized gene expression and protein expression.

        This is denoted as :math:`\rho_n` in the totalVI paper for genes, and TODO
        for proteins, :math:`(1-\pi_{nt})\alpha_{nt}\beta_{nt}`.

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
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
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
        adata = self._validate_anndata(adata)
        post = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = _get_var_names_from_setup_anndata(adata)
            gene_mask = [True if gene in gene_list else False for gene in all_genes]
        if protein_list is None:
            protein_mask = slice(None)
        else:
            all_proteins = self.scvi_setup_dict_["protein_names"]
            protein_mask = [True if p in protein_list else False for p in all_proteins]
        if indices is None:
            indices = np.arange(adata.n_obs)

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                logger.warning(
                    "return_numpy must be True if n_samples > 1 and return_mean is False, returning np.ndarray"
                )
            return_numpy = True

        if not isinstance(transform_batch, IterableClass):
            transform_batch = [transform_batch]

        transform_batch = _get_batch_code_from_category(adata, transform_batch)

        scale_list_gene = []
        scale_list_pro = []

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
                if library_size == "latent":
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
                index=adata.obs_names[indices],
            )
            pro_df = pd.DataFrame(
                scale_list_pro,
                columns=self.scvi_setup_dict_["protein_names"][protein_mask],
                index=adata.obs_names[indices],
            )

            return gene_df, pro_df
        else:
            return scale_list_gene, scale_list_pro

    @torch.no_grad()
    def get_protein_foreground_probability(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        transform_batch: Optional[Sequence[Union[Number, str]]] = None,
        protein_list: Optional[Sequence[str]] = None,
        n_samples: int = 1,
        batch_size: Optional[int] = None,
        return_mean: bool = True,
        return_numpy: Optional[bool] = None,
    ):
        r"""
        Returns the foreground probability for proteins.

        This is denoted as :math:`(1 - \pi_{nt})` in the totalVI paper.

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

            - None, then real observed batch is used
            - int, then batch transform_batch is used
            - List[int], then average over batches in list
        protein_list
            Return protein expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
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
        - **foreground_probability** - probability foreground for each protein

        If `n_samples` > 1 and `return_mean` is False, then the shape is `(samples, cells, genes)`.
        Otherwise, shape is `(cells, genes)`. In this case, return type is :class:`~pandas.DataFrame` unless `return_numpy` is True.
        """
        adata = self._validate_anndata(adata)
        post = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)

        if protein_list is None:
            protein_mask = slice(None)
        else:
            all_proteins = self.scvi_setup_dict_["protein_names"]
            protein_mask = [True if p in protein_list else False for p in all_proteins]

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                logger.warning(
                    "return_numpy must be True if n_samples > 1 and return_mean is False, returning np.ndarray"
                )
            return_numpy = True
        if indices is None:
            indices = np.arange(adata.n_obs)

        py_mixings = []
        if not isinstance(transform_batch, IterableClass):
            transform_batch = [transform_batch]

        transform_batch = _get_batch_code_from_category(adata, transform_batch)
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            y = tensors[_CONSTANTS.PROTEIN_EXP_KEY]
            batch_index = tensors[_CONSTANTS.BATCH_KEY]
            label = tensors[_CONSTANTS.LABELS_KEY]
            py_mixing = torch.zeros_like(y[..., protein_mask])
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
            pro_names = self.scvi_setup_dict_["protein_names"]
            foreground_prob = pd.DataFrame(
                1 - py_mixings,
                columns=pro_names[protein_mask],
                index=adata.obs_names[indices],
            )
            return foreground_prob

    def _expression_for_de(
        self,
        adata=None,
        indices=None,
        transform_batch: Optional[Sequence[Union[Number, str]]] = None,
        scale_protein=False,
        batch_size: Optional[int] = None,
        sample_protein_mixing=False,
        include_protein_background=False,
        protein_prior_count=0.5,
    ):
        rna, protein = self.get_normalized_expression(
            adata=adata,
            indices=indices,
            transform_batch=transform_batch,
            return_numpy=True,
            n_samples=1,
            batch_size=batch_size,
            scale_protein=scale_protein,
            sample_protein_mixing=sample_protein_mixing,
            include_protein_background=include_protein_background,
        )
        protein += protein_prior_count

        joint = np.concatenate([rna, protein], axis=1)
        return joint

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
        protein_prior_count: float = 0.1,
        scale_protein: bool = False,
        sample_protein_mixing: bool = False,
        include_protein_background: bool = False,
        **kwargs,
    ) -> pd.DataFrame:
        r"""
        A unified method for differential expression analysis.

        Implements `"vanilla"` DE [Lopez18]_ and `"change"` mode DE [Boyeau19]_.

        Parameters
        ----------
        {doc_differential_expression}
        protein_prior_count
            Prior count added to protein expression before LFC computation
        scale_protein
            Force protein values to sum to one in every single cell (post-hoc normalization)
        sample_protein_mixing
            Sample the protein mixture component, i.e., use the parameter to sample a Bernoulli
            that determines if expression is from foreground/background.
        include_protein_background
            Include the protein background component as part of the protein expression
        **kwargs
            Keyword args for :func:`scvi.core.utils.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential expression DataFrame.
        """
        adata = self._validate_anndata(adata)
        model_fn = partial(
            self._expression_for_de,
            scale_protein=scale_protein,
            sample_protein_mixing=sample_protein_mixing,
            include_protein_background=include_protein_background,
            protein_prior_count=protein_prior_count,
            batch_size=batch_size,
        )
        col_names = np.concatenate(
            [
                np.asarray(_get_var_names_from_setup_anndata(adata)),
                self.scvi_setup_dict_["protein_names"],
            ]
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
            cite_seq_raw_counts_properties,
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
        batch_size: Optional[int] = None,
        gene_list: Optional[Sequence[str]] = None,
        protein_list: Optional[Sequence[str]] = None,
    ) -> np.ndarray:
        r"""
        Generate observation samples from the posterior predictive distribution.

        The posterior predictive distribution is written as :math:`p(\hat{x}, \hat{y} \mid x, y)`.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        n_samples
            Number of required samples for each cell
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        gene_list
            Names of genes of interest
        protein_list
            Names of proteins of interest

        Returns
        -------
        x_new : :class:`~numpy.ndarray`
            tensor with shape (n_cells, n_genes, n_samples)
        """
        if self.model.gene_likelihood not in ["nb"]:
            raise ValueError("Invalid gene_likelihood")

        adata = self._validate_anndata(adata)
        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = _get_var_names_from_setup_anndata(adata)
            gene_mask = [True if gene in gene_list else False for gene in all_genes]
        if protein_list is None:
            protein_mask = slice(None)
        else:
            all_proteins = self.scvi_setup_dict_["protein_names"]
            protein_mask = [True if p in protein_list else False for p in all_proteins]

        post = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)

        scdl_list = []
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

            rna_dist = NegativeBinomial(mu=px_["rate"], theta=px_["r"])
            protein_dist = NegativeBinomialMixture(
                mu1=py_["rate_back"],
                mu2=py_["rate_fore"],
                theta1=py_["r"],
                mixture_logits=py_["mixing"],
            )

            rna_sample = rna_dist.sample().cpu()[..., gene_mask]
            protein_sample = protein_dist.sample().cpu()[..., protein_mask]
            data = torch.cat([rna_sample, protein_sample], dim=-1).numpy()

            scdl_list += [data]
            if n_samples > 1:
                scdl_list[-1] = np.transpose(scdl_list[-1], (1, 2, 0))
        scdl_list = np.concatenate(scdl_list, axis=0)

        return scdl_list

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
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            indices of `adata` to use
        n_samples
            How may samples per cell
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        rna_size_factor
            size factor for RNA prior to sampling gamma distribution
        transform_batch
            int of which batch to condition on for all cells
        """
        adata = self._validate_anndata(adata)
        post = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)

        scdl_list = []
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
            scdl_list += [data]

            scdl_list[-1] = np.transpose(scdl_list[-1], (1, 2, 0))

        return np.concatenate(scdl_list, axis=0)

    @torch.no_grad()
    def get_feature_correlation_matrix(
        self,
        adata=None,
        indices=None,
        n_samples: int = 10,
        batch_size: int = 64,
        rna_size_factor: int = 1000,
        transform_batch: Optional[Sequence[Union[Number, str]]] = None,
        correlation_type: Literal["spearman", "pearson"] = "spearman",
        log_transform: bool = False,
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
            size factor for RNA prior to sampling gamma distribution
        transform_batch
            Batches to condition on.
            If transform_batch is:

            - None, then real observed batch is used
            - int, then batch transform_batch is used
            - list of int, then values are averaged over provided batches.
        correlation_type
            One of "pearson", "spearman".
        log_transform
            Whether to log transform denoised values prior to correlation calculation.

        Returns
        -------
        Gene-protein-gene-protein correlation matrix
        """
        from scipy.stats import spearmanr

        adata = self._validate_anndata(adata)

        if not isinstance(transform_batch, IterableClass):
            transform_batch = [transform_batch]

        transform_batch = _get_batch_code_from_category(adata, transform_batch)

        corr_mats = []
        for b in transform_batch:
            denoised_data = self._get_denoised_samples(
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
                flattened[:, : self.n_genes] = np.log(
                    flattened[:, : self.n_genes] + 1e-8
                )
                flattened[:, self.n_genes :] = np.log1p(flattened[:, self.n_genes :])
            if correlation_type == "pearson":
                corr_matrix = np.corrcoef(flattened, rowvar=False)
            else:
                corr_matrix, _ = spearmanr(flattened, axis=0)
            corr_mats.append(corr_matrix)

        corr_matrix = np.mean(np.stack(corr_mats), axis=0)
        var_names = _get_var_names_from_setup_anndata(adata)
        names = np.concatenate(
            [np.asarray(var_names), self.scvi_setup_dict_["protein_names"]]
        )
        return pd.DataFrame(corr_matrix, index=names, columns=names)

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
        Estimates for the parameters of the likelihood :math:`p(x, y \mid z)`.

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
        raise NotImplementedError

    def _validate_anndata(
        self, adata: Optional[AnnData] = None, copy_if_view: bool = True
    ):
        adata = super()._validate_anndata(adata, copy_if_view)
        error_msg = "Number of {} in anndata different from when setup_anndata was run. Please rerun setup_anndata."
        if _CONSTANTS.PROTEIN_EXP_KEY in adata.uns["_scvi"]["data_registry"].keys():
            pro_exp = get_from_registry(adata, _CONSTANTS.PROTEIN_EXP_KEY)
            if self.summary_stats["n_proteins"] != pro_exp.shape[1]:
                raise ValueError(error_msg.format("proteins"))
            is_nonneg_int = _check_nonnegative_integers(pro_exp)
            if not is_nonneg_int:
                logger.warning(
                    "Make sure the registered protein expression in anndata contains unnormalized count data."
                )
        else:
            raise ValueError("No protein data found, please setup or transfer anndata")

        return adata

    @property
    def _trainer_class(self):
        return TotalTrainer

    @property
    def _scvi_dl_class(self):
        return TotalDataLoader


def _get_totalvi_protein_priors(adata, n_cells=100):
    """Compute an empirical prior for protein background."""
    import warnings
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.mixture import GaussianMixture

    warnings.filterwarnings("error")

    batch = get_from_registry(adata, _CONSTANTS.BATCH_KEY).ravel()
    cats = adata.uns["_scvi"]["categorical_mappings"]["_scvi_batch"]["mapping"]
    codes = np.arange(len(cats))

    batch_avg_mus, batch_avg_scales = [], []
    for b in np.unique(codes):
        # can happen during online updates
        # the values of these batches will not be used
        num_in_batch = np.sum(batch == b)
        if num_in_batch == 0:
            batch_avg_mus.append(0)
            batch_avg_scales.append(1)
            continue
        pro_exp = get_from_registry(adata, _CONSTANTS.PROTEIN_EXP_KEY)[batch == b]

        # for missing batches, put dummy values -- scarches case, will be replaced anyway
        if pro_exp.shape[0] == 0:
            batch_avg_mus.append(0.0)
            batch_avg_scales.append(0.05)

        cells = np.random.choice(np.arange(pro_exp.shape[0]), size=n_cells)
        if isinstance(pro_exp, pd.DataFrame):
            pro_exp = pro_exp.to_numpy()
        pro_exp = pro_exp[cells]
        gmm = GaussianMixture(n_components=2)
        mus, scales = [], []
        # fit per cell GMM
        for c in pro_exp:
            try:
                gmm.fit(np.log1p(c.reshape(-1, 1)))
            # when cell is all 0
            except ConvergenceWarning:
                mus.append(0)
                scales.append(0.05)
                continue

            means = gmm.means_.ravel()
            sorted_fg_bg = np.argsort(means)
            mu = means[sorted_fg_bg].ravel()[0]
            covariances = gmm.covariances_[sorted_fg_bg].ravel()[0]
            scale = np.sqrt(covariances)
            mus.append(mu)
            scales.append(scale)

        # average distribution over cells
        batch_avg_mu = np.mean(mus)
        batch_avg_scale = np.sqrt(np.sum(np.square(scales)) / (n_cells ** 2))

        batch_avg_mus.append(batch_avg_mu)
        batch_avg_scales.append(batch_avg_scale)

    # repeat prior for each protein
    batch_avg_mus = np.array(batch_avg_mus, dtype=np.float32).reshape(1, -1)
    batch_avg_scales = np.array(batch_avg_scales, dtype=np.float32).reshape(1, -1)
    batch_avg_mus = np.tile(batch_avg_mus, (pro_exp.shape[1], 1))
    batch_avg_scales = np.tile(batch_avg_scales, (pro_exp.shape[1], 1))

    warnings.resetwarnings()

    return batch_avg_mus, batch_avg_scales
