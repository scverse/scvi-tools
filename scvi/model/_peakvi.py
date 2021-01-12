import logging
from functools import partial
from typing import Iterable, Optional, Sequence, Union

import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from scipy.sparse import csr_matrix, vstack

from scvi._compat import Literal
from scvi._docs import doc_differential_expression
from scvi._utils import _doc_params
from scvi.dataloaders import ScviDataLoader
from scvi.lightning import VAETask
from scvi.modules import PEAKVAE
from scvi.model._utils import (
    _get_batch_code_from_category,
    _get_var_names_from_setup_anndata,
    scrna_raw_counts_properties,
)

from .base import BaseModelClass, VAEMixin
from .base._utils import _de_core

logger = logging.getLogger(__name__)


class PEAKVI(VAEMixin, BaseModelClass):
    """
    PeakVI.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    n_hidden
        Number of nodes per hidden layer. If `None`, defaults to square root
        of number of regions.
    n_latent
        Dimensionality of the latent space. If `None`, defaults to square root
        of `n_hidden`.
    n_layers_encoder
        Number of hidden layers used for encoder NN.
    n_layers_decoder
        Number of hidden layers used for decoder NN.
    dropout_rate
        Dropout rate for neural networks
    model_depth
        Model sequencing depth / library size (default: True)
    region_factors
        Include region-specific factors in the model (default: True)
    latent_distribution
        One of

        * ``'normal'`` - Normal distribution (Default)
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.dataset.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.models.PeakVI(adata)
    >>> vae.train()
    """

    def __init__(
        self,
        adata: AnnData,
        n_hidden: Optional[int] = 128,
        n_latent: Optional[int] = None,
        n_layers_encoder: int = 2,
        n_layers_decoder: int = 2,
        dropout_rate: float = 0.1,
        model_depth: bool = True,
        region_factors: bool = True,
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        latent_distribution: Literal["normal", "ln"] = "normal",
        use_gpu: bool = True,
        **model_kwargs,
    ):
        super(PEAKVI, self).__init__(adata, use_gpu=use_gpu)

        n_cats_per_cov = (
            self.scvi_setup_dict_["extra_categoricals"]["n_cats_per_key"]
            if "extra_categoricals" in self.scvi_setup_dict_
            else []
        )

        self.model = PEAKVAE(
            n_input_regions=self.summary_stats["n_vars"],
            n_batch=self.summary_stats["n_batch"],
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers_encoder=n_layers_encoder,
            n_layers_decoder=n_layers_decoder,
            n_continuous_cov=self.summary_stats["n_continuous_covs"],
            n_cats_per_cov=n_cats_per_cov,
            dropout_rate=dropout_rate,
            model_depth=model_depth,
            region_factors=region_factors,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            latent_distribution=latent_distribution,
            **model_kwargs,
        )
        self._model_summary_string = (
            "PeakVI Model with params: \nn_hidden: {}, n_latent: {}, n_layers_encoder: {}, "
            "n_layers_decoder: {} , dropout_rate: {}, latent_distribution: {}"
        ).format(
            self.model.n_hidden,
            self.model.n_latent,
            n_layers_encoder,
            n_layers_decoder,
            dropout_rate,
            latent_distribution,
        )
        self.n_latent = n_latent
        self.init_params_ = self._get_init_params(locals())

    @property
    def _data_loader_cls(self):
        return ScviDataLoader

    @property
    def _task_class(self):
        return VAETask

    def train(
        self,
        max_epochs: int = 500,
        lr: float = 1e-3,
        use_gpu: Optional[bool] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 256,
        weight_decay: float = 1e-3,
        early_stopping: bool = True,
        check_val_every_n_epoch: Optional[int] = None,
        n_steps_kl_warmup: Union[int, None] = None,
        n_epochs_kl_warmup: Union[int, None] = 50,
        vae_task_kwargs: Optional[dict] = None,
        **kwargs,
    ):
        """
        Trains the model using amortized variational inference.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset.
        lr
            Learning rate for optimization.
        use_gpu
            If `True`, use the GPU if available.
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        weight_decay
            weight decay regularization term for optimization
        early_stopping
            Whether to perform early stopping with respect to the validation set.
        check_val_every_n_epoch
            Check val every n train epochs. By default, val is not checked, unless `early_stopping` is `True`.
            If so, val is checked every epoch.
        n_steps_kl_warmup
            Number of training steps (minibatches) to scale weight on KL divergences from 0 to 1.
            Only activated when `n_epochs_kl_warmup` is set to None. If `None`, defaults
            to `floor(0.75 * adata.n_obs)`.
        n_epochs_kl_warmup
            Number of epochs to scale weight on KL divergences from 0 to 1.
            Overrides `n_steps_kl_warmup` when both are not `None`.
        vae_task_kwargs
            Keyword args for :class:`~scvi.lightning.VAETask`. Keyword arguments passed to
            `train()` will overwrite values present in `vae_task_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.lightning.Trainer`.
        """
        update_dict = dict(
            lr=lr,
            weight_decay=weight_decay,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
            n_steps_kl_warmup=n_steps_kl_warmup,
        )
        if vae_task_kwargs is not None:
            vae_task_kwargs.update(update_dict)
        else:
            vae_task_kwargs = update_dict
        super().train(
            max_epochs=max_epochs,
            train_size=train_size,
            use_gpu=use_gpu,
            validation_size=validation_size,
            early_stopping=True,
            early_stopping_monitor="reconstruction_loss_validation",
            early_stopping_patience=50,
            vae_task_kwargs=vae_task_kwargs,
            **kwargs,
        )

    @torch.no_grad()
    def get_library_size_factors(
        self,
        adata: Optional[AnnData] = None,
        indices: Sequence[int] = None,
        batch_size: int = 128,
    ):
        adata = self._validate_anndata(adata)
        scdl = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)

        library_sizes = []
        for tensors in scdl:
            inference_inputs = self.model._get_inference_input(tensors)
            outputs = self.model.inference(**inference_inputs)
            library_sizes.append(outputs["d"].cpu())

        return torch.cat(library_sizes).numpy().squeeze()

    @torch.no_grad()
    def get_region_factors(self):
        if self.region_factors is None:
            raise RuntimeError("region factors were not included in this model")
        return torch.sigmoid(self.region_factors).cpu().numpy()

    @torch.no_grad()
    def get_imputed_values(
        self,
        adata: Optional[AnnData] = None,
        indices: Sequence[int] = None,
        transform_batch: Optional[Union[str, int]] = None,
        use_z_mean: bool = True,
        threshold: Optional[float] = None,
        scale_cells: bool = False,
        scale_regions: bool = False,
        batch_size: int = 128,
    ) -> Union[np.ndarray, csr_matrix]:
        """
        Impute the full accessibility matrix.

        Returns a matrix of accessibility probabilities for each cell and genomic region in the input
        (for return matrix A, A[i,j] is the probability that region j is accessible in cell i).

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
        use_z_mean
            If True (default), use the distribution mean. Otherwise, sample from the distribution.
        threshold
            If provided, the matrix is thresholded and a sparse binary accessibility
            matrix is returned instead. This is recommended if the matrix is very big
            and needs to be saved to file.
        normalize_cells
            Whether to reintroduce library size factors to scale the normalized probabilities.
            This makes the estimates closer to the input, but removes the library size correction.
            False by default.
        normalize_regions
            Whether to reintroduce region factors to scale the normalized probabilities. This makes
            the estimates closer to the input, but removes the region-level bias correction. False by 
            default.
        batch_size
            Minibatch size for data loading into model

        """
        adata = self._validate_anndata(adata)
        post = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)
        transform_batch = _get_batch_code_from_category(adata, transform_batch)

        if threshold is not None:
            assert 0 <= threshold <= 1

        imputed = []
        for tensors in post:
            # TODO implement iteration over multiple batches like in RNAMixin
            get_generative_input_kwargs = dict(transform_batch=transform_batch[0])
            generative_kwargs = dict(use_z_mean=use_z_mean)
            _, generative_outputs = self.model.forward(
                tensors=tensors,
                get_generative_input_kwargs=get_generative_input_kwargs,
                generative_kwargs=generative_kwargs,
                compute_loss=False,
            )
            p = generative_outputs["p"].cpu()
            if threshold:
                p = csr_matrix((p >= threshold).numpy())
            imputed.append(p)

        if threshold:  # imputed is a list of csr_matrix objects
            imputed = vstack(imputed, format="csr")
        else:  # imputed is a list of tensors
            imputed = torch.cat(imputed).numpy()

        if scale_cells:
            imputed = imputed * self.get_library_size_factors(adata, indices, batch_size)
        if scale_regions:
            imputed = imputed * self.get_region_factors()
        
        return imputed

    @_doc_params(
        doc_differential_expression=doc_differential_expression,
    )
    def differential_accessibility(
        self,
        adata: Optional[AnnData] = None,
        groupby: Optional[str] = None,
        group1: Optional[Iterable[str]] = None,
        group2: Optional[str] = None,
        idx1: Optional[Union[Sequence[int], Sequence[bool]]] = None,
        idx2: Optional[Union[Sequence[int], Sequence[bool]]] = None,
        mode: Literal["vanilla", "change"] = "change",
        delta: float = 0.1,
        batch_size: Optional[int] = None,
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: Optional[Iterable[str]] = None,
        batchid2: Optional[Iterable[str]] = None,
        fdr_target: float = 0.05,
        two_sided: bool = True,
        **kwargs,
    ) -> pd.DataFrame:
        r"""
        A unified method for differential accessibility analysis.

        Implements `"vanilla"` DE [Lopez18]_ and `"change"` mode DE [Boyeau19]_.

        Parameters
        ----------
        {doc_differential_expression}
        two_sided
            Whether to perform a two-sided test, or a one-sided test. 
        **kwargs
            Keyword args for :func:`scvi.utils.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential accessibility DataFrame.
        """
        adata = self._validate_anndata(adata)
        col_names = _get_var_names_from_setup_anndata(adata)
        model_fn = partial(
            self.get_imputed_values, use_z_mean=False, batch_size=batch_size
        )

        # TODO check if change_fn in kwargs and raise error if so
        def change_fn(a, b):
            return a - b

        if two_sided:
            def m1_domain_fn(samples):
                return np.abs(samples) >= delta
        else:
            def m1_domain_fn(samples):
                return samples >= delta

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
            change_fn=change_fn,
            m1_domain_fn=m1_domain_fn,
            **kwargs,
        )

        return result
