import logging
import pandas as pd
import torch
from anndata import AnnData
import numpy as np
from functools import partial
from typing import Optional, Sequence, Union

from scvi._compat import Literal
from scvi.modules import PEAKVAE
from scvi.model._utils import (
    _get_batch_code_from_category,
    _get_var_names_from_setup_anndata,
)
from .base import BaseModelClass, VAEMixin

from scipy.sparse import csr_matrix, vstack
from scvi.dataloaders import ScviDataLoader
from scvi.lightning import VAETask
from scvi.utils import DifferentialComputation


logger = logging.getLogger(__name__)


class PEAKVI(VAEMixin, BaseModelClass):
    """
    PeakVI.

    Parameters
    ----------
    adata
        AnnData object that has been registered with scvi
    n_hidden
        Number of nodes per hidden layer. If `None`, defaults to square root
        of number of regions.
    n_latent
        Dimensionality of the latent space. If `None`, defaults to square root
        of `n_hidden`.
    n_layers
        Number of hidden layers used for encoder NN
    dropout_rate
        Dropout rate for neural networks
    model_depth
        Model sequencing depth / library size or not.
    region_factors
        Include region-specific factors in the model
    latent_distribution
        One of

        * ``'normal'`` - Normal distribution
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
        n_layers: int = 2,
        dropout_rate: float = 0.1,
        model_depth: bool = True,
        region_factors: bool = True,
        latent_distribution: Literal["normal", "ln"] = "normal",
        use_gpu: bool = True,
        **model_kwargs,
    ):
        super(PEAKVI, self).__init__(adata, use_gpu=use_gpu)
        self.model = PEAKVAE(
            n_input=self.summary_stats["n_vars"],
            n_batch=self.summary_stats["n_batch"],
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            latent_distribution=latent_distribution,
            model_depth=model_depth,
            region_factors=region_factors,
            **model_kwargs,
        )
        self._model_summary_string = (
            "PeakVI Model with params: \nn_hidden: {}, n_latent: {}, n_layers: {}, dropout_rate: "
            "{}, latent_distribution: {}"
        ).format(
            n_hidden,
            n_latent,
            n_layers,
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
        max_epochs=2000,
        train_size=0.9,
        validation_size=None,
        lr=1e-4,
        n_steps_kl_warmup=None,
        n_epochs_kl_warmup=50,
        **kwargs,
    ):
        task_kwargs = dict(
            weight_decay=1e-3,
            lr=lr,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
            n_steps_kl_warmup=n_steps_kl_warmup,
        )
        super().train(
            max_epochs=max_epochs,
            train_size=train_size,
            validation_size=validation_size,
            early_stopping=True,
            early_stopping_monitor="reconstruction_loss_validation",
            early_stopping_patience=100,
            vae_task_kwargs=task_kwargs,
            **kwargs,
        )
        self.is_trained_ = True

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

    def get_region_factors(self):
        if self.region_factors is None:
            raise RuntimeError("region factors were not included in this model")
        return torch.sigmoid(self.region_factors).detach().numpy()

    @torch.no_grad()
    def get_imputed_values(
        self,
        adata: Optional[AnnData] = None,
        indices: Sequence[int] = None,
        transform_batch: Optional[Union[str, int]] = None,
        use_z_mean: bool = True,
        threshold: Optional[float] = None,
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
            Should the model sample from the latent space or use the distribution means.
        threshold
            If provided, the matrix is thresholded and a binary accessibility
            matrix is returned instead. This is recommended if the matrix is very big
            and needs to be saved to file.
        batch_size
            Minibatch size for data loading into model

        """
        adata = self._validate_anndata(adata)
        post = self._make_scvi_dl(adata=adata, indices=indices, batch_size=batch_size)
        transform_batch = _get_batch_code_from_category(adata, transform_batch)

        if threshold is not None:
            assert 0 <= threshold <= 1

        # TODO: enable depth correction
        # TODO: remove technical effects? Zero-out some parts of the s tensor.
        # TODO: user-defined batch size
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
            p = generative_outputs["p"]
            if threshold:
                p = csr_matrix((p >= threshold).numpy())
            imputed.append(p)

        if threshold:  # imputed is a list of csr_matrix objects
            imputed = vstack(imputed, format="csr")
        else:  # imputed is a list of tensors
            imputed = torch.cat(imputed).numpy()
        return imputed

    def differential_accessibility(
        self,
        groupby,
        group1,
        group2=None,
        adata=None,
        mode="vanilla",
        use_permutation=False,
    ):
        adata = self._validate_anndata(adata)
        cell_idx1 = adata.obs[groupby] == group1
        if group2 is None:
            cell_idx2 = ~cell_idx1
        else:
            cell_idx2 = adata.obs[groupby] == group2

        model_fn = partial(self.get_imputed_values, use_z_mean=False)
        dc = DifferentialComputation(model_fn, adata)
        all_info = dc.get_bayes_factors(
            cell_idx1, cell_idx2, mode=mode, use_permutation=use_permutation
        )

        region_names = _get_var_names_from_setup_anndata(adata)
        res = pd.DataFrame(all_info, index=region_names)
        sort_key = "proba_de" if mode == "change" else "bayes_factor"
        res = res.sort_values(by=sort_key, ascending=False)
        return res
