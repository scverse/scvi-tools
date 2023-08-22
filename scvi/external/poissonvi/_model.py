import logging
from typing import Dict, Iterable, List, Literal, Optional, Sequence, Union

import numpy as np
from anndata import AnnData
import pandas as pd
import torch

from scvi import REGISTRY_KEYS
from scvi._types import Number
from scvi.data import AnnDataManager
from scvi.data.fields import (
    BaseAnnDataField,
    CategoricalJointObsField,
    CategoricalObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
    ObsmField,
    StringUnsField,
)
from scvi.model import PEAKVI
from scvi.model._utils import _init_library_size
from scvi.model.base import (
    RNASeqMixin,
    UnsupervisedTrainingMixin,
    VAEMixin,
    ArchesMixin,
    BaseModelClass,
)
from scvi.module import VAE
from scvi.train._callbacks import SaveBestState
from scvi.utils import setup_anndata_dsp
from scvi.utils._docstrings import devices_dsp

logger = logging.getLogger(__name__)


class POISSONVI(
    RNASeqMixin,
    VAEMixin,
    ArchesMixin,
    UnsupervisedTrainingMixin,
    BaseModelClass,
):
    """
    Peak Variational Inference using a Poisson distribution :cite:p:`Martens22`.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.POISSONVI.setup_anndata`.
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
    deeply_inject_covariates
        Whether to deeply inject covariates into all layers of the decoder. If False (default),
        covariates will only be included in the input layer.
    **model_kwargs
        Keyword args for :class:`~scvi.module.POISSONVAE`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.POISSINVI.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.model.POISSONVI(adata)
    >>> vae.train()

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/PeakVI`
    """

    # TODO: change tutorial link
    # TODO: change citation
    # TODO: change model_kwargs docstring
    _module_cls = VAE

    def __init__(
        self,
        adata: AnnData,
        n_hidden: Optional[int] = None,
        n_latent: Optional[int] = None,
        n_layers: int = 2,
        dropout_rate: float = 0.1,
        latent_distribution: Literal["normal", "ln"] = "normal",
        **model_kwargs,
    ):
        super().__init__(adata)

        n_cats_per_cov = (
            self.adata_manager.get_state_registry(
                REGISTRY_KEYS.CAT_COVS_KEY
            ).n_cats_per_key
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )
        n_batch = self.summary_stats.n_batch
        use_size_factor_key = (
            REGISTRY_KEYS.SIZE_FACTOR_KEY in self.adata_manager.data_registry
        )
        library_log_means, library_log_vars = None, None
        if not use_size_factor_key is None:
            library_log_means, library_log_vars = _init_library_size(
                self.adata_manager, n_batch
            )

        # to be consitent with PEAKVI architecture
        n_hidden = (
            int(np.sqrt(self.summary_stats.n_vars)) if n_hidden is None else n_hidden
        )
        n_latent = int(np.sqrt(n_hidden)) if n_latent is None else n_latent

        self.module = self._module_cls(
            n_input=self.summary_stats.n_vars,
            n_batch=n_batch,
            n_labels=self.summary_stats.n_labels,
            n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
            n_cats_per_cov=n_cats_per_cov,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion="gene",  # not needed here
            gene_likelihood="poisson",  # fixed value for now, but we could think of also allowing nb
            latent_distribution=latent_distribution,
            use_size_factor_key=use_size_factor_key,
            library_log_means=library_log_means,
            library_log_vars=library_log_vars,
            use_batch_norm="none",  # to be consitent with PEAKVI architecture
            use_layer_norm="both",  # to be consitent with PEAKVI architecture
            extra_encoder_kwargs={
                "activation_fn": torch.nn.LeakyReLU
            },  # to be consitent with PEAKVI architecture
            extra_decoder_kwargs={
                "activation_fn": torch.nn.LeakyReLU
            },  # to be consitent with PEAKVI architecture
            **model_kwargs,
        )
        self._model_summary_string = (
            "PoissonVI Model with the following params: \nn_hidden: {}, n_latent: {}, n_layers: {}, dropout_rate: "
            "{}, peak_likelihood: {}, latent_distribution: {}"
        ).format(
            n_hidden,
            n_latent,
            n_layers,
            dropout_rate,
            "poisson",
            latent_distribution,
        )
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: Optional[str] = None,
        batch_key: Optional[str] = None,
        labels_key: Optional[str] = None,
        size_factor_key: Optional[str] = None,
        categorical_covariate_keys: Optional[List[str]] = None,
        continuous_covariate_keys: Optional[List[str]] = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        %(param_batch_key)s
        %(param_labels_key)s
        %(param_size_factor_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        """
        # TODO: where should we check that we are using fragment counts?
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
            NumericalObsField(
                REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False
            ),
            CategoricalJointObsField(
                REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys
            ),
            NumericalJointObsField(
                REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys
            ),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    # to be consitent with PEAKVI training
    @devices_dsp.dedent
    def train(
        self,
        max_epochs: int = 500,
        lr: float = 1e-4,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float = 0.9,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        batch_size: int = 128,
        weight_decay: float = 1e-3,
        eps: float = 1e-08,
        early_stopping: bool = True,
        early_stopping_patience: int = 50,
        save_best: bool = True,
        check_val_every_n_epoch: int | None = None,
        n_steps_kl_warmup: int | None = None,
        n_epochs_kl_warmup: int | None = 50,
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        **kwargs,
    ):
        """Trains the model using amortized variational inference.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset.
        lr
            Learning rate for optimization.
        %(param_accelerator)s
        %(param_devices)s
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        shuffle_set_split
            Whether to shuffle indices before splitting. If `False`, the val, train, and test set are split in the
            sequential order of the data according to `validation_size` and `train_size` percentages.
        batch_size
            Minibatch size to use during training.
        weight_decay
            weight decay regularization term for optimization
        eps
            Optimizer eps
        early_stopping
            Whether to perform early stopping with respect to the validation set.
        early_stopping_patience
            How many epochs to wait for improvement before early stopping
        save_best
            Save the best model state with respect to the validation loss (default), or use the final
            state in the training procedure
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
        datasplitter_kwargs
            Additional keyword arguments passed into :class:`~scvi.dataloaders.DataSplitter`.
        plan_kwargs
            Keyword args for :class:`~scvi.train.TrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        update_dict = {
            "lr": lr,
            "weight_decay": weight_decay,
            "eps": eps,
            "n_epochs_kl_warmup": n_epochs_kl_warmup,
            "n_steps_kl_warmup": n_steps_kl_warmup,
            "optimizer": "AdamW",
        }
        if plan_kwargs is not None:
            plan_kwargs.update(update_dict)
        else:
            plan_kwargs = update_dict
        if save_best:
            if "callbacks" not in kwargs.keys():
                kwargs["callbacks"] = []
            kwargs["callbacks"].append(
                SaveBestState(monitor="reconstruction_loss_validation")
            )

        super().train(
            max_epochs=max_epochs,
            train_size=train_size,
            accelerator=accelerator,
            devices=devices,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            early_stopping=early_stopping,
            early_stopping_monitor="reconstruction_loss_validation",
            early_stopping_patience=early_stopping_patience,
            datasplitter_kwargs=datasplitter_kwargs,
            plan_kwargs=plan_kwargs,
            check_val_every_n_epoch=check_val_every_n_epoch,
            batch_size=batch_size,
            **kwargs,
        )

    @torch.inference_mode()
    def get_accessibility_estimates(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] = None,
        transform_batch: str | int | None = None,
        region_list: Sequence[str] | None = None,
        library_size: Union[float, Literal["latent"]] = 1,
        normalize_regions: bool = False,
        n_samples: int = 1,
        n_samples_overall: int = None,
        weights: Optional[Literal["uniform", "importance"]] = None,
        batch_size: int = 128,
        return_mean: bool = True,
        return_numpy: bool = False,
        **importance_weighting_kwargs,
    ) -> pd.DataFrame | np.ndarray | csr_matrix:
        """Returns the normalized accessibility matrix.

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
        region_list
            Return frequencies of accessibility for a subset of regions.
            This can save memory when working with large datasets and few regions are
            of interest.
        library_size
            Scale the accessibility frequencies to a common library size.
            This allows accessibility counts to be interpreted on a common scale of relevant
            magnitude. If set to `"latent"`, use the latent library size.
        normalize_regions
            Whether to reintroduce region factors to scale the normalized accessibility. This makes
            the estimates closer to the input, but removes the region-level bias correction. False by
            default.
        n_samples
            Number of posterior samples to use for estimation.
        n_samples_overall
            Number of posterior samples to use for estimation. Overrides `n_samples`.
        weights
            Weights to use for sampling. If `None`, defaults to `"uniform"`.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`. DataFrame includes
            region names as columns. If either `n_samples=1` or `return_mean=True`, defaults to `False`.
            Otherwise, it defaults to `True`.
        importance_weighting_kwargs
            Keyword arguments passed into :meth:`~scvi.model.base.RNASeqMixin._get_importance_weights`.

         Returns
        -------
        If `n_samples` is provided and `return_mean` is False,
        this method returns a 3d tensor of shape (n_samples, n_cells, n_regions).
        If `n_samples` is provided and `return_mean` is True, it returns a 2d tensor
        of shape (n_cells, n_regions).
        In this case, return type is :class:`~pandas.DataFrame` unless `return_numpy` is True.
        Otherwise, the method expects `n_samples_overall` to be provided and returns a 2d tensor
        of shape (n_samples_overall, n_regions).
        """
        if not normalize_regions:
            region_factors = self.module.decoder.px_scale_decoder[-2].bias
            # set region_factors (bias) to 0
            self.module.decoder.px_scale_decoder[-2].bias = torch.nn.Parameter(
                torch.zeros_like(region_factors)
            )
        accs = super().get_normalized_expression(
            adata=adata,
            indices=indices,
            transform_batch=transform_batch,
            gene_list=region_list,
            library_size=library_size,
            n_samples=n_samples,
            n_samples_overall=n_samples_overall,
            weights=weights,
            batch_size=batch_size,
            return_mean=return_mean,
            return_numpy=return_numpy,
            **importance_weighting_kwargs,
        )
        if not normalize_regions:
            # reset region_factors (bias)
            self.module.decoder.px_scale_decoder[-2].bias = torch.nn.Parameter(
                region_factors
            )
        return accs

    @torch.inference_mode()
    def get_normalized_expression(
        self,
    ):
        # Refer to function get_accessibility_estimates
        print(
            "get_normalized_expression is not implemented for POISSONVI, please use get_accessibility_estimates"
        )
        return None

    @torch.inference_mode()
    def differential_accessibility(self, groupby: str, group1: str, group2: str):
        """
        Parameters
        ----------
        groupby
            The key of the observation grouping to consider.
        group1
            Name of first group
        group2
            Name of second group

        Returns
        -------
        Differential accessibility scores
        """
        return self.module.differential_accessibility(groupby, group1, group2)

    @torch.inference_mode()
    def differential_expression(
        self,
    ):
        # Refer to function differential_accessibility
        print(
            "differential_expression is not implemented for POISSONVI, please use differential_accessibility"
        )
        return None
