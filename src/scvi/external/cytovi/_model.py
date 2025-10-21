from __future__ import annotations

import logging
import warnings
from functools import partial
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import rich
import torch
import torch.distributions as dist
from anndata import AnnData
from tqdm import tqdm

from scvi import settings
from scvi.data import AnnDataManager
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LayerField,
    NumericalJointObsField,
)
from scvi.dataloaders import DataSplitter
from scvi.distributions._utils import DistributionConcatenator
from scvi.model._utils import _get_batch_code_from_category, scrna_raw_counts_properties
from scvi.model.base import (
    ArchesMixin,
    BaseModelClass,
    RNASeqMixin,
    UnsupervisedTrainingMixin,
    VAEMixin,
)
from scvi.model.base._de_core import _de_core
from scvi.train import AdversarialTrainingPlan, TrainRunner
from scvi.utils import de_dsp, setup_anndata_dsp
from scvi.utils._docstrings import devices_dsp

from ._constants import CYTOVI_DEFAULT_REP, CYTOVI_REGISTRY_KEYS
from ._module import CytoVAE
from ._utils import (
    clip_lfc_factory,
    encode_categories,
    get_balanced_sample_indices,
    get_n_latent_heuristic,
    impute_cats_with_neighbors,
    impute_expr_with_neighbors,
    log_median,
    validate_expression_range,
    validate_layer_key,
    validate_marker,
    validate_obs_keys,
    validate_obsm_keys,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence
    from typing import Literal

    from scvi._types import Number

logger = logging.getLogger(__name__)


class CYTOVI(
    RNASeqMixin,
    VAEMixin,
    ArchesMixin,
    UnsupervisedTrainingMixin,
    BaseModelClass,
):
    """Variational inference for cytometry (CytoVI).

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.external.CYTOVI.setup_anndata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space. If None, will be set using a heuristic based on
        number of input features.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    dropout_rate
        Dropout rate for neural networks.
    protein_likelihood
        Likelihood function used for modeling protein expression. One of:

        * ``'normal'`` - Normal distribution
        * ``'beta'`` - Beta distribution (requires expression in (0, 1))

    latent_distribution
        Distribution of the latent space. One of:

        * ``'normal'`` - Normal distribution
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)

    encode_backbone_only
        If True, only encode backbone markers (i.e., those present in all samples). This is
        required when analyzing overlapping panels with missing values.
    encoder_marker_list
        Optional list of markers to use for encoding. Must be a subset of backbone markers
          if `encode_backbone_only` is True.
    prior_mixture
        If True, uses a mixture of Gaussians as a prior in the latent space (MoG prior).
    prior_mixture_k
        Number of mixture components in the MoG prior. Defaults to `n_latent` if None.
    **model_kwargs
        Keyword arguments passed to :class:`~scvi.external.cytovi.CytoVAE`.

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.external.CYTOVI.setup_anndata(adata, batch_key="batch")
    >>> model = scvi.external.CYTOVI(adata)
    >>> model.train()
    >>> adata.obsm["X_CytoVI"] = model.get_latent_representation()
    >>> adata.layers["imputed"] = model.get_normalized_expression()

    Notes
    -----
    When analyzing overlapping cytometry panels (i.e., samples with partially shared markers),
    CytoVI uses the shared "backbone" markers for encoding and reconstructs the full set.
    An adversarial classifier loss can be used to encourage batch-invariance in the latent space.
    If the data includes missing values, ensure that `nan_layer` is correctly registered using
    :meth:`~scvi.external.CYTOVI.setup_anndata`. This is handled automatically when using
     scvi.external.cytovi.merge_batches().

    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/cytometry/CytoVI_batch_correction_tutorial`
    2. :doc:`/tutorials/notebooks/cytometry/CytoVI_advanced_tutorial`
    """

    _module_cls = CytoVAE
    _training_plan_cls = AdversarialTrainingPlan
    _data_splitter_cls = DataSplitter
    _train_runner_cls = TrainRunner

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int | None = None,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        protein_likelihood: Literal["normal", "beta"] = "normal",
        latent_distribution: Literal["normal", "ln"] = "normal",
        encode_backbone_only: bool | None = None,
        encoder_marker_list: list | None = None,
        prior_mixture: bool | None = True,
        prior_mixture_k: int | None = None,
        **model_kwargs,
    ):
        super().__init__(adata)

        n_cats_per_cov = (
            self.adata_manager.get_state_registry(CYTOVI_REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
            if CYTOVI_REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )
        n_batch = self.summary_stats.n_batch
        all_markers = adata.var_names

        if encoder_marker_list is not None:
            validate_marker(adata, encoder_marker_list)
            encoder_marker_mask = all_markers.isin(encoder_marker_list)
        else:
            encoder_marker_mask = None

        if CYTOVI_REGISTRY_KEYS.PROTEIN_NAN_MASK in self.adata_manager.data_registry:
            nan_layer = self.adata_manager.get_from_registry("nan_layer")

            backbone_markers = list(all_markers[~np.any(nan_layer == 0, axis=0)])
            self.backbone_markers = backbone_markers
            self.nan_imputation = True
            self.backbone_marker_mask = all_markers.isin(backbone_markers)
            backbone_str = ", ".join(backbone_markers)
            self._use_adversarial_classifier = True

            if encode_backbone_only is None and encoder_marker_list is None:
                encode_backbone_only = True
            elif encode_backbone_only is False:
                raise NotImplementedError(
                    "When analyzing overlapping panels, only encoding of the backbone markers is "
                    "currently supported."
                )

            if encoder_marker_mask is not None:
                enc_marker_intersection = [
                    marker in backbone_markers for marker in encoder_marker_list
                ]
                probl_markers = [
                    marker
                    for marker, intersection in zip(
                        encoder_marker_list, enc_marker_intersection, strict=False
                    )
                    if not intersection
                ]
                probl_markers_str = ", ".join(probl_markers)
                if not all(enc_marker_intersection):
                    raise ValueError(
                        f"{probl_markers_str} are in 'encoder_marker_list' but not in backbone "
                        "marker list. When analyzing overlapping panels, only encoding of the "
                        "backbone markers is currently supported."
                    )

            if encode_backbone_only:
                encoder_marker_mask = self.backbone_marker_mask

        else:
            self.backbone_markers = None
            self.backbone_marker_mask = None
            self.nan_imputation = False
            self._use_adversarial_classifier = False

        if n_latent is None:
            if encoder_marker_mask is not None:
                n_vars_encoded = encoder_marker_mask.sum()
            else:
                n_vars_encoded = self.summary_stats.n_vars
            n_latent = get_n_latent_heuristic(n_vars_encoded)

        if prior_mixture_k is None:
            prior_mixture_k = n_latent

        if protein_likelihood == "beta":
            expr = self.adata_manager.get_from_registry("X")
            corr_range = validate_expression_range(expr, 0, 1)
            if not corr_range:
                raise ValueError(
                    "Protein expression must be in the range (0, 1) for beta likelihood. "
                    "Perform scaling or choose other likelihood."
                )

        self.sample_key = self.adata_manager.get_state_registry(
            CYTOVI_REGISTRY_KEYS.SAMPLE_KEY
        ).original_key

        self.batch_key = self.adata_manager.get_state_registry(
            CYTOVI_REGISTRY_KEYS.BATCH_KEY
        ).original_key

        self._model_summary_string = (  # noqa: UP032
            "CytoVI Model with the following params: \nn_hidden: {}, n_latent: {}, n_layers: {}, "
            "dropout_rate: {}, \nprotein_likelihood: {}, latent_distribution: {}, \nMoG prior: {},"
            " n_labels {}, n_proteins: {}, \nImpute missing markers: {}"
        ).format(
            n_hidden,
            n_latent,
            n_layers,
            dropout_rate,
            protein_likelihood,
            latent_distribution,
            prior_mixture,
            self.summary_stats.n_labels,
            self.summary_stats.n_vars,
            self.nan_imputation,
        )

        if self.nan_imputation is True:
            self._model_summary_string += f", \nBackbone markers: {backbone_str}"

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
            protein_likelihood=protein_likelihood,
            latent_distribution=latent_distribution,
            encoder_marker_mask=encoder_marker_mask,
            prior_mixture=prior_mixture,
            prior_mixture_k=prior_mixture_k,
            **model_kwargs,
        )

        self.init_params_ = self._get_init_params(locals())

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        batch_key: str | None = None,
        labels_key: str | None = None,
        sample_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        nan_layer: str | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        layer
            if not `None`, uses this as the key in `adata.layers` for the transformed
            protein expression data.
        %(param_batch_key)s
        %(param_labels_key)s
        %(param_sample_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        nan_layer
            Optional layer key containing binary NaN feature mask to handle overlapping
            antibody panels.
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(CYTOVI_REGISTRY_KEYS.X_KEY, layer, is_count_data=False),
            CategoricalObsField(CYTOVI_REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(CYTOVI_REGISTRY_KEYS.LABELS_KEY, labels_key),
            CategoricalObsField(CYTOVI_REGISTRY_KEYS.SAMPLE_KEY, sample_key),
            CategoricalJointObsField(
                CYTOVI_REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys
            ),
            NumericalJointObsField(CYTOVI_REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
        ]

        if nan_layer is None and "_nan_mask" in adata.layers:
            msg = (
                "Found nan_layer in adata. Will register nan_layer for missing marker imputation."
            )
            warnings.warn(msg, UserWarning, stacklevel=settings.warnings_stacklevel)
            nan_layer = "_nan_mask"

        if nan_layer is not None:
            anndata_fields.append(LayerField(CYTOVI_REGISTRY_KEYS.PROTEIN_NAN_MASK, nan_layer))

        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    def __repr__(
        self,
    ):
        summary_string = self._model_summary_string
        summary_string += "\nTraining status: {}".format(
            "Trained" if self.is_trained_ else "Not Trained"
        )
        rich.print(summary_string)
        return ""

    @devices_dsp.dedent
    def train(
        self,
        max_epochs: int | None = 1000,
        lr: float = 1e-3,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float = 0.9,
        validation_size: float | None = None,
        batch_size: int = 128,
        early_stopping: bool = True,
        check_val_every_n_epoch: int | None = None,
        n_steps_kl_warmup: int | None = None,
        n_epochs_kl_warmup: int | None = 400,
        adversarial_classifier: bool | None = None,
        plan_kwargs: dict | None = None,
        early_stopping_patience: int | None = 30,
        **kwargs,
    ):
        """
        Trains the model using amortized variational inference.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset, by default 1000.
        lr
            Learning rate for optimization, by default 1e-3.
        accelerator
            Accelerator to use for training, by default "auto".
        devices
            Devices to use for training, by default "auto".
        train_size
            Size of the training set in the range [0.0, 1.0], by default 0.9.
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training, by default 128.
        early_stopping
            Whether to perform early stopping with respect to the validation set, by default True.
        check_val_every_n_epoch
            Check validation set every n train epochs. By default, the validation set is not
            checked, unless `early_stopping` is `True` or `reduce_lr_on_plateau` is `True`.
            If either of the latter conditions are met, the validation set is checked
            every epoch.
        n_steps_kl_warmup
            Number of training steps (minibatches) to scale weight on KL divergences from 0 to 1.
            Only activated when `n_epochs_kl_warmup` is set to None. If `None`, defaults
            to `floor(0.75 * adata.n_obs)`.
        n_epochs_kl_warmup
            Number of epochs to scale weight on KL divergences from 0 to 1.
            Overrides `n_steps_kl_warmup` when both are not `None`, by default 400.
        adversarial_classifier
            Whether to use an adversarial classifier in the latent space. This helps mixing when
            there are missing proteins in any of the batches. Defaults to `True` if missing
            proteins are detected.
        plan_kwargs
            Keyword arguments for the `AdversarialTrainingPlan` class. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        early_stopping_patience
            Number of epochs to wait before early stopping, by default 30.
        **kwargs
            Other keyword arguments for the `Trainer` class.

        Returns
        -------
        runner : object
            The runner object used for training.
        """
        if adversarial_classifier is None:
            adversarial_classifier = self._use_adversarial_classifier

        update_dict = {
            "lr": lr,
            "adversarial_classifier": adversarial_classifier,
            "n_epochs_kl_warmup": n_epochs_kl_warmup,
            "n_steps_kl_warmup": n_steps_kl_warmup,
        }
        if plan_kwargs is not None:
            plan_kwargs.update(update_dict)
        else:
            plan_kwargs = update_dict

        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else {}

        data_splitter = self._data_splitter_cls(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
        )
        training_plan = self._training_plan_cls(self.module, **plan_kwargs)
        runner = self._train_runner_cls(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=devices,
            early_stopping=early_stopping,
            check_val_every_n_epoch=check_val_every_n_epoch,
            early_stopping_patience=early_stopping_patience,
            **kwargs,
        )
        return runner()

    @torch.inference_mode()
    def posterior_predictive_sample(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        n_samples: int = 1,
        protein_list: Sequence[str] | None = None,
        batch_size: int | None = None,
    ):
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
        protein_list
            Names of proteins of interest.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        x_new : :py:class:`torch.Tensor`
            Tensor with shape (n_cells, n_proteins, n_samples)
        """
        if self.module.protein_likelihood not in ["beta", "normal"]:
            raise ValueError("Invalid protein_likelihood.")

        adata = self._validate_anndata(adata)

        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        if indices is None:
            indices = np.arange(adata.n_obs)

        if protein_list is None:
            protein_mask = slice(None)
        else:
            all_proteins = adata.var_names
            protein_mask = [True if protein in protein_list else False for protein in all_proteins]

        x_new = []
        for tensors in scdl:
            samples = self.module.sample(
                tensors,
                n_samples=n_samples,
            )
            if protein_list is not None:
                samples = samples[:, protein_mask, ...]
            x_new.append(samples)

        x_new = torch.cat(x_new)  # Shape (n_cells, n_genes, n_samples)

        return x_new.numpy()

    @torch.inference_mode()
    def get_normalized_expression(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        transform_batch: Sequence[Number | str] | None = "all",
        protein_list: Sequence[str] | None = None,
        n_samples: int = 1,
        n_samples_overall: int = None,
        weights: Literal["uniform", "importance"] | None = None,
        batch_size: int | None = None,
        return_mean: bool = True,
        return_numpy: bool | None = None,
        nan_warning: bool | None = True,
        **importance_weighting_kwargs,
    ) -> np.ndarray | pd.DataFrame:
        r"""
        Returns the normalized (decoded) protein expression.

        The model's reconstructed (normalized) expression is written as
        :math:\hat{x} = p(x \mid z).

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
            - 'all', then the mean across batches is used
            - None, then real observed batch is used.
            - int, then batch transform_batch is used.
            This behaviour affects only proteins that are detected across multiple batches.
            Unobserved proteins are decoded in the batch(es), in which they were measured.
        protein_list
            Return frequencies of expression for a subset of protein.
            This can save memory when working with large datasets and few proteins are
            of interest.
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
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`. DataFrame
            includes gene names as columns. If either `n_samples=1` or `return_mean=True`,
            defaults to `False`. Otherwise, it defaults to `True`.
        nan_warning
            Whether to show a warning if missing proteins are detected between batches.
        **importance_weighting_kwargs
            Additional keyword arguments for importance weighting.

        Returns
        -------
        If `n_samples` > 1 and `return_mean` is False, then the shape is
        `(samples, cells, genes)`. Otherwise, shape is `(cells, genes)`. In this case,
        return type is :class:`~pandas.DataFrame` unless `return_numpy` is True.
        """
        adata = self._validate_anndata(adata)
        all_batches = list(np.unique(self.adata_manager.get_from_registry("batch")))

        if self.nan_imputation is True:
            if nan_warning is True:
                msg = "detected missing proteins between batches - will impute missing markers"
                warnings.warn(msg, UserWarning, stacklevel=settings.warnings_stacklevel)

        if indices is None:
            indices = np.arange(adata.n_obs)
        if n_samples_overall is not None:
            assert n_samples == 1  # default value
            n_samples = n_samples_overall // len(indices) + 1
        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        if protein_list is None:
            protein_mask = slice(None)
        else:
            protein_mask = [
                True if protein in protein_list else False for protein in adata.var_names
            ]

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                msg = "return_numpy must be True if n_samples > 1 and return_mean is False, "
                "returning np.ndarray"
                warnings.warn(msg, UserWarning, stacklevel=settings.warnings_stacklevel)
            return_numpy = True

        if transform_batch == "all":
            transform_batch = all_batches
        else:
            transform_batch = _get_batch_code_from_category(
                self.get_anndata_manager(adata, required=True), transform_batch
            )

        store_distributions = weights == "importance"
        if store_distributions and len(transform_batch) > 1:
            raise NotImplementedError(
                "Importance weights cannot be computed when expression levels are averaged "
                "across batches."
            )

        exprs = []
        zs = []
        qz_store = DistributionConcatenator()
        px_store = DistributionConcatenator()
        for tensors in scdl:
            per_batch_exprs = []
            for batch in transform_batch:
                generative_kwargs = self._get_transform_batch_gen_kwargs(batch)
                inference_kwargs = {"n_samples": n_samples}
                inference_outputs, generative_outputs = self.module.forward(
                    tensors=tensors,
                    inference_kwargs=inference_kwargs,
                    generative_kwargs=generative_kwargs,
                    compute_loss=False,
                )

                output = generative_outputs["px"].mean
                output = output.cpu().numpy()

                output = output[..., protein_mask]

                per_batch_exprs.append(output)
                if store_distributions:
                    qz_store.store_distribution(inference_outputs["qz"])
                    px_store.store_distribution(generative_outputs["px"])
            zs.append(inference_outputs["z"].cpu())
            per_batch_exprs = np.stack(per_batch_exprs)

            exprs += [per_batch_exprs.mean(axis=0)]

        if n_samples > 1:
            # The -2 axis correspond to cells.
            exprs = np.concatenate(exprs, axis=-2)
            zs = torch.concat(zs, dim=-2)
        else:
            exprs = np.concatenate(exprs, axis=0)
            zs = torch.concat(zs, dim=0)

        if n_samples_overall is not None:
            # Converts the 3d tensor to a 2d tensor
            exprs = exprs.reshape(-1, exprs.shape[-1])
            n_samples_ = exprs.shape[0]
            if (weights is None) or weights == "uniform":
                p = None
            else:
                qz = qz_store.get_concatenated_distributions(axis=0)
                x_axis = 0 if n_samples == 1 else 1
                px = px_store.get_concatenated_distributions(axis=x_axis)
                p = self.get_importance_weights(
                    adata,
                    indices,
                    qz=qz,
                    px=px,
                    zs=zs,
                    **importance_weighting_kwargs,
                )
            ind_ = np.random.choice(n_samples_, n_samples_overall, p=p, replace=True)
            exprs = exprs[ind_]
        elif n_samples > 1 and return_mean:
            exprs = exprs.mean(axis=0)

        if return_numpy is None or return_numpy is False:
            return pd.DataFrame(
                exprs,
                columns=adata.var_names[protein_mask],
                index=adata.obs_names[indices],
            )
        else:
            return exprs

    @de_dsp.dedent
    def differential_expression(
        self,
        adata: AnnData | None = None,
        groupby: str | None = None,
        group1: list[str] | None = None,
        group2: list[str] | None = None,
        idx1: list[int] | list[bool] | str | None = None,
        idx2: list[int] | list[bool] | str | None = None,
        mode: Literal["vanilla", "change"] = "change",
        test_mode: Literal["two", "three"] = "two",
        delta: float = 0.25,
        batch_size: int | None = None,
        all_stats: bool = False,
        batch_correction: bool = False,
        batchid1: list[str] | None = None,
        batchid2: list[str] | None = None,
        fdr_target: float = 0.05,
        silent: bool = False,
        weights: Literal["uniform", "importance"] | None = "uniform",
        filter_outlier_cells: bool = False,
        lfc_clipping: bool = True,
        clipping_range: tuple = (0, 1),
        balance_samples: bool | None = None,
        importance_weighting_kwargs: dict | None = None,
        **kwargs,
    ) -> pd.DataFrame:
        r"""A unified method for differential expression analysis.

        Implements ``'vanilla'`` DE :cite:p:`Lopez18` and ``'change'`` mode DE :cite:p:`Boyeau19`.

        Parameters
        ----------
        adata
            Annotated data matrix, by default None
        groupby
            Key in `adata.obs` containing the groups, by default None
        group1
            List of group names for group 1, by default None
        group2
            List of group names for group 2, by default None
        idx1
            Indices or boolean mask for group 1, by default None
        idx2
            Indices or boolean mask for group 2, by default None
        mode
            Differential expression mode, by default "change"
        test_mode
            Test mode for differential expression, by default "two"
        delta
            Minimum fold change for differential expression, by default 0.25
        batch_size
            Batch size for computation, by default None
        all_stats
            Whether to compute all statistics, by default False
        batch_correction
            Whether to perform batch correction, by default False
        batchid1
            List of batch IDs for group 1, by default None
        batchid2
            List of batch IDs for group 2, by default None
        fdr_target
            Target false discovery rate, by default 0.05
        silent
            Whether to suppress progress bar and messages, by default False
        weights
            Weights to use for sampling, by default "uniform"
        filter_outlier_cells
            Whether to filter outlier cells, by default False
        lfc_clipping
            Whether to clip log fold change values, by default True
        clipping_range
            Range for clipping log fold change values, by default (0, 1)
        balance_samples
            Whether to subsample equal amount of cells per sample based on `sample_key`.
            If `None`, defaults to `True` if more than one sample is present in the data.
        importance_weighting_kwargs
            Keyword arguments for importance weighting, by default None
        **kwargs
            Additional keyword arguments for differential expression computation

        Returns
        -------
        Differential expression DataFrame.
        """
        adata = self._validate_anndata(adata)
        col_names = adata.var_names
        importance_weighting_kwargs = importance_weighting_kwargs or {}
        model_fn = partial(
            self.get_normalized_expression,
            return_numpy=True,
            n_samples=1,
            batch_size=batch_size,
            weights=weights,
            nan_warning=False,
            **importance_weighting_kwargs,
        )
        representation_fn = self.get_latent_representation if filter_outlier_cells else None

        if lfc_clipping is True:
            eps = 1e-6
            clip_min = clipping_range[0] + eps
            clip_max = clipping_range[1] - eps

            expr = self.adata_manager.get_from_registry("X")
            corr_range = validate_expression_range(expr, clipping_range[0], clipping_range[1])
            if not corr_range:
                msg = "Protein expression exceeds clipping range, which can lead to poor "
                "DE results. Please adjust clipping range to data range."
                warnings.warn(msg, UserWarning, stacklevel=settings.warnings_stacklevel)
            change_fn_clp = clip_lfc_factory(clip_min, clip_max)

            if kwargs is None:
                kwargs = {}
                kwargs["change_fn"] = change_fn_clp
                kwargs["test_mode"] = test_mode
            else:
                kwargs["change_fn"] = change_fn_clp
                kwargs["test_mode"] = test_mode

        if self.registry_["setup_args"]["sample_key"] and balance_samples is not False:
            subset_idx = get_balanced_sample_indices(adata, self.sample_key)
            kwargs["subset_idx"] = subset_idx

        result = _de_core(
            self.get_anndata_manager(adata, required=True),
            model_fn,
            representation_fn,
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
            silent,
            **kwargs,
        )

        return result

    def get_aggregated_posterior(
        self,
        adata: AnnData = None,
        sample: int | str = None,
        indices: Sequence[int] = None,
        batch_size: int = None,
        dof: float | None = 3.0,
    ) -> dist.Distribution:
        """Compute the aggregated posterior over the ``z`` latent representations.

        Parameters
        ----------
        adata
            AnnData object to use. Defaults to the AnnData object used to initialize the model.
        sample
            Name or index of the sample to filter on. If ``None``, uses all cells.
        indices
            Indices of cells to use.
        batch_size
            Batch size to use for computing the latent representation.
        dof
            Degrees of freedom for the Student's t-distribution components. If ``None``,
            components are Normal.

        Returns
        -------
        A mixture distribution of the aggregated posterior.
        """
        self._check_if_trained(warn=False)
        adata = self._validate_anndata(adata)

        if indices is None:
            indices = np.arange(self.adata.n_obs)
        if sample is not None:
            indices = np.intersect1d(
                np.array(indices), np.where(adata.obs[self.sample_key] == sample)[0]
            )

        dataloader = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        qu_loc, qu_scale = self.get_latent_representation(
            batch_size=batch_size, return_dist=True, dataloader=dataloader, give_mean=True
        )

        qu_loc = torch.tensor(qu_loc, device=self.device).T
        qu_scale = torch.tensor(qu_scale, device=self.device).T

        if dof is None:
            components = dist.Normal(qu_loc, qu_scale)
        else:
            components = dist.StudentT(dof, qu_loc, qu_scale)
        return dist.MixtureSameFamily(
            dist.Categorical(logits=torch.ones(qu_loc.shape[1], device=qu_loc.device)), components
        )

    def get_sample_logprobs(
        self,
        adata: AnnData | None = None,
        batch_size: int = 128,
        downsample_cells: int | None = None,
        dof: float | None = None,
    ) -> pd.DataFrame:
        """Compute the differential abundance log probabilities for each sample.

        Computes the logarithm of the ratio of the probabilities of each sample conditioned on the
        estimated aggregate posterior distribution of each cell.

        Parameters
        ----------
        adata
            The data object to compute the differential abundance for.
        batch_size
            Minibatch size for computing the differential abundance.
        downsample_cells
            Number of cells to subset to before computing the differential abundance.
        dof
            Degrees of freedom for the Student's t-distribution components for aggregated
            posterior. If ``None``, components are Normal.

        Returns
        -------
        DataFrame of shape (n_cells, n_samples) containing the log probabilities
        for each cell across samples. The rows correspond to cell names from `adata.obs_names`,
        and the columns correspond to unique sample identifiers.
        """
        adata = self._validate_anndata(adata)

        zs = self.get_latent_representation(
            batch_size=batch_size, return_dist=False, give_mean=True
        )

        unique_samples = adata.obs[self.sample_key].unique()
        dataloader = torch.utils.data.DataLoader(zs, batch_size=batch_size)
        log_probs = []
        for sample_name in tqdm(unique_samples):
            indices = np.where(adata.obs[self.sample_key] == sample_name)[0]
            if downsample_cells is not None and downsample_cells < indices.shape[0]:
                indices = np.random.choice(indices, downsample_cells, replace=False)

            ap = self.get_aggregated_posterior(adata=adata, indices=indices, dof=dof)
            log_probs_ = []
            for z_rep in dataloader:
                z_rep = z_rep.to(self.device)
                log_probs_.append(ap.log_prob(z_rep).sum(-1, keepdims=True))
            log_probs.append(torch.cat(log_probs_, axis=0).cpu().numpy())

        log_probs = np.concatenate(log_probs, 1)
        log_probs_df = pd.DataFrame(
            data=log_probs, index=adata.obs_names.to_numpy(), columns=unique_samples
        )
        return log_probs_df

    def differential_abundance(
        self,
        adata: AnnData | None = None,
        groupby: str | None = None,
        batch_size: int = 128,
        downsample_cells: int | None = None,
        dof: float | None = None,
        aggregation_fn: Callable[[np.ndarray], float] = log_median,
        return_log_probs: bool = False,
    ) -> pd.DataFrame:
        """
        Compute differential abundance (DA) scores across experimental conditions.

        This function estimates differential abundance by comparing aggregated
        log-sample-probabilities between groups defined by `groupby`. The aggregation
        is performed using the specified `aggregation_fn`, typically `log_median` or
        `logsumexp`, across posterior sample log-probabilities.

        Parameters
        ----------
        adata
            AnnData object to use. If `None`, defaults to the model's internal AnnData.
        groupby
            Key in `adata.obs` that contains condition or group labels.
            If not provided, returns log-probabilities per sample without aggregation.
        batch_size
            Mini-batch size for computing log-probabilities. Default: 128.
        downsample_cells
            If provided, randomly subsample this many cells before computing
            log-probabilities.
        dof
            Degrees of freedom to use in the sampling distribution, if applicable.
        aggregation_fn
            Function used to aggregate log-probabilities across samples.
            Common choices include `log_median` or `logsumexp`. Default: `log_median`.
        return_log_probs
            If `True`, skip aggregation and return per-sample log-probabilities directly.
            Default: False.

        Returns
        -------
        If `return_log_probs=True` or `groupby=None`, returns a DataFrame of
        shape `(n_cells, n_samples)` containing log-sample-probabilities.

        Otherwise, returns a DataFrame of shape `(n_cells, n_conditions)` where
        each column corresponds to the differential abundance log-ratio for one
        condition compared to all others.

        Examples
        --------
        >>> da_scores = model.differential_abundance(adata, groupby="condition")
        >>> da_scores.head()

        >>> log_probs = model.differential_abundance(adata, return_log_probs=True)
        """
        adata = self._validate_anndata(adata)

        log_probs = self.get_sample_logprobs(
            adata, batch_size=batch_size, downsample_cells=downsample_cells, dof=dof
        )

        if return_log_probs:
            return log_probs

        if groupby is None:
            warnings.warn(
                "`groupby` is not specified. Will return the log probabilities per sample "
                "without aggregation.",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )
            return log_probs
        else:
            validate_obs_keys(adata, groupby)

        md = adata.obs[[self.sample_key, groupby]].drop_duplicates()

        da_dict = {}
        for cond in set(md[groupby].values):
            is_case = (md[groupby] == cond).values
            log_probs_cond = aggregation_fn(log_probs.loc[:, is_case], 1)
            log_prop_controls = aggregation_fn(log_probs.loc[:, ~is_case], 1)
            log_ratios = log_probs_cond - log_prop_controls
            da_dict[f"DA_{cond}"] = log_ratios

        da_df = pd.DataFrame(da_dict)

        return da_df

    def impute_categories_from_reference(
        self,
        adata_reference: AnnData,
        cat_key: str,
        use_rep: str | None = None,
        n_neighbors: int = 20,
        return_uncertainty: bool = False,
    ):
        """
        Impute missing categories from a reference dataset using a shared representation.

        Parameters
        ----------
        adata_reference
            Annotated data matrix for the reference dataset. This dataset contains the categories
            to be imputed onto the query data.
        cat_key
            The key in the `.obs` attribute of `adata_reference` that specifies the categorical
            variable (e.g., cell types or clusters) to impute.
        use_rep
            The key in the `.obsm` attribute to use as the representation space (e.g., latent
            space). If `None`, the function will attempt to use a default latent
            representation `X_CytoVI`.
        n_neighbors
            The number of nearest neighbors to use for imputation. The imputation is based on
            similarity in the chosen representation space. Default: 20.
        return_uncertainty
            If `True`, the function will also return the uncertainty of the imputation.
            Default: False.

        Returns
        -------
        Array of imputed categories for the query dataset, corresponding to the categorical
        variable specified by `cat_key`.

        Notes
        -----
        This function assumes that both the query and reference datasets have a precomputed
        representation in `.obsm` (typically CytoVI latent space). If not, you must either
        provide a common representation manually or ensure that one is generated.
        """
        adata_query = self.adata

        if use_rep is None:
            ref_obsm_keys = adata_reference.obsm.keys()
            query_obsm_keys = adata_query.obsm.keys()
            shared_obsm_keys = [key for key in ref_obsm_keys if key in query_obsm_keys]

            if CYTOVI_DEFAULT_REP in shared_obsm_keys:
                use_rep = CYTOVI_DEFAULT_REP
            elif CYTOVI_DEFAULT_REP in ref_obsm_keys:
                adata_query.obsm[CYTOVI_DEFAULT_REP] = self.get_latent_representation()
                use_rep = CYTOVI_DEFAULT_REP
            else:
                raise ValueError(
                    "No shared representation found between reference and query data. Please "
                    "specify a representation to use."
                )

        # Validate input keys
        validate_obsm_keys(adata_query, use_rep)
        validate_obsm_keys(adata_reference, use_rep)
        validate_obs_keys(adata_reference, cat_key)

        # One-Hot Encode the reference categories
        cat_encoded_ref, ohe = encode_categories(adata_reference, cat_key)
        n_cats = cat_encoded_ref.shape[1]

        # Get representations
        rep_ref = adata_reference.obsm[use_rep]
        rep_query = adata_query.obsm[use_rep]

        # Impute missing categories for the query data
        imputed_query_cat_indices, uncertainty = impute_cats_with_neighbors(
            rep_query,
            rep_ref,
            cat_encoded_ref,
            n_neighbors=n_neighbors,
            compute_uncertainty=return_uncertainty,
        )

        # Convert imputed indices back to category labels
        imputed_query_cat = ohe.inverse_transform(
            np.eye(n_cats)[imputed_query_cat_indices]
        ).reshape(-1)

        if return_uncertainty:
            return imputed_query_cat, uncertainty
        else:
            return imputed_query_cat

    def impute_rna_from_reference(
        self: AnnData,
        reference_batch: str,
        adata_rna: AnnData,
        layer_key: str,
        use_rep: str | None = None,
        n_neighbors: int = 20,
        compute_uncertainty: bool = False,
        return_query_only: bool = False,
    ):
        """
        Impute RNA expression in a query dataset using a reference and a shared representation.

        Parameters
        ----------
        reference_batch
            Identifier for the reference batch in ``adata.obs[batch_key]``.
        adata_rna
            Annotated data matrix containing the expression data to impute.
        layer_key
            Key in ``adata_rna.layers`` for the reference expression data.
        use_rep
            Key in ``.obsm`` to use as the representation space (e.g., latent space).
            If ``None``, defaults to ``X_CytoVI``.
        n_neighbors
            Number of nearest neighbors to use for imputation. Default: 20.
        compute_uncertainty
            If ``True``, also computes the uncertainty of the imputation. Default: False.
        return_query_only
            If ``True``, return only the imputed query dataset as an ``AnnData`` object.
            Default: False.

        Returns
        -------
        Imputed ``AnnData`` object. If ``return_query_only`` is ``True``, only the query
        dataset is returned. If ``compute_uncertainty`` is ``True``, an uncertainty matrix
        is also returned.
        """
        adata = self.adata
        batch_key = self.batch_key

        # validate input
        validate_obsm_keys(adata, use_rep)
        validate_layer_key(adata_rna, layer_key)

        # retrieve reference and query indices
        reference_indices = adata.obs_names[adata.obs[batch_key] == reference_batch]
        query_indices = adata.obs_names[adata.obs[batch_key] != reference_batch]

        # validate that query indices are in to impute adata
        if not all(idx in adata_rna.obs_names for idx in reference_indices):
            raise ValueError("Some query indices are not present in `adata_to_impute`.")

        # get representations
        if use_rep is None:
            obsm_keys = adata.obsm.keys()

            if CYTOVI_DEFAULT_REP in obsm_keys:
                use_rep = CYTOVI_DEFAULT_REP
            else:
                adata.obsm[CYTOVI_DEFAULT_REP] = self.get_latent_representation()
                use_rep = CYTOVI_DEFAULT_REP

        # Get representations and reference expression
        rep_ref = adata[reference_indices, :].obsm[use_rep]
        rep_query = adata[query_indices, :].obsm[use_rep]
        expr_data_ref = adata_rna[reference_indices, :].layers[layer_key]

        # Impute expression in query
        imputed_expr_query, uncertainty = impute_expr_with_neighbors(
            rep_query,
            rep_ref,
            expr_data_ref,
            n_neighbors=n_neighbors,
            compute_uncertainty=compute_uncertainty,
        )

        # create anndata for imputed query dataset
        adata_imputed_query = AnnData(
            X=imputed_expr_query,
            obs=adata[query_indices, :].obs,
            obsm=adata[query_indices, :].obsm,
            var=adata_rna.var,
            layers={layer_key: imputed_expr_query},
        )

        if return_query_only:
            return adata_imputed_query

        # assemble new anndata with imputed expression
        expr_comb = np.concatenate([expr_data_ref, imputed_expr_query], axis=0)
        obs_comb = adata.obs.loc[np.concatenate([reference_indices, query_indices]), :]

        # restore original indices and add metadata
        adata_combined = AnnData(X=expr_comb, obs=obs_comb, var=adata_rna.var)

        adata_imputed = AnnData(
            X=adata_combined[adata.obs_names].X,
            obs=adata.obs,
            var=adata_rna.var,
            obsm=adata.obsm,
            layers={layer_key: adata_combined[adata.obs_names].X},
        )

        return adata_imputed
