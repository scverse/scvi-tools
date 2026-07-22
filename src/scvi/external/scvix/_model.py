from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data._utils import _get_adata_minify_type
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LabelsWithUnlabeledObsField,
    LayerField,
    NumericalJointObsField,
)
from scvi.dataloaders import DataSplitter
from scvi.model._utils import (
    get_max_epochs_heuristic,
)
from scvi.model.base import (
    ArchesMixin,
    BaseMinifiedModeModelClass,
    EmbeddingMixin,
    RNASeqMixin,
    VAEMixin,
)
from scvi.train import AdversarialTrainingPlan, TrainRunner
from scvi.utils._docstrings import devices_dsp, setup_anndata_dsp

from ._module import VAEX

if TYPE_CHECKING:
    from typing import Literal

    import pandas as pd
    from anndata import AnnData

logger = logging.getLogger(__name__)


class SCVIX(
    EmbeddingMixin,
    RNASeqMixin,
    VAEMixin,
    ArchesMixin,
    BaseMinifiedModeModelClass,
):
    """single-cell Variational Inference :cite:p:`Lopez18`.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.SCVI.setup_anndata`. If
        ``None``, then the underlying module will not be initialized until training, and a
        :class:`~lightning.pytorch.core.LightningDataModule` must be passed in during training
        (``EXPERIMENTAL``).
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    dropout_rate
        Dropout rate for neural networks.
    dispersion
        One of the following:

        * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - dispersion can differ between different batches
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
    gene_likelihood
        One of:

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
        * ``'poisson'`` - Poisson distribution
        * ``'normal'`` - ``EXPERIMENTAL`` Normal distribution
    prior
        One of:
        * ``'gaussian'`` - Gaussian prior
        * ``'mog'`` - Mixture of Gaussians prior
        * ``'vamp'`` - Variational Amortized Mixture of Posteriors prior
        * ``'mog_celltype'`` - Mixture of Gaussians prior with cell-type bias
    pseudoinputs_data_indices
        Indices of cells to use as pseudoinputs for the VAMP prior. If ``None``, a random sample of
        ``n_prior_components`` cells will be used.
    n_prior_components
        Number of components to use for the VAMP and MOG priors. This is the number of pseudoinputs
        used for the VAMP prior, and the number of components in the MOG prior.
        Defaults to 50.
    **kwargs
        Additional keyword arguments for :class:`~scvi.module.VAE`.

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.model.SCVI(adata)
    >>> vae.train()
    >>> adata.obsm["X_scVI"] = vae.get_latent_representation()
    >>> adata.obsm["X_normalized_scVI"] = vae.get_normalized_expression()

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/quick_start/api_overview`
    2. :doc:`/tutorials/notebooks/scrna/harmonization`
    3. :doc:`/tutorials/notebooks/scrna/scarches_scvi_tools`
    4. :doc:`/tutorials/notebooks/scrna/scvi_in_R`

    See Also
    --------
    :class:`~scvi.module.VAE`
    """

    _module_cls = VAEX
    _LATENT_QZM_KEY = "scvix_latent_qzm"
    _LATENT_QZV_KEY = "scvix_latent_qzv"
    _data_splitter_cls = DataSplitter
    _training_plan_cls = AdversarialTrainingPlan
    _train_runner_cls = TrainRunner

    def __init__(
        self,
        adata: AnnData | None = None,
        n_hidden: int = 512,
        n_latent: int = 20,
        n_layers: int = 2,
        dropout_rate: float = 0.05,
        dispersion: Literal["gene", "gene-batch", "gene-assay", "gene-cell"] = "gene-assay",
        gene_likelihood: Literal["zinb", "nb", "poisson", "normal"] = "nb",
        prior: Literal["gaussian", "mog", "vamp", "mog_celltype"] = "gaussian",
        batch_representation: Literal["one-hot", "embedding"] = "embedding",
        batch_embedding_kwargs: dict | None = None,
        pseudoinputs_data_indices: np.array | None = None,
        n_prior_components: int = 50,
        **kwargs,
    ):
        super().__init__(adata)

        if batch_embedding_kwargs is None:
            batch_embedding_kwargs = {"variational": True, "embedding_dim": n_latent}

        self._module_kwargs = {
            "n_hidden": n_hidden,
            "n_latent": n_latent,
            "n_layers": n_layers,
            "dropout_rate": dropout_rate,
            "dispersion": dispersion,
            "gene_likelihood": gene_likelihood,
            "batch_representation": batch_representation,
            "batch_embedding_kwargs": batch_embedding_kwargs,
            **kwargs,
        }
        self._model_summary_string = (
            "SCVI model with the following parameters: \n"
            f"n_hidden: {n_hidden}, n_latent: {n_latent}, n_layers: {n_layers}, "
            f"dropout_rate: {dropout_rate}, dispersion: {dispersion}, "
            f"gene_likelihood: {gene_likelihood}, "
            f"batch_representation: {batch_representation}."
        )

        if prior == "vamp":
            if pseudoinputs_data_indices is None:
                pseudoinputs_data_indices = np.random.randint(
                    0, self.summary_stats.n_cells, n_prior_components
                )
            assert pseudoinputs_data_indices.shape[0] == n_prior_components
            assert pseudoinputs_data_indices.ndim == 1
            pseudoinput_data = next(
                iter(
                    self._make_data_loader(
                        adata=adata,
                        indices=pseudoinputs_data_indices,
                        batch_size=n_prior_components,
                        shuffle=False,
                    )
                )
            )
        else:
            pseudoinput_data = None

        n_cats_per_cov = (
            self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )
        self.module = self._module_cls(
            n_input=self.summary_stats.n_vars,
            n_batch=self.summary_stats.n_batch,
            n_assay=self.summary_stats.n_assay,
            n_adversarial_group=self.summary_stats.get("n_adversarial_group", 1),
            n_labels=self.summary_stats.get("n_labels", 1),
            n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
            n_cats_per_cov=n_cats_per_cov,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            batch_representation=batch_representation,
            batch_embedding_kwargs=batch_embedding_kwargs,
            prior=prior,
            pseudoinput_data=pseudoinput_data,
            n_prior_components=n_prior_components,
            **kwargs,
        )
        self.module.minified_data_type = self.minified_data_type

        self.init_params_ = self._get_init_params(locals())

    @devices_dsp.dedent
    def train(
        self,
        max_epochs: int | None = None,
        lr: float = 4e-3,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float | None = None,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        batch_size: int = 1024,
        early_stopping: bool = False,
        n_epochs_kl_warmup: int | None = 5,
        adversarial_classifier: bool | None = None,
        adversarial_key: str = "assay",
        scale_adversarial_loss: float | str = 5.0,
        adversarial_steps: int = 5,
        weight_kl_sample: float = 1e-5,
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        external_indexing: list[np.array] = None,
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
            Whether to shuffle indices before splitting. If `False`, the val, train, and test set
            are split in the sequential order of the data according to `validation_size` and
            `train_size` percentages.
        batch_size
            Minibatch size to use during training.
        early_stopping
            Whether to perform early stopping with respect to the validation set.
        n_epochs_kl_warmup
            Number of epochs to scale weight on KL divergences from 0 to 1.
        adversarial_key
            Key in `adata.obs` that corresponds to the batch or assay key to use for adversarial
            training. Defaults to the assay key.
        scale_adversarial_loss
            Scaling factor for the adversarial loss component. Higher values enforce stronger
            assay mixing. Use ``"auto"`` to scale automatically with the KL warmup weight.
        adversarial_steps
            Number of adversarial classifier gradient steps per training step.
        weight_kl_sample
            Weight on the KL divergence of the variational batch embeddings.
        datasplitter_kwargs
            Additional keyword arguments passed into :class:`~scvi.dataloaders.DataSplitter`.
        plan_kwargs
            Keyword args for :class:`~scvi.train.AdversarialTrainingPlan`. Keyword arguments passed
            to `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        external_indexing
            A list of data split indices in the order of training, validation, and test sets.
            Validation and test set are not required and can be left empty.
        **kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        if adversarial_classifier is None:
            if self.module.n_assay > 1:
                adversarial_classifier = True
            else:
                adversarial_classifier = False
        n_epochs_kl_warmup = (
            n_epochs_kl_warmup if n_epochs_kl_warmup is not None else max_epochs // 2
        )

        update_dict = {
            "lr": lr,
            "adversarial_classifier": adversarial_classifier,
            "adversarial_key": adversarial_key,
            "n_epochs_kl_warmup": n_epochs_kl_warmup,
            "scale_adversarial_loss": scale_adversarial_loss,
            "adversarial_steps": adversarial_steps,
            "weight_kl_sample": weight_kl_sample,
        }
        if plan_kwargs is not None:
            plan_kwargs.update(update_dict)
        else:
            plan_kwargs = update_dict

        if max_epochs is None:
            max_epochs = get_max_epochs_heuristic(self.adata.n_obs)

        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else {}
        datasplitter_kwargs = datasplitter_kwargs or {}

        data_splitter = self._data_splitter_cls(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            batch_size=batch_size,
            external_indexing=external_indexing,
            **datasplitter_kwargs,
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
            **kwargs,
        )
        return runner()

    def _get_transform_batch_gen_kwargs(self, batch):
        kwargs = super()._get_transform_batch_gen_kwargs(batch)
        if hasattr(self, "_transform_assay_override"):
            kwargs["transform_assay"] = self._transform_assay_override
        return kwargs

    def get_normalized_expression(
        self,
        adata: AnnData | None = None,
        indices: list[int] | None = None,
        transform_batch: list[int | str] | None = None,
        transform_assay: int | str | None = None,
        gene_list: list[str] | None = None,
        library_size: float | Literal["latent"] = 1,
        n_samples: int = 1,
        n_samples_overall: int | None = None,
        weights: Literal["uniform", "importance"] | None = None,
        batch_size: int | None = None,
        return_mean: bool = True,
        return_numpy: bool | None = None,
        silent: bool = True,
        dataloader=None,
        data_loader_kwargs: dict | None = None,
        **importance_weighting_kwargs,
    ) -> np.ndarray | pd.DataFrame:
        """Returns the normalized (decoded) gene expression.

        Extends the base ``get_normalized_expression`` with assay conditioning via
        ``transform_assay``.

        Parameters
        ----------
        transform_assay
            Assay to condition on for all cells when decoding. If a string, must be a valid
            assay category registered via :meth:`~scvi.external.SCVIX.setup_anndata`. If an
            integer, used directly as the assay index. If ``None`` (default), each cell is
            decoded using its own observed assay.
        transform_batch
            Batch to condition on. See base class for details.

        Notes
        -----
        All other parameters are identical to
        :meth:`~scvi.model.base.RNASeqMixin.get_normalized_expression`.
        """
        if transform_assay is not None:
            adata_ = self._validate_anndata(adata)
            assay_mappings = (
                self.get_anndata_manager(adata_, required=True)
                .get_state_registry(REGISTRY_KEYS.ASSAY_KEY)
                .categorical_mapping
            )
            # Hack for transform_assay to be passed in get_normalized_expression.
            if isinstance(transform_assay, str):
                if transform_assay not in assay_mappings:
                    raise ValueError(f'"{transform_assay}" is not a valid assay category.')
                self._transform_assay_override = int(
                    np.where(assay_mappings == transform_assay)[0][0]
                )
            else:
                self._transform_assay_override = int(transform_assay)

        try:
            result = super().get_normalized_expression(
                adata=adata,
                indices=indices,
                transform_batch=transform_batch,
                gene_list=gene_list,
                library_size=library_size,
                n_samples=n_samples,
                n_samples_overall=n_samples_overall,
                weights=weights,
                batch_size=batch_size,
                return_mean=return_mean,
                return_numpy=return_numpy,
                silent=silent,
                dataloader=dataloader,
                data_loader_kwargs=data_loader_kwargs,
                **importance_weighting_kwargs,
            )
        finally:
            if hasattr(self, "_transform_assay_override"):
                del self._transform_assay_override

        return result

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        batch_key: str | None = None,
        assay_key: str | None = None,
        labels_key: str | None = None,
        adversarial_group_key: str | None = None,
        unlabeled_category: str = "unlabeled",
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        %(param_batch_key)s
        assay_key
            Key in ``adata.obs`` that corresponds to the assay of the data.
        %(param_labels_key)s
        adversarial_group_key
            Key in ``adata.obs`` that corresponds to the adversarial group for adversarial
            training. If ``None``, performs no conditional adversarial training.
        %(param_unlabeled_category)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.ASSAY_KEY, assay_key),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
        ]
        if labels_key is not None:
            anndata_fields.append(
                LabelsWithUnlabeledObsField(
                    REGISTRY_KEYS.LABELS_KEY, labels_key, unlabeled_category
                )
            )
        if adversarial_group_key is not None:
            anndata_fields.append(
                CategoricalObsField(REGISTRY_KEYS.ADVERSARIAL_GROUP_KEY, adversarial_group_key)
            )
        # register new fields if the adata is minified
        adata_minify_type = _get_adata_minify_type(adata)
        if adata_minify_type is not None:
            anndata_fields += cls._get_fields_for_adata_minification(adata_minify_type)
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
