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
from scvi.utils import setup_anndata_dsp
from scvi.utils._docstrings import devices_dsp, setup_anndata_dsp

from ._module import ASSAYVAE

if TYPE_CHECKING:
    from typing import Literal

    import numpy as np
    from anndata import AnnData

logger = logging.getLogger(__name__)
print(2)


class ASSAYVI(
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
    latent_distribution
        One of:

        * ``'normal'`` - Normal distribution
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)
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

    _module_cls = ASSAYVAE
    _LATENT_QZM_KEY = "assayvi_latent_qzm"
    _LATENT_QZV_KEY = "assayvi_latent_qzv"
    _data_splitter_cls = DataSplitter
    _training_plan_cls = AdversarialTrainingPlan
    _train_runner_cls = TrainRunner

    def __init__(
        self,
        adata: AnnData | None = None,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.05,
        dispersion: Literal["gene", "gene-batch", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson", "normal"] = "nb",
        prior: Literal["normal", "mog", "vamp"] = "normal",
        pseudoinputs_data_indices: np.array | None = None,
        n_prior_components: int = 50,
        **kwargs,
    ):
        super().__init__(adata)

        self._module_kwargs = {
            "n_hidden": n_hidden,
            "n_latent": n_latent,
            "n_layers": n_layers,
            "dropout_rate": dropout_rate,
            "dispersion": dispersion,
            "gene_likelihood": gene_likelihood,
            **kwargs,
        }
        self._model_summary_string = (
            "SCVI model with the following parameters: \n"
            f"n_hidden: {n_hidden}, n_latent: {n_latent}, n_layers: {n_layers}, "
            f"dropout_rate: {dropout_rate}, dispersion: {dispersion}, "
            f"gene_likelihood: {gene_likelihood}."
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
            n_labels=self.summary_stats.get("n_labels", 1),
            n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
            n_cats_per_cov=n_cats_per_cov,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
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
        batch_size: int = 256,
        early_stopping: bool = True,
        check_val_every_n_epoch: int | None = None,
        reduce_lr_on_plateau: bool = True,
        n_steps_kl_warmup: int | None = None,
        n_epochs_kl_warmup: int | None = None,
        adversarial_classifier: bool | None = None,
        adversarial_key: str = "assay",
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
        check_val_every_n_epoch
            Check val every n train epochs. By default, val is not checked, unless `early_stopping`
            is `True` or `reduce_lr_on_plateau` is `True`. If either of the latter conditions are
            met, val is checked every epoch.
        reduce_lr_on_plateau
            Reduce learning rate on plateau of validation metric (default is ELBO).
        n_epochs_kl_warmup
            Number of epochs to scale weight on KL divergences from 0 to 1.
        scale_adversarial_classifier
            How to weight adversarial classifier in the latent space. This helps mixing when
            there are multiple assays. Defaults to `1`.
        adversarial_key
            Key in `adata.obs` that corresponds to batch or assay key to use for adversarial
            training. If `None`, defaults to the assay key.
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
            n_epochs_kl_warmup if n_epochs_kl_warmup is not None else max_epochs//2
        )
        if reduce_lr_on_plateau:
            check_val_every_n_epoch = 1

        update_dict = {
            "lr": lr,
            "adversarial_classifier": adversarial_classifier,
            "adversarial_key": adversarial_key,
            "reduce_lr_on_plateau": reduce_lr_on_plateau,
            "n_epochs_kl_warmup": n_epochs_kl_warmup,
            "n_steps_kl_warmup": n_steps_kl_warmup,
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
            check_val_every_n_epoch=check_val_every_n_epoch,
            **kwargs,
        )
        return runner()

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        batch_key: str | None = None,
        assay_key: str | None = None,
        labels_key: str | None = None,
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
                    REGISTRY_KEYS.LABELS_KEY, labels_key, unlabeled_category))
        # register new fields if the adata is minified
        adata_minify_type = _get_adata_minify_type(adata)
        if adata_minify_type is not None:
            anndata_fields += cls._get_fields_for_adata_minification(adata_minify_type)
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
