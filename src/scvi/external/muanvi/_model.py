import logging
import warnings
from collections.abc import Sequence
from copy import deepcopy
from typing import Literal

import numpy as np
import pandas as pd
import torch
from anndata import AnnData

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data._constants import _SETUP_ARGS_KEY
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LabelsWithUnlabeledObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
)
from scvi.dataloaders import SemiSupervisedDataSplitter
from scvi.model import SCVI
from scvi.model._utils import get_max_epochs_heuristic, parse_device_args
from scvi.model.base import ArchesMixin, BaseModelClass, RNASeqMixin, VAEMixin
from scvi.model.base._archesmixin import _get_loaded_data, _set_params_online_update
from scvi.model.base._save_load import (
    _initialize_model,
    _validate_var_names,
)
from scvi.train import SemiSupervisedTrainingPlan, TrainRunner
from scvi.train._callbacks import SubSampleLabels
from scvi.utils import setup_anndata_dsp
from scvi.utils._docstrings import devices_dsp

from ._module import MUANVAE
from ._utils import LabelsWithUnlabeledJointObsField, _get_site_code_from_category

logger = logging.getLogger(__name__)


class MUANVI(RNASeqMixin, VAEMixin, ArchesMixin, BaseModelClass):
    """
    Hierarchical multi-annotator Variational Inference [Xu21]_.

    Inspired from M1 + M2 model, as described in (https://arxiv.org/pdf/1406.5298.pdf).

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.MUANVI.setup_anndata`.
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
        * ``'gene-label'`` - dispersion can differ between different labels
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
    gene_likelihood
        One of:
        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
        * ``'poisson'`` - Poisson distribution
    update_yprior
        Whether to perform the hierarchical update of the y prior parameter in the loss
    batches_to_harmonize
        List of two indices of the two batches to label-harmonize. Has to be defined if the dataset has more than 2 batches.
    **model_kwargs
        Keyword args for :class:`~scvi.module.MUANVAE`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.external.MUANVI.setup_anndata(adata, labels=["labels_0", "labels_1"], unknown_categories=["unknown, "unknown"])
    >>> model = scvi.external.MUANVI(adata)
    >>> model.train()
    >>> adata.obsm["X_muanvi"] = model.get_latent_representation()
    >>> adata.obs["pred_label_coarse"] = model.predict()[0]
    >>> adata.obs["pred_label_fine"] = model.predict()[1]

    """

    _module_cls = MUANVAE
    _training_plan_cls = SemiSupervisedTrainingPlan

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "nb",
        update_yprior: bool = True,
        eps_yprior: float = 1e-4,
        **model_kwargs,
    ):
        super().__init__(adata)
        muanvae_model_kwargs = dict(model_kwargs)
        self._set_indices_and_labels()

        n_batch = self.summary_stats.n_batch
        n_fine_labels = self.summary_stats.n_labels - 1
        n_site = self.summary_stats.n_site
        n_assay = self.summary_stats.n_assay

        self.hierarchy_dict, self.num_classes, self.hierarchy_matrix = self.extract_hierarchy(
            n_site=n_site, eps_yprior=eps_yprior
        )
        n_cats_per_cov = (
            self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )

        use_size_factor_key = REGISTRY_KEYS.SIZE_FACTOR_KEY in self.adata_manager.data_registry

        self.module = self._module_cls(
            n_input=self.summary_stats.n_vars,
            n_batch=n_batch,
            n_site=n_site,
            n_assay=n_assay,
            n_fine_labels = n_fine_labels,
            num_classes=self.num_classes,
            n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
            n_cats_per_cov=n_cats_per_cov,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            use_size_factor_key=use_size_factor_key,
            hierarchy_dict=self.hierarchy_dict,
            hierarchy_matrix=self.hierarchy_matrix,
            update_yprior=update_yprior,
            **muanvae_model_kwargs,
        )

        self.unsupervised_history_ = None
        self.semisupervised_history_ = None

        self._model_summary_string = (
            f"muANVI Model with the following params: \nunlabeled_category: {self.unlabeled_category}, n_hidden: {n_hidden}, n_latent: {n_latent}"
            f", n_layers: {n_layers}, dropout_rate: {dropout_rate}, dispersion: {dispersion}, gene_likelihood: {gene_likelihood}"
        )
        self.init_params_ = self._get_init_params(locals())
        self.was_pretrained = False
        self.n_fine_labels = n_fine_labels

    @classmethod
    def from_scvi_model(
        cls,
        scvi_model: SCVI,
        unlabeled_category: list[str],
        fine_labels: str | None = None,
        label_hierarchy: list[str] | None = None,
        adata: AnnData | None = None,
        **muanvi_kwargs,
    ):
        """
        Initialize scHANVI model with weights from pretrained :class:`~scvi.model.SCVI` model.

        Parameters
        ----------
        scvi_model
            Pretrained scvi model
        fine_labels
            key in `adata.obs` for label information. If this value is not None, the key will
            overwrite the `labels_key` used to setup AnnData with scvi.
        label_hierarchy
            List of strings, with levels of the cell-type hierarchy. Full hierarchy is inferred by
            concatenating label_hierarchy and fine_labels.
        unlabeled_category
            Value used for unlabeled cells in `labels_key`.
        adata
            AnnData object that has been registered via :meth:`~scvi.model.MUANVI.setup_anndata`.
        muanvi_kwargs
            kwargs for muANVI model
        """
        scvi_model._check_if_trained(message="Passed in scvi model hasn't been trained yet.")

        muanvi_kwargs = dict(muanvi_kwargs)
        init_params = scvi_model.init_params_
        non_kwargs = init_params["non_kwargs"]
        kwargs = init_params["kwargs"]
        kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        for k, v in {**non_kwargs, **kwargs}.items():
            if k in muanvi_kwargs.keys():
                warnings.warn(
                    f"Ignoring param '{k}' as it was already passed in to "
                    + f"pretrained scvi model with value {v}.",
                    stacklevel=2,
                )
                del muanvi_kwargs[k]

        if adata is None:
            adata = scvi_model.adata
        else:
            # validate new anndata against old model
            scvi_model._validate_anndata(adata)

        scvi_setup_args = deepcopy(scvi_model.adata_manager.registry[_SETUP_ARGS_KEY])
        scvi_labels_key = scvi_setup_args["labels_key"]
        if fine_labels is None and scvi_labels_key is None:
            raise ValueError(
                "A `labels_key` list is necessary as the SCVI model was initialized without one."
            )
        if fine_labels is not None:
            scvi_setup_args.update({"fine_labels": fine_labels})
            scvi_setup_args.pop("labels_key", None)
        else:
            scvi_setup_args["fine_labels"] = scvi_setup_args.pop("labels_key")

        cls.setup_anndata(
            adata,
            unlabeled_category=unlabeled_category,
            label_hierarchy=label_hierarchy,
            **scvi_setup_args,
        )
        muanvi_model = cls(adata, **non_kwargs, **kwargs, **muanvi_kwargs)
        scvi_state_dict = scvi_model.module.state_dict()
        muanvi_model.module.load_state_dict(scvi_state_dict, strict=False)
        muanvi_model.was_pretrained = True

        return muanvi_model

    def _set_indices_and_labels(self):
        """Set indices for labeled and unlabeled cells."""
        labels_state_registry = self.adata_manager.get_state_registry("label_hierarchy")
        self.original_label_keys = labels_state_registry.field_keys
        self.unlabeled_category = labels_state_registry.unlabeled_category

        # Dataframe of 2 columns for the 2 layers of labels
        labels = {field: self.adata.obs[field] for field in self.original_label_keys}
        self.labels = pd.DataFrame(labels)
        self._label_mapping = labels_state_registry.mappings.to_dict()
        # a cell is unlabeled if it is not labeled at the finest state
        labeled_indices_list = [
            set(np.where(self.labels.iloc[:, idx] != self.unlabeled_category)[0])
            for idx, _ in enumerate(self.original_label_keys)
        ]
        self._labeled_indices = list(set.intersection(*labeled_indices_list))
        self._unlabeled_indices = list(
            set(np.arange(self.adata.n_obs)) - set(self._labeled_indices)
        )
        self._code_to_label = {
            layer: dict(enumerate(self._label_mapping[layer])) for layer in self._label_mapping
        }

    def extract_hierarchy(self, n_site, eps_yprior):
        """
        Method to extract automatically the intrinsic hierarchy in the data.

        Parameters
        ----------
        n_site
            Number of different vocabulary sites in the data
        eps_yprior
            Epsilon value to add to the y prior parameter to avoid overconfidence
        """
        labels_state_registry = self.adata_manager.get_state_registry("label_hierarchy")
        label_keys = labels_state_registry.field_keys

        def fixed_depth_groupby(df, label_keys):
            if len(label_keys)==1:
                return list(df[label_keys[-1]].unique()) if not df.empty else []
            else:
                # Recursive case: build dictionaries up to the fixed depth
                return {
                    key: fixed_depth_groupby(sub_df, label_keys[1:])
                    for key, sub_df in df.groupby(label_keys[0], observed=True)
                }

        num_classes = labels_state_registry.n_cats_per_key
        hierarchy_dict = fixed_depth_groupby(self.labels, label_keys)
        hierarchy_matrix = []

        for n_label in range(1, len(num_classes)):
            curr = pd.DataFrame(
                0,
                index=np.arange(num_classes[n_label - 1]),
                columns=np.arange(num_classes[n_label]),
            )
            # Site specific last layer.
            if n_label == len(num_classes) - 1:
                hierarchy_matrix_ = torch.zeros(
                    num_classes[n_label - 1], num_classes[n_label], n_site
                )
                for site in range(n_site):
                    adata = self.adata[self.adata.obs["_scvi_site"] == site]

                    if adata.n_obs > 0:
                        curr_ = pd.crosstab(
                            adata.obsm["_scvi_label_hierarchy"][
                                label_keys[n_label - 1]
                            ],
                            adata.obsm["_scvi_label_hierarchy"][label_keys[n_label]],
                        )
                        curr_[curr_ > 0] = 1
                        curr_ = curr_.loc[curr.index, curr.columns]
                        curr = curr_.div(curr_.sum(axis=1), axis=0).fillna(0)
                        hierarchy_matrix_[:, :, site] = torch.tensor(curr.values)
            else:
                curr_ = pd.crosstab(
                    self.adata.obsm["_scvi_label_hierarchy"][label_keys[n_label - 1]],
                    self.adata.obsm["_scvi_label_hierarchy"][label_keys[n_label]],
                )
                curr_[curr_ > 0] = 1
                curr_ = curr_.loc[curr.index, curr.columns]
                curr = curr_.div(curr_.sum(axis=1), axis=0).fillna(0)
                hierarchy_matrix_ = torch.tensor(curr.values)

            hierarchy_matrix.append(hierarchy_matrix_)

        return (hierarchy_dict, num_classes, hierarchy_matrix)

    def predict(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        soft: bool = False,
        batch_size: int | None = None,
        level: int = -1,
        sites_to_predict: int | str | str | None = None,
    ) -> np.ndarray | pd.DataFrame:
        """
        Return cell label predictions.

        Parameters
        ----------
        adata
            AnnData object that has been registered via :meth:`~scvi.model.SCANVI.setup_anndata`.
        indices
            indices for which to return probabilities.
        soft
            If True, returns per class probabilities
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        site_to_predict
            For cross prediction purposes : sites to use when cross-predicting labels.
            If None, normal prediction occurs and each cell is labeled accordingly to
            its own site-specific classifier. Can be list of strings to use multiple sites,
            or a string to use a single site.
        level
            Level of the hierarchy to predict. If -1, predicts at the finest level.
        """
        adata = self._validate_anndata(adata)

        if indices is None:
            indices = np.arange(adata.n_obs)

        scdl = self._make_data_loader(
            adata=adata,
            indices=indices,
            batch_size=batch_size,
        )
        # total depth of the hierarchy
        total_level = len(self.module.num_classes)
        class_labels = self.adata_manager.get_state_registry("label_hierarchy").field_keys

        sites_to_predict_, site_mappings_ = _get_site_code_from_category(
            self.get_anndata_manager(adata, required=True), sites_to_predict
        )
        pred = {str(site_mappings_[i]) + "_fine": [] for i in sites_to_predict_ if i is not None}
        if None in sites_to_predict_:
            for i in class_labels:
                pred[i] = []

        for _, tensors in enumerate(scdl):
            for site_to_predict_ in sites_to_predict_:
                probs_, _ = self.module.classification(tensors, site_to_predict=site_to_predict_)
                if site_to_predict_ is not None:
                    if not soft:
                        pred_ = probs_.argmax(dim=1)
                    else:
                        pred_ = probs_
                    pred[site_mappings_[site_to_predict_] + "_fine"].append(pred_.detach().cpu())

                else:
                    # pred is a tuple of probabilities and logits for each layer
                    probs = probs_

                    for i in range(total_level):
                        if not soft:
                            pred[class_labels[i]].append(probs[i].argmax(dim=1).cpu())
                        else:
                            pred[class_labels[i]].append(probs[i].detach().cpu())

        for key in pred.keys():
            pred[key] = torch.cat(pred[key]).numpy()
            if not soft:
                pred[key] = [self._code_to_label[key][ct] for ct in pred[key]]

        if not soft:
            pred = pd.DataFrame.from_dict(pred)
            pred.index = adata.obs_names[indices]
            return pred
        else:
            for key in pred.keys():
                columns = list(self._code_to_label[key].values())[:-1]

                pred[key] = pd.DataFrame(
                    pred[key],
                    columns=columns,
                    index=adata.obs_names[indices],
                )
            return pred

    @classmethod
    @devices_dsp.dedent
    def load_query_data(
        cls,
        adata: AnnData,
        reference_model: str | BaseModelClass,
        inplace_subset_query_vars: bool = False,
        accelerator: str = "auto",
        device: int | str = "auto",
        unfrozen: bool = False,
        freeze_dropout: bool = False,
        freeze_expression: bool = True,
        freeze_decoder_first_layer: bool = True,
        freeze_batchnorm_encoder: bool = True,
        freeze_batchnorm_decoder: bool = False,
        freeze_classifier: bool = True,
    ):
        """Online update of a reference model with scArches algorithm :cite:p:`Lotfollahi21`.

        Parameters
        ----------
        adata
            AnnData organized in the same way as data used to train model.
            It is not necessary to run setup_anndata,
            as AnnData is validated against the ``registry``.
        reference_model
            Either an already instantiated model of the same class, or a path to
            saved outputs for reference model.
        inplace_subset_query_vars
            Whether to subset and rearrange query vars inplace based on vars used to
            train reference model.
        %(param_accelerator)s
        %(param_device)s
        unfrozen
            Override all other freeze options for a fully unfrozen model
        freeze_dropout
            Whether to freeze dropout during training
        freeze_expression
            Freeze neurons corersponding to expression in first layer
        freeze_decoder_first_layer
            Freeze neurons corresponding to first layer in decoder
        freeze_batchnorm_encoder
            Whether to freeze batchnorm weight and bias during training for encoder
        freeze_batchnorm_decoder
            Whether to freeze batchnorm weight and bias during training for decoder
        freeze_classifier
            Whether to freeze classifier completely. Only applies to `SCANVI`.
        """
        _, _, device = parse_device_args(
            accelerator=accelerator,
            devices=device,
            return_device="torch",
            validate_single_device=True,
        )

        attr_dict, var_names, load_state_dict = _get_loaded_data(reference_model, device=device)

        if inplace_subset_query_vars:
            logger.debug("Subsetting query vars to reference vars.")
            adata._inplace_subset_var(var_names)
        _validate_var_names(adata, var_names)

        registry = attr_dict.pop("registry_")
        if _SETUP_ARGS_KEY not in registry:
            raise ValueError(
                "Saved model does not contain original setup inputs. "
                "Cannot load the original setup."
            )

        cls.setup_anndata(
            adata,
            source_registry=registry,
            extend_categories=True,
            allow_missing_labels=True,
            **registry[_SETUP_ARGS_KEY],
        )

        model = _initialize_model(cls, adata, attr_dict)
        adata_manager = model.get_anndata_manager(adata, required=True)

        if REGISTRY_KEYS.CAT_COVS_KEY in adata_manager.data_registry:
            raise NotImplementedError(
                "scArches currently does not support models with extra categorical covariates."
            )

        model.to_device(device)

        # model tweaking
        new_state_dict = model.module.state_dict()
        additional_parameters = set()
        for key, new_ten in new_state_dict.items():
            load_ten = load_state_dict.get(key, None)
            if load_ten is None:
                # Picks up that additional site classifier was added, makes it trainable by default.
                if "y_prior_fine" not in key:
                    additional_parameters.add(key)  # TODO check this.
                load_state_dict[key] = new_ten
                continue
            if new_ten.size() == load_ten.size():
                continue
            # new categoricals changed size
            else:
                if new_ten.size()[0] != load_ten.size()[0]:
                    new_ten = new_ten.to(load_ten.device)
                    dim_diff = new_ten.size()[0] - load_ten.size()[0]
                    fixed_ten = torch.cat([load_ten, new_ten[-dim_diff:, ...]], dim=0)
                    load_state_dict[key] = fixed_ten
                else:
                    new_ten = new_ten.to(load_ten.device)
                    dim_diff = new_ten.size()[-1] - load_ten.size()[-1]
                    fixed_ten = torch.cat([load_ten, new_ten[..., -dim_diff:]], dim=-1)
                    load_state_dict[key] = fixed_ten

        model.module.load_state_dict(load_state_dict)
        model.module.eval()

        _set_params_online_update(
            model.module,
            unfrozen=unfrozen,
            freeze_decoder_first_layer=freeze_decoder_first_layer,
            freeze_batchnorm_encoder=freeze_batchnorm_encoder,
            freeze_batchnorm_decoder=freeze_batchnorm_decoder,
            freeze_dropout=freeze_dropout,
            freeze_expression=freeze_expression,
            freeze_classifier=freeze_classifier,
            parameters_yes_grad=additional_parameters,
        )
        model.is_trained_ = False

        return model

    def train(
        self,
        max_epochs: int | None = None,
        n_samples_per_label: float | None = None,
        check_val_every_n_epoch: int | None = None,
        train_size: float = 0.9,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        batch_size: int = 128,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        **trainer_kwargs,
    ):
        """Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset for semisupervised training.
        n_samples_per_label
            Number of subsamples for each label class to sample per epoch. By default, there
            is no label subsampling.
        check_val_every_n_epoch
            Frequency with which metrics are computed on the data for validation set for both
            the unsupervised and semisupervised trainers. If you'd like a different frequency for
            the semisupervised trainer, set check_val_every_n_epoch in semisupervised_train_kwargs.
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
        %(param_accelerator)s
        %(param_devices)s
        datasplitter_kwargs
            Additional keyword arguments passed into
            :class:`~scvi.dataloaders.SemiSupervisedDataSplitter`.
        plan_kwargs
            Keyword args for :class:`~scvi.train.SemiSupervisedTrainingPlan`. Keyword arguments
            passed to `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        if max_epochs is None:
            max_epochs = get_max_epochs_heuristic(self.adata.n_obs)

            if self.was_pretrained:
                max_epochs = int(np.min([10, np.max([2, round(max_epochs / 3.0)])]))

        plan_kwargs = {} if plan_kwargs is None else plan_kwargs
        datasplitter_kwargs = datasplitter_kwargs or {}

        # if we have labeled cells, we want to subsample labels each epoch
        sampler_callback = [SubSampleLabels()] if len(self._labeled_indices) != 0 else []

        data_splitter = SemiSupervisedDataSplitter(
            adata_manager=self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            n_samples_per_label=n_samples_per_label,
            batch_size=batch_size,
            **datasplitter_kwargs,
        )

        warmup_epochs = plan_kwargs.pop("warmup_epochs", None)

        if warmup_epochs is not None and warmup_epochs > 0:
            logger.info(f"Pretraining for {max_epochs} epochs.")

            plan_kwargs_pre = plan_kwargs.copy()
            plan_kwargs_pre["warmup_model"] = True
            plan_kwargs_pre["n_epochs_kl_warmup"] = warmup_epochs

            training_plan = self._training_plan_cls(
                self.module, n_classes=self.n_fine_labels, **plan_kwargs_pre
            )  # n_classes set at the finest level to track accuracy at that level.
            runner_pre = TrainRunner(
                self,
                training_plan=training_plan,
                data_splitter=data_splitter,
                max_epochs=warmup_epochs,
                accelerator=accelerator,
                devices=devices,
                check_val_every_n_epoch=check_val_every_n_epoch,
                **trainer_kwargs,
            )
            runner_pre()
            self.was_pretrained = True

        logger.info(f"Training for {max_epochs} epochs.")

        if "callbacks" in trainer_kwargs.keys():
            trainer_kwargs["callbacks"] + [sampler_callback]
        else:
            trainer_kwargs["callbacks"] = sampler_callback
        training_plan = self._training_plan_cls(
            self.module, n_classes=self.n_fine_labels, **plan_kwargs
        )

        runner = TrainRunner(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=devices,
            check_val_every_n_epoch=check_val_every_n_epoch,
            **trainer_kwargs,
        )

        return runner()

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        fine_labels_key: str,
        unlabeled_category: list[str | int | float],
        layer: str | None = None,
        site_key: str | None = None,
        assay_key: str | None = None,
        batch_key: str | None = None,
        size_factor_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        label_hierarchy: list[str] | None = None,
        **kwargs,
    ):
        """
        %(summary)s.

        Parameters
        ----------
        %(param_layer)s
        %(param_batch_key)s
        %(param_site_key)s
        %(param_assay_key)s
        fine_labels_key
            key in `adata.obs` for fine label information. Categories will automatically be
            converted into integer categories and saved to `adata.obs['_scvi_labels']`.
            If `None`, assigns the same label to all the data. This information can be
            site-specific. In this case we expect all the labels in a single obs column.
            This is analogous to the label key in scANVI.
        %(param_size_factor_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        label_hierarchy
            List of strings, with levels of the cell-type hierarchy.
            The first list is the root level, the last list is the second finest level.
            The full hierarchy is inferred by concatenating label_hierarchy and fine_labels.
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.ASSAY_KEY, assay_key),
            CategoricalObsField(REGISTRY_KEYS.SITE_KEY, site_key),
            LabelsWithUnlabeledObsField(REGISTRY_KEYS.LABELS_KEY, fine_labels_key, unlabeled_category),
            NumericalObsField(REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
            LabelsWithUnlabeledJointObsField(
                "label_hierarchy", label_hierarchy+[fine_labels_key], unlabeled_category),
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
