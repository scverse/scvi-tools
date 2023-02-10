import io
import logging
import warnings
from contextlib import redirect_stdout
from typing import List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import torch
import anndata
from anndata import AnnData
from sklearn.preprocessing import OneHotEncoder

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import CategoricalObsField, LayerField, ObsmField
from scvi.dataloaders import DataSplitter
from scvi.model import SCVI
from scvi.model.base import BaseModelClass
from scvi.module import Classifier
from scvi.module.base import auto_move_data
from scvi.train import (
    ClassifierTrainingPlan,
    MultiBinaryClassifierTrainingPlan,
    LoudEarlyStopping,
    TrainRunner,
)
from scvi.utils import setup_anndata_dsp

logger = logging.getLogger(__name__)

LABELS_KEY = "_cpxvi_sim"
ONE_HOT_CLASS_KEY = "_cpxvi_one_hot"
SEP_STR = ":_:"


class CPXVI(BaseModelClass):
    """
    Multi barcode detection in single cell genomics data.

    Most users will initialize the model using the class method
    :meth:`~scvi.external.CPXVI.from_scvi_model`, which takes as
    input a pre-trained :class:`~scvi.model.SCVI` object.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.SCVI.setup_anndata`.
        Object should contain latent representation of real cells and multiplets as `adata.X`.
        Object should also be registered, using `.X` and `labels_key="_cpxvi_sim"`.
    **classifier_kwargs
        Keyword args for :class:`~scvi.module.Classifier`

    Examples
    --------
    In the case of scVI trained with multiple batches:

    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.model.SCVI(adata)
    >>> vae.train()
    >>> cpxvi_batch_1 = scvi.external.CPXVI.from_scvi_model(vae, restrict_to_batch="batch 1")
    >>> cpxvi_batch_1.train()
    >>> cpxvi_batch_1.predict()

    Otherwise:

    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.SCVI.setup_anndata(adata)
    >>> vae = scvi.model.SCVI(adata)
    >>> vae.train()
    >>> solo = scvi.external.CPXVI.from_scvi_model(vae)
    >>> CPXVI.train()
    >>> CPXVI.predict()

    Notes
    -----
    Solo should be trained on one lane of data at a time. An
    :class:`~scvi.model.SCVI` instance that was trained with multiple
    batches can be used as input, but Solo should be created and run
    multiple times, each with a new `restrict_to_batch` in
    :meth:`~scvi.external.CPXVI.from_scvi_model`.
    """

    def __init__(
        self,
        adata: AnnData,
        **classifier_kwargs,
    ):
        # TODO, catch user warning here and logger warning
        # about non count data
        super().__init__(adata)

        if ONE_HOT_CLASS_KEY not in adata.obsm_keys():
            msg = f"{ONE_HOT_CLASS_KEY} not in `adata.obsm`.\n"
            msg += "Try preparing with `CPXVI.from_scvi_model(...)`."
            raise ValueError(msg)

        # infer the number of labels from the one-hot classes
        n_labels = adata.obsm[ONE_HOT_CLASS_KEY].shape[1]

        # `logits=True` means the module outputs raw logits
        # rather than running a final activation internally
        self.module = Classifier(
            n_input=self.summary_stats.n_vars,
            n_labels=n_labels,
            logits=True,
            **classifier_kwargs,
        )
        self._model_summary_string = "cpxVI model"
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    def from_scvi_model(
        cls,
        scvi_model: SCVI,
        adata: Optional[AnnData] = None,
        restrict_to_batch: Optional[str] = None,
        multiplet_ratio: np.ndarray = np.array([2] * 2),
        **classifier_kwargs,
    ):
        """
        Instantiate a CPXVI model from an scvi model.

        Parameters
        ----------
        scvi_model
            Pre-trained model of :class:`~scvi.model.SCVI`. The
            adata object used to initialize this model should have only
            been setup with count data, and optionally a `batch_key`;
            i.e., no extra covariates or labels, etc.
        adata
            Optional anndata to use that is compatible with scvi_model.
        restrict_to_batch
            Batch category in `batch_key` used to setup adata for scvi_model
            to restrict cpxVI model to. This allows to train a CPXVI model on
            one batch of a scvi_model that was trained on multiple batches.
            AnnData object setup with setup_anndata. [cells, features].
        multiplet_ratio : np.ndarray
            [n_complexity,] ratio of generated multiplets to produce relative
            to number of cells in adata or length of indices, if not `None`.
            index 0 corresponds to doublets to generate, index 1 to triplets,
            and so on. maximum complexity of multiplets is inferred as
            `len(multiplet_ratio) + 1`.
        **classifier_kwargs
            Keyword args for :class:`~scvi.module.Classifier`

        Returns
        -------
        CPXVI model
        """
        _validate_scvi_model(scvi_model, restrict_to_batch=restrict_to_batch)
        orig_adata_manager = scvi_model.adata_manager
        orig_batch_key = orig_adata_manager.get_state_registry(
            REGISTRY_KEYS.BATCH_KEY
        ).original_key
        orig_labels_key = orig_adata_manager.get_state_registry(
            REGISTRY_KEYS.LABELS_KEY
        ).original_key

        if adata is not None:
            adata_manager = orig_adata_manager.transfer_fields(adata)
            cls.register_manager(adata_manager)
        else:
            adata_manager = orig_adata_manager
        adata = adata_manager.adata

        if restrict_to_batch is not None:
            batch_mask = adata.obs[orig_batch_key] == restrict_to_batch
            if np.sum(batch_mask) == 0:
                raise ValueError(
                    "Batch category given to restrict_to_batch not found.\n"
                    + "Available categories: {}".format(
                        adata.obs[orig_batch_key].astype("category").cat.categories
                    )
                )
            # indices in adata with restrict_to_batch category
            batch_indices = np.where(batch_mask)[0]
        else:
            # use all indices
            batch_indices = None

        # anndata with only generated multiplets
        multiplet_adata = cls.create_multiplets(
            adata_manager,
            indices=batch_indices,
            multiplet_ratio=multiplet_ratio,
            class_labels_key=orig_labels_key,
        )
        # if scvi wasn't trained with batch correction having the
        # zeros here does nothing.
        multiplet_adata.obs[orig_batch_key] = (
            restrict_to_batch if restrict_to_batch is not None else 0
        )

        # TODO: Remove. This breaks .get_latent_representation below because it adds
        # new categories to the label registry and `.validate_anndata(..., extend_categories=False)`
        # by default in the SCVI model operations
        # revert to dummy mode!

        # Create cat names in orig labels column in adata (does not affect inference).
        # cat_labels = []
        # for i in range(multiplet_adata.shape[0]):
        #     cat_labels.append(
        #         SEP_STR.join(
        #             multiplet_adata.uns[ONE_HOT_CLASS_KEY][
        #                 multiplet_adata.obsm[ONE_HOT_CLASS_KEY][i, :].astype(bool)
        #             ]
        #         )
        #     )
        # multiplet_adata.obs[orig_labels_key] = cat_labels
        # Create dummy labels column set to first label in adata (does not affect inference).
        dummy_label = orig_adata_manager.get_state_registry(
            REGISTRY_KEYS.LABELS_KEY
        ).categorical_mapping[0]
        multiplet_adata.obs[orig_labels_key] = dummy_label

        # save one hot's to add back later
        pre_merge_one_hots = np.concatenate(
            [adata.obsm[ONE_HOT_CLASS_KEY], multiplet_adata.obsm[ONE_HOT_CLASS_KEY]],
            axis=0,
        ).copy()

        # if model is using observed lib size, needs to get lib sample
        # which is just observed lib size on log scale
        give_mean_lib = not scvi_model.module.use_observed_lib_size

        # get latent representations and make input anndata
        latent_rep = scvi_model.get_latent_representation(adata, indices=batch_indices)
        lib_size = scvi_model.get_latent_library_size(
            adata, indices=batch_indices, give_mean=give_mean_lib
        )
        latent_adata = AnnData(np.concatenate([latent_rep, np.log(lib_size)], axis=1))
        latent_adata.obs[LABELS_KEY] = "singlet"
        orig_obs_names = adata.obs_names
        latent_adata.obs_names = (
            orig_obs_names[batch_indices]
            if batch_indices is not None
            else orig_obs_names
        )

        logger.info("Creating synthetic multiples, preparing CPXVI model.")
        f = io.StringIO()
        with redirect_stdout(f):
            multiplet_latent_rep = scvi_model.get_latent_representation(multiplet_adata)
            multiplet_lib_size = scvi_model.get_latent_library_size(
                multiplet_adata, give_mean=give_mean_lib
            )
            multiplet_adata = AnnData(
                np.concatenate(
                    [multiplet_latent_rep, np.log(multiplet_lib_size)], axis=1
                )
            )
            # TODO figure out how we can use multiple one hot labels rather than a categorical!
            multiplet_adata.obs[LABELS_KEY] = "multiplet"

            full_adata = latent_adata.concatenate(multiplet_adata)
            # add back one hots
            full_adata.obsm[ONE_HOT_CLASS_KEY] = pre_merge_one_hots
            full_adata.uns[ONE_HOT_CLASS_KEY] = adata.uns[ONE_HOT_CLASS_KEY]
            cls.setup_anndata(full_adata, labels_key=ONE_HOT_CLASS_KEY)
        return cls(full_adata, **classifier_kwargs)

    @staticmethod
    def labels2hotvec(class_labels: np.ndarray) -> Tuple[np.ndarray]:
        """Convert a vector of categorical labels to a one-hot binary matrix
        of class labels"""
        ohe = OneHotEncoder()
        class_one_hot = ohe.fit_transform(class_labels.reshape(-1, 1))
        class_one_hot = class_one_hot.toarray()
        return class_one_hot, ohe.categories_[0]

    @classmethod
    def make_one_hot_labels(cls, adata: AnnData, class_labels_key: str) -> None:
        """Adds one-hot labels for classes in `.obs[class_labels_key]` to a
        matrix in `.obsm[ONE_HOT_CLASS_KEY]` in-place"""
        true_class_labels = np.array(adata.obs[class_labels_key])
        one_hot_labels, one_hot_class_names = cls.labels2hotvec(true_class_labels)
        adata.obsm[ONE_HOT_CLASS_KEY] = one_hot_labels
        adata.uns[ONE_HOT_CLASS_KEY] = one_hot_class_names
        return

    @classmethod
    def create_multiplets(
        cls,
        adata_manager: AnnDataManager,
        multiplet_ratio: np.ndarray,
        class_labels_key: str,
        indices: Optional[Sequence[int]] = None,
        seed: int = 1,
        norm_lib_size: bool = True,
    ) -> AnnData:
        """Simulate multiplets.

        Parameters
        ----------
        adata
            AnnData object setup with setup_anndata. [cells, features].
        multiplet_ratio : np.ndarray
            [n_complexity,] ratio of generated multiplets to produce relative
            to number of cells in adata or length of indices, if not `None`.
            index 0 corresponds to doublets to generate, index 1 to triplets,
            and so on. maximum complexity of multiplets is inferred as
            `len(multiplet_ratio) + 1`.
        class_labels_key : str
            key for the class labels to be used.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        seed
            Seed for reproducibility.
        norm_lib_size : bool
            normalize the library size based on the number of cells used to create
            the synthetic multiplet. note, `solo` doesn't do this because the larger
            lib size is informative. here, it may not be.

        Returns
        -------
        multiplet_adata : AnnData
            [sum(n_obs*multiple_ratio), features] where `n_obs` is `adata.shape[0]`
            if `indices` is `None`, otherwise `n_obs` is `len(indices)`.
            `.obs_names` are `sim_multiplet_{count}_n_bcs_{n_cells_per_multiplet}`.
            `.var_names` match `adata` in `adata_manager`.

        Notes
        -----
        Currently include multiplets that have only a single label (synthetisized from
        >1 cells with the same true class). Could remove in the future.
        """
        adata = adata_manager.adata
        n_obs = adata.n_obs if indices is None else len(indices)
        max_complexity = len(multiplet_ratio) + 1

        # counts can be in many locations, this uses where it was registered in setup
        x = adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)
        if indices is not None:
            x = x[indices]

        # get the class labels provided
        cls.make_one_hot_labels(adata, class_labels_key)
        one_hot_labels = adata.obsm[ONE_HOT_CLASS_KEY].copy()
        one_hot_class_names = adata.uns[ONE_HOT_CLASS_KEY]

        multiplet_adatas = []
        for i_complexity in range(len(multiplet_ratio)):
            # each index value is two less than the number of cells for that ratio
            n_cells_per_multiplet = i_complexity + 2
            logger.info(f"Generating multiples with {n_cells_per_multiplet} cells.")
            num_multiplets = multiplet_ratio[i_complexity] * n_obs

            random_state = np.random.RandomState(seed=seed)
            # generate an [n_obs, n_cells_per_multiplet] array of randomized indices
            # each row is a random combination of cells to generate
            parent_inds = random_state.choice(
                n_obs, size=(num_multiplets, n_cells_per_multiplet)
            )
            multiplets = x[parent_inds[:, 0]]
            multiplets_classes = one_hot_labels[parent_inds[:, 0], :]
            for j_inds in range(1, parent_inds.shape[1]):
                multiplets += x[parent_inds[:, j_inds]]
                multiplets_classes += one_hot_labels[parent_inds[:, j_inds], :]
            # convert back to binary, now a "multihot" label
            multiplets_classes = multiplets_classes > 0

            if norm_lib_size:
                # normalize the total counts by the number of cells
                # used to synthesize the multiplet
                multiplets /= n_cells_per_multiplet
                # round back to integer counts
                multiplets = np.round(multiplets)

            multiplets_ad = AnnData(multiplets)
            multiplets_ad.var_names = adata.var_names
            multiplets_ad.obs_names = [
                f"sim_multiplet_{i}_n_bcs_{n_cells_per_multiplet}"
                for i in range(num_multiplets)
            ]

            # if adata setup with a layer, need to add layer to multiplets adata
            layer = adata_manager.data_registry[REGISTRY_KEYS.X_KEY].attr_key
            if layer is not None:
                multiplets_ad.layers[layer] = multiplets

            # add one hot classes
            multiplets_ad.obsm[ONE_HOT_CLASS_KEY] = multiplets_classes

            multiplet_adatas.append(multiplets_ad)

        joint_multiplets_ad = anndata.concat(multiplet_adatas)
        joint_multiplets_ad.uns[ONE_HOT_CLASS_KEY] = one_hot_class_names
        return joint_multiplets_ad

    def train(
        self,
        max_epochs: int = 400,
        lr: float = 1e-3,
        use_gpu: Optional[Union[str, int, bool]] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        plan_kwargs: Optional[dict] = None,
        early_stopping: bool = True,
        early_stopping_patience: int = 30,
        early_stopping_min_delta: float = 0.0,
        **kwargs,
    ):
        """
        Trains the model.

        Parameters
        ----------
        max_epochs
            Number of epochs to train for
        lr
            Learning rate for optimization.
        use_gpu
            Use default GPU if available (if None or True), or index of GPU to use (if int),
            or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        plan_kwargs
            Keyword args for :class:`~scvi.train.ClassifierTrainingPlan`. Keyword arguments passed to
        early_stopping
            Adds callback for early stopping on validation_loss
        early_stopping_patience
            Number of times early stopping metric can not improve over early_stopping_min_delta
        early_stopping_min_delta
            Threshold for counting an epoch torwards patience
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        update_dict = {
            "lr": lr,
            # "labels_key": ONE_HOT_CLASS_KEY,
        }
        if plan_kwargs is not None:
            plan_kwargs.update(update_dict)
        else:
            plan_kwargs = update_dict

        if early_stopping:
            early_stopping_callback = [
                LoudEarlyStopping(
                    monitor="validation_loss" if train_size != 1.0 else "train_loss",
                    min_delta=early_stopping_min_delta,
                    patience=early_stopping_patience,
                    mode="min",
                )
            ]
            if "callbacks" in kwargs:
                kwargs["callbacks"] += early_stopping_callback
            else:
                kwargs["callbacks"] = early_stopping_callback
            kwargs["check_val_every_n_epoch"] = 1

        if max_epochs is None:
            n_cells = self.adata.n_obs
            max_epochs = int(np.min([round((20000 / n_cells) * 400), 400]))

        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else dict()
        plan_kwargs["loss"] = torch.nn.BCEWithLogitsLoss

        # note that `.adata_manager` here is a [Cells+SyntheticCells, Latents] object
        # created in the `.from_scvi_model` method.
        logger.info("Creating DataSplitter")
        data_splitter = DataSplitter(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            use_gpu=use_gpu,
        )
        # classifier training plan will use `labels_key` as the labels to predict
        # `labels_key` is extracted from the `adata_manager` registry, rather than
        # directly querying any part of the `AnnData` object.
        logger.info("Creating TrainingPlan")
        training_plan = MultiBinaryClassifierTrainingPlan(self.module, **plan_kwargs)
        # use `BCEWithLogitsLoss`
        runner = TrainRunner(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            **kwargs,
        )
        return runner()

    @torch.inference_mode()
    def predict(
        self, soft: bool = True, include_simulated_multiplets: bool = False
    ) -> pd.DataFrame:
        """
        Return multiplet predictions.

        Parameters
        ----------
        soft
            Return probabilities instead of class label
        include_simulated_multiplets
            Return probabilities for simulated multiplets as well.

        Returns
        -------
        DataFrame with prediction, index corresponding to cell barcode.
        """
        adata = self._validate_anndata(None)

        scdl = self._make_data_loader(
            adata=adata,
        )

        @auto_move_data
        def auto_forward(module, x):
            return module(x)

        y_pred = []
        for _, tensors in enumerate(scdl):
            x = tensors[REGISTRY_KEYS.X_KEY]
            pred = auto_forward(self.module, x)
            y_pred.append(pred.cpu())

        y_pred = torch.cat(y_pred).numpy()

        label = self.adata.obs["_cpxvi_sim"].values.ravel()
        mask = label == "singlet" if not include_simulated_multiplets else slice(None)

        preds = y_pred[mask]

        cols = adata.uns[ONE_HOT_CLASS_KEY]
        preds_df = pd.DataFrame(preds, columns=cols, index=self.adata.obs_names[mask])

        if not soft:
            preds_df = preds_df.idxmax(axis=1)

        return preds_df

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        labels_key: Optional[str] = None,
        layer: Optional[str] = None,
        **kwargs,
    ):
        """
        Setup the AnnData of latent variable means used to train
        the classifier.

        %(summary)s.

        Parameters
        ----------
        %(param_labels_key)s
        %(param_layer)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=False),
            ObsmField(REGISTRY_KEYS.LABELS_KEY, labels_key),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)


def _validate_scvi_model(scvi_model: SCVI, restrict_to_batch: str):
    if scvi_model.summary_stats.n_batch > 1 and restrict_to_batch is None:
        warnings.warn(
            "Solo should only be trained on one lane of data using `restrict_to_batch`. Performance may suffer.",
            UserWarning,
        )
