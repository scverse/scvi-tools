"""scPoli model: population-level integration with cell-type prototypes."""

from __future__ import annotations

import logging
from typing import Literal

import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from lightning import Callback

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data._utils import _get_adata_minify_type, get_anndata_attribute
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LabelsWithUnlabeledObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
)
from scvi.model._utils import _init_library_size, get_max_epochs_heuristic
from scvi.model.base import (
    ArchesMixin,
    BaseMinifiedModeModelClass,
    EmbeddingMixin,
    UnsupervisedTrainingMixin,
)
from scvi.model.base._rnamixin import RNASeqMixin
from scvi.model.base._vaemixin import VAEMixin
from scvi.utils import setup_anndata_dsp

from ._module import ScPoliVAE

logger = logging.getLogger(__name__)


class ScPoliPrototypeCallback(Callback):
    """Manage prototype initialisation and epoch-end updates.

    Handles both labeled and unlabeled prototypes:

    - **Labeled**: At epoch ``pretrain_epochs`` (the phase boundary), each
      labeled prototype is seeded from the per-cell-type mean of the current
      encoder output (training cells only).  At the end of every subsequent
      epoch the means are recomputed — no gradient, matching the non-gradient
      update rule of the original scArches implementation.
    - **Unlabeled**: At the same phase boundary, training-cell latents are
      clustered (KMeans or Leiden) to find unlabeled prototype positions.
      After each epoch a dedicated Adam optimizer takes one gradient step on
      the prototype positions while keeping the encoder fixed (alternating
      optimization, equivalent to the approach in :cite:p:`Lotfollahi23`).

    Parameters
    ----------
    scpoli_model
        The parent :class:`~scvi.external.ScPoli` model instance.
    pretrain_epochs
        Number of ELBO-only warm-up epochs before prototype loss is activated.
    unlabeled_prototype_training
        Whether to fit unlabeled prototypes via clustering + gradient updates.
    clustering
        Clustering algorithm for unlabeled prototype initialisation.
        ``"leiden"`` (default) or ``"kmeans"``.
    clustering_res
        Leiden resolution.  Larger values produce finer clusters.
    n_clusters
        Number of KMeans clusters.  Required when ``clustering="kmeans"``;
        ignored for Leiden (cluster count is determined automatically).
    lr
        Learning rate for the unlabeled prototype Adam optimizer.
    eps
        Epsilon for the prototype optimizer.
    weight_decay
        Weight decay for the prototype optimizer.
    """

    def __init__(
        self,
        scpoli_model: ScPoli,
        pretrain_epochs: int,
        unlabeled_prototype_training: bool = False,
        clustering: str = "leiden",
        clustering_res: float = 2.0,
        n_clusters: int | None = None,
        lr: float = 1e-3,
        eps: float = 0.01,
        weight_decay: float = 0.04,
    ) -> None:
        super().__init__()
        self.scpoli_model = scpoli_model
        self.pretrain_epochs = pretrain_epochs
        self.unlabeled_prototype_training = unlabeled_prototype_training
        self.clustering = clustering
        self.clustering_res = clustering_res
        self.n_clusters = n_clusters
        self._lr = lr
        self._eps = eps
        self._weight_decay = weight_decay
        self.prototype_optim: torch.optim.Optimizer | None = None

    def _device(self) -> torch.device:
        """Device of the module (safe regardless of prototype state)."""
        return next(self.scpoli_model.module.parameters()).device

    def _get_train_indices(self, trainer) -> np.ndarray | None:
        """Return training-split cell indices from the Lightning datamodule."""
        dm = getattr(trainer, "datamodule", None)
        if dm is not None and hasattr(dm, "train_idx"):
            return np.asarray(dm.train_idx)
        return None  # fall back to all cells

    def _get_latent_and_labels(
        self, indices: np.ndarray | None = None
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Return latent z and integer label codes.

        Parameters
        ----------
        indices
            Subset of cell indices (e.g. training split).  ``None`` means all
            cells — used as a fallback when the datamodule is not available.
        """
        model = self.scpoli_model
        module = model.module

        # Temporarily mark as trained so get_latent_representation() passes its
        # _check_if_trained() guard.  Restored in the finally block so that
        # _update_history() in the TrainingRunner still sees is_trained_=False
        # and correctly initialises model.history_ (rather than trying to merge
        # into the still-None history_ dict and silently returning).
        was_trained = model.is_trained_
        was_training = module.training
        model.is_trained_ = True
        module.eval()
        try:
            with torch.inference_mode():
                latent_np = model.get_latent_representation(indices=indices)
        finally:
            model.is_trained_ = was_trained
            module.train(was_training)

        labels_raw = model.adata_manager.get_from_registry(REGISTRY_KEYS.LABELS_KEY)
        if hasattr(labels_raw, "values"):
            labels_raw = labels_raw.values
        labels_np = np.asarray(labels_raw).ravel().astype(np.int64)
        if indices is not None:
            labels_np = labels_np[indices]

        device = self._device()
        latent_t = torch.tensor(latent_np, dtype=torch.float32, device=device)
        labels_t = torch.tensor(labels_np, dtype=torch.long, device=device)
        return latent_t, labels_t

    def _initialize_unlabeled_prototypes(self, latent: torch.Tensor) -> None:
        """Cluster training-cell latents and seed unlabeled prototypes.

        Only runs when ``unlabeled_prototype_training=True`` and the module's
        ``unlabeled_weight > 0``.
        """
        module = self.scpoli_model.module

        if not self.unlabeled_prototype_training:
            return

        lat_np = latent.detach().cpu().numpy()
        device = self._device()

        if self.clustering == "kmeans" and self.n_clusters is not None:
            logger.info(
                f"Initializing unlabeled prototypes with KMeans (n_clusters={self.n_clusters})."
            )
            from sklearn.cluster import KMeans

            km = KMeans(n_clusters=self.n_clusters, n_init="auto").fit(lat_np)
            centers = torch.tensor(km.cluster_centers_, dtype=torch.float32, device=device)
        else:
            # Default to Leiden (auto-determines cluster count)
            if self.clustering == "kmeans":
                logger.info("n_clusters not provided; falling back to Leiden clustering.")
            else:
                logger.info(
                    f"Initializing unlabeled prototypes with Leiden "
                    f"(resolution={self.clustering_res})."
                )
            import scanpy as sc

            lat_adata = sc.AnnData(lat_np)
            sc.pp.neighbors(lat_adata)
            sc.tl.leiden(lat_adata, resolution=self.clustering_res)

            cluster_ids = np.asarray(lat_adata.obs["leiden"], dtype=int)
            n_clusters = cluster_ids.max() + 1
            cluster_centers = np.stack(
                [lat_np[cluster_ids == k].mean(axis=0) for k in range(n_clusters)]
            )
            self.n_clusters = n_clusters
            logger.info(f"Leiden clustering found {self.n_clusters} unlabeled clusters.")
            centers = torch.tensor(cluster_centers, dtype=torch.float32, device=device)

        # Allocate the buffer (replaces the empty (0, n_latent) placeholder)
        module.set_unlabeled_prototypes(centers)

        # Create the prototype optimizer AFTER set_unlabeled_prototypes so
        # it holds a reference to the freshly allocated tensor.
        self.prototype_optim = torch.optim.Adam(
            [module.prototypes_unlabeled],
            lr=self._lr,
            eps=self._eps,
            weight_decay=self._weight_decay,
        )
        logger.info(f"Unlabeled prototype optimizer created (lr={self._lr}, eps={self._eps}).")

    def _update_unlabeled_prototypes(self, latent: torch.Tensor) -> None:
        """Take one Adam step on unlabeled prototype positions.

        The ``latent`` tensor must already be detached from the computation
        graph (i.e. produced by :meth:`get_latent_representation`) so that
        gradients flow only to the prototype positions and not back through
        the encoder.  This mirrors the alternating optimization in scArches.
        """
        if self.prototype_optim is None:
            return

        module = self.scpoli_model.module
        prototypes = module.prototypes_unlabeled

        # Temporarily enable gradients on the prototype tensor so autograd
        # can compute ∂loss/∂prototypes.  The loss is recomputed from scratch
        # each iteration so that (a) the graph is fresh after each .backward()
        # call, and (b) each step uses the updated prototype positions.
        prototypes.requires_grad_(True)
        for _ in range(10):
            dists = torch.cdist(latent, prototypes, p=2)
            min_dist, y_hat = torch.min(dists, dim=1)
            unique_clusters = y_hat.unique()
            loss = torch.stack([min_dist[y_hat == c].mean() for c in unique_clusters]).mean()
            self.prototype_optim.zero_grad()
            loss.backward()
            self.prototype_optim.step()

        # Re-disable gradients so that the buffer is treated as a constant
        # during the next mini-batch forward pass.
        prototypes.requires_grad_(False)

    def on_train_epoch_start(self, trainer, pl_module) -> None:
        """Initialise prototypes at the phase boundary."""
        if getattr(self.scpoli_model, "is_query_model_", False):
            return
        if trainer.current_epoch != self.pretrain_epochs:
            return

        train_idx = self._get_train_indices(trainer)
        latent, labels = self._get_latent_and_labels(indices=train_idx)

        # Labeled prototypes — skip only if already seeded (e.g. loaded from a
        # reference model via load_query_data).
        if self.scpoli_model.module._prototypes_initialized.item():
            logger.info(
                f"Epoch {trainer.current_epoch}: labeled prototypes already "
                "initialised — skipping reinitialisation."
            )
        else:
            logger.info(
                f"Epoch {trainer.current_epoch}: initialising labeled prototypes "
                "from training-cell encoder means."
            )
            self.scpoli_model.module.initialize_prototypes(latent, labels)

        # Unlabeled prototypes — independent check so re-training with a
        # non-zero unlabeled_weight still seeds them even when labeled
        # prototypes were already present.
        if self.scpoli_model.module._unlabeled_prototypes_initialized.item():
            logger.info(
                f"Epoch {trainer.current_epoch}: unlabeled prototypes already "
                "initialised — skipping reinitialisation."
            )
        else:
            self._initialize_unlabeled_prototypes(latent)

    def on_train_epoch_end(self, trainer, pl_module) -> None:
        """Update prototypes at the end of every prototype-training epoch."""
        if getattr(self.scpoli_model, "is_query_model_", False):
            return
        if trainer.current_epoch >= self.pretrain_epochs:
            train_idx = self._get_train_indices(trainer)
            latent, labels = self._get_latent_and_labels(indices=train_idx)
            # Labeled prototype update: non-gradient mean reassignment
            self.scpoli_model.module.update_prototypes(latent, labels)
            # Unlabeled prototype update: one Adam gradient step
            self._update_unlabeled_prototypes(latent.detach())


class ScPoli(
    RNASeqMixin,
    VAEMixin,
    ArchesMixin,
    EmbeddingMixin,
    UnsupervisedTrainingMixin,
    BaseMinifiedModeModelClass,
):
    r"""Population-level integration with prototype-based cell-type mapping.

    .. note::

        Storage keys used when minifying the adata (written to ``obsm`` / ``obs``):

        - ``_LATENT_QZM_KEY = "_scpoli_latent_qzm"``
        - ``_LATENT_QZV_KEY = "_scpoli_latent_qzv"``
        - ``_OBSERVED_LIB_SIZE_KEY`` is inherited from :class:`~scvi.model.base.BaseModelClass`
          (``"observed_lib_size"``).


    scPoli :cite:p:`Lotfollahi23` extends a standard VAE (:cite:p:`Lopez18`)
    with cell-type prototypes in latent space.  Training proceeds in two
    phases managed by a single training loop:

    1. **Pre-training** (first ``pretrain_epochs`` epochs, default 90 % of
       total): the VAE is trained with the standard ELBO only.  No prototype
       loss is applied.
    2. **Prototype training** (remaining epochs): at the epoch boundary,
       labeled prototypes are initialised as per-cell-type latent means of
       training cells, and unlabeled prototypes are seeded by clustering.
       After each subsequent epoch the labeled prototypes are updated to the
       new latent means (no gradient), and unlabeled prototype positions take
       one Adam step (gradient only on the centroids).

    After reference training, new query batches can be mapped via
    :meth:`~scvi.model.base.ArchesMixin.load_query_data`, which freezes
    reference encoder/decoder weights and only updates batch embeddings.

    Parameters
    ----------
    adata
        AnnData object registered via :meth:`~scvi.external.ScPoli.setup_anndata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers in encoder and decoder.
    dropout_rate
        Dropout rate applied to encoder layers.
    gene_likelihood
        Reconstruction distribution. One of ``"nb"`` (default), ``"zinb"``,
        or ``"poisson"``.
    eta
        Weight :math:`\\eta` on the labeled prototype loss. Matches the ``eta``
        parameter in the original scArches scPoli implementation.
    batch_embedding_dim
        Dimensionality of the :class:`~torch.nn.Embedding` layer used to
        embed batch labels before passing them to the encoder.  Set to ``0``
        to disable embedding and use standard one-hot encoding instead.
        The default of ``10`` matches the scArches conditioning approach.
    unlabeled_weight
        Weight on the unlabeled prototype loss.  ``0.0`` (default) disables
        it entirely.  A common non-zero value used in scArches is ``0.01``.
    **model_kwargs
        Additional keyword arguments forwarded to
        :class:`~scvi.external.scpoli.ScPoliVAE`.

    Examples
    --------
    Reference model:

    >>> import scvi
    >>> scvi.external.ScPoli.setup_anndata(
    ...     adata_ref,
    ...     labels_key="cell_type",
    ...     unlabeled_category="Unknown",
    ...     batch_key="batch",
    ... )
    >>> model = scvi.external.ScPoli(adata_ref)
    >>> model.train(max_epochs=400)  # 360 pretrain + 40 prototype epochs

    Query mapping:

    >>> scvi.external.ScPoli.setup_anndata(
    ...     adata_query,
    ...     labels_key="cell_type",
    ...     unlabeled_category="Unknown",
    ...     batch_key="batch",
    ... )
    >>> query_model = scvi.external.ScPoli.load_query_data(adata_query, model)
    >>> query_model.train(max_epochs=50, pretrain_epochs=0)

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/scrna/scpoli`
    """

    _module_cls = ScPoliVAE
    _LATENT_QZM_KEY = "_scpoli_latent_qzm"
    _LATENT_QZV_KEY = "_scpoli_latent_qzv"

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 2,
        dropout_rate: float = 0.0,
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "nb",
        batch_embedding_dim: int = 10,
        **model_kwargs,
    ):
        super().__init__(adata)
        self._set_indices_and_labels()
        n_labels = self.summary_stats.n_labels - 1

        n_cats_per_cov = (
            self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )
        n_batch = self.summary_stats.n_batch

        use_size_factor_key = REGISTRY_KEYS.SIZE_FACTOR_KEY in self.adata_manager.data_registry
        library_log_means, library_log_vars = None, None
        if not use_size_factor_key:
            library_log_means, library_log_vars = _init_library_size(self.adata_manager, n_batch)

        _batch_embedding_kwargs: dict | None = None
        if batch_embedding_dim > 0 and n_batch > 1:
            _batch_embedding_kwargs = {"embedding_dim": batch_embedding_dim}

        self.module = ScPoliVAE(
            n_input=self.summary_stats.n_vars,
            n_batch=n_batch,
            n_labels=n_labels,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
            n_cats_per_cov=n_cats_per_cov,
            dropout_rate=dropout_rate,
            gene_likelihood=gene_likelihood,
            use_size_factor_key=use_size_factor_key,
            library_log_means=library_log_means,
            library_log_vars=library_log_vars,
            batch_embedding_kwargs=_batch_embedding_kwargs,
            **model_kwargs,
        )
        self._model_summary_string = (
            f"ScPoli Model\n"
            f"n_vars: {self.summary_stats.n_vars}, "
            f"n_batch: {n_batch}, "
            f"n_labels (known): {n_labels}\n"
            f"n_hidden: {n_hidden}, n_latent: {n_latent}, n_layers: {n_layers}\n"
            f"batch_embedding_dim: {batch_embedding_dim if n_batch > 1 else 'n/a (single batch)'}"
        )
        self.init_params_ = self._get_init_params(locals())
        self.n_labels = n_labels
        self.module.minified_data_type = self.minified_data_type

    def train(
        self,
        max_epochs: int | None = None,
        pretrain_epochs: int | None = None,
        unlabeled_weight: float = 0.01,
        n_epochs_kl_warmup: int = 100,
        eta: float = 0.0,
        clustering: str = "leiden",
        clustering_res: float = 2.0,
        n_clusters: int | None = None,
        plan_kwargs: dict | None = None,
        **kwargs,
    ):
        """Train the model in two phases inside a single training loop.

        **Phase 1 — ELBO pre-training** (epochs ``0 … pretrain_epochs - 1``):
        The VAE is optimised with the standard ELBO; prototypes are not used
        and add zero computational overhead.

        **Phase 2 — Prototype training** (epochs ``pretrain_epochs … max_epochs - 1``):
        At the start of epoch ``pretrain_epochs``,
        :meth:`~scvi.external.scpoli.ScPoliVAE.initialize_prototypes` seeds each
        labeled prototype from the per-cell-type mean of the current encoder
        output over **training cells only**.  Unlabeled cells are simultaneously
        clustered to initialise unlabeled prototypes (when
        ``unlabeled_prototype_training=True``.
        At the end of every subsequent epoch:

        - Labeled prototypes are recomputed as per-cell-type means (no gradient).
        - Unlabeled prototype positions take one Adam gradient step that minimises
          the mean nearest-centroid distance (soft k-means objective).

        Parameters
        ----------
        max_epochs
            Total number of training epochs. If ``None``, determined by the
            dataset-size heuristic.
        pretrain_epochs
            Number of ELBO-only warm-up epochs. If ``None``, defaults to
            ``floor(max_epochs * 0.9)`` — matching the scArches default.
            Set to ``0`` to skip pre-training (e.g. for query fine-tuning).
        n_epochs_kl_warmup
            Number of epochs over which the KL weight is annealed from 0 to 1.
            Defaults to ``100``, matching the scArches ``alpha_epoch_anneal``
            default.
        eta
            Weight eta on the labeled prototype loss.
        unlabeled_weight
            Weight on unlabeled prototypes.
        clustering
            Clustering algorithm for unlabeled prototype initialisation.
            ``"leiden"`` (default) or ``"kmeans"``.
        clustering_res
            Leiden resolution parameter.
        n_clusters
            Number of clusters for KMeans.  Required if
            ``clustering="kmeans"``.
        plan_kwargs
            Keyword arguments for :class:`~scvi.train.TrainingPlan`.
            ScPoli sets ``weight_decay=0.04`` and ``eps=0.01`` by default
            (matching scArches); pass ``plan_kwargs={"weight_decay": 0.0}``
            for query fine-tuning to prevent reference embedding drift.
        **kwargs
            Additional keyword arguments forwarded to
            :meth:`~scvi.model.base.UnsupervisedTrainingMixin.train`
            (``accelerator``, ``devices``, ``batch_size``, ``train_size``,
            ``early_stopping``, ``datasplitter_kwargs``, etc.).
        """
        if max_epochs is None:
            max_epochs = get_max_epochs_heuristic(self.adata.n_obs)

        if pretrain_epochs is None:
            pretrain_epochs = int(np.floor(max_epochs * 0.9))

        logger.info(
            f"Training ScPoli for {max_epochs} epochs "
            f"({pretrain_epochs} pretrain, "
            f"{max_epochs - pretrain_epochs} prototype)."
        )

        # scArches optimizer defaults; user plan_kwargs take priority.
        _plan_kwargs: dict = {
            "weight_decay": 0.04,
            "eps": 0.01,
            "n_epochs_kl_warmup": n_epochs_kl_warmup,
            "eta": eta,
            "unlabeled_weight": unlabeled_weight,
        }
        _plan_kwargs.update(plan_kwargs or {})

        proto_callback = ScPoliPrototypeCallback(
            scpoli_model=self,
            pretrain_epochs=pretrain_epochs,
            unlabeled_prototype_training=unlabeled_weight > 0,
            clustering=clustering,
            clustering_res=clustering_res,
            n_clusters=n_clusters,
        )
        callbacks = list(kwargs.pop("callbacks", None) or [])
        callbacks.append(proto_callback)

        super().train(
            max_epochs=max_epochs,
            plan_kwargs=_plan_kwargs,
            callbacks=callbacks,
            **kwargs,
        )

    def get_prototypes(self) -> pd.DataFrame:
        """Return labeled prototype embeddings as a :class:`~pandas.DataFrame`.

        Returns
        -------
        :class:`~pandas.DataFrame` of shape ``(n_labels, n_latent)`` with cell-type
        names as the index and latent dimension indices as columns.  Raises
        ``RuntimeError`` if prototypes have not yet been initialised.
        """
        if self.module.n_prototypes == 0:
            # No labeled cell types — return an empty DataFrame.
            return pd.DataFrame(columns=[f"dim_{i}" for i in range(self.module._n_latent)])
        if not self.module._prototypes_initialized.item():
            raise RuntimeError(
                "Prototypes have not been initialised yet. "
                "Train the model first (or set pretrain_epochs < max_epochs)."
            )
        with torch.inference_mode():
            protos = self.module.prototypes_labeled.detach().cpu().numpy()
        label_names = [self._code_to_label[i] for i in range(self.n_labels)]
        return pd.DataFrame(
            protos,
            index=label_names,
            columns=[f"dim_{i}" for i in range(protos.shape[1])],
        )

    def get_unlabeled_prototypes(self) -> np.ndarray:
        """Return unlabeled prototype embeddings as a numpy array.

        Returns
        -------
        Array of shape ``(n_clusters, n_latent)``.  Raises ``RuntimeError``
        if unlabeled prototypes have not been initialised (e.g. because
        ``unlabeled_weight=0`` or ``unlabeled_prototype_training=False``).
        """
        if not self.module._unlabeled_prototypes_initialized.item():
            raise RuntimeError(
                "Unlabeled prototypes have not been initialised. "
                "Set unlabeled_weight > 0 and train the model."
            )
        with torch.inference_mode():
            return self.module.prototypes_unlabeled.detach().cpu().numpy()

    def _set_indices_and_labels(self, datamodule=None):
        """Set indices for labeled and unlabeled cells."""
        labels_state_registry = self.get_state_registry(REGISTRY_KEYS.LABELS_KEY)
        self.original_label_key = labels_state_registry.original_key
        self.unlabeled_category_ = labels_state_registry.unlabeled_category

        if datamodule is None:
            self.labels_ = get_anndata_attribute(
                self.adata,
                self.adata_manager.data_registry.labels.attr_name,
                self.original_label_key,
                mod_key=getattr(self.adata_manager.data_registry.labels, "mod_key", None),
            ).ravel()
        else:
            if datamodule.registry["setup_method_name"] == "setup_datamodule":
                self.labels_ = datamodule.labels_.ravel()
            else:
                self.labels_ = datamodule.labels.ravel()
        self._label_mapping = labels_state_registry.categorical_mapping

        # set unlabeled and labeled indices
        self._unlabeled_indices = np.argwhere(self.labels_ == self.unlabeled_category_).ravel()
        self._labeled_indices = np.argwhere(self.labels_ != self.unlabeled_category_).ravel()
        self._code_to_label = dict(enumerate(self._label_mapping))

    @torch.inference_mode()
    def predict(
        self,
        adata: AnnData | None = None,
        indices: list[int] | None = None,
        soft: bool = False,
        batch_size: int = 256,
    ) -> np.ndarray | pd.DataFrame:
        """Predict cell-type labels by nearest labeled prototype.

        Assigns each cell to the cell type whose prototype is closest in L2
        distance to the cell's latent representation, matching the
        ``classify`` method of the original scArches scPoli.

        Parameters
        ----------
        adata
            AnnData to predict on. Uses the model's ``adata`` if ``None``.
        indices
            Optional integer subset of cells.
        soft
            If ``True``, return a :class:`~pandas.DataFrame` of softmax
            proximity scores.  If ``False`` (default), return string labels.
        batch_size
            Minibatch size for computing latent representations.

        Returns
        -------
        If ``soft=False``:
            :class:`~numpy.ndarray` of shape ``(n_obs,)`` with string labels.
        If ``soft=True``:
            :class:`~pandas.DataFrame` of shape ``(n_obs, n_labels)`` with
            per-cell-type softmax proximity scores.
        """
        adata = self._validate_anndata(adata)
        latent = self.get_latent_representation(
            adata=adata, indices=indices, batch_size=batch_size
        )
        prototypes = self.get_prototypes().values  # (n_prototypes, n_latent) numpy array

        # L2 distances — consistent with torch.cdist(p=2) used in loss
        sq_dist = np.sum(
            (latent[:, np.newaxis, :] - prototypes[np.newaxis, :, :]) ** 2, axis=-1
        )  # (n_obs, n_prototypes)

        if soft:
            neg_dist = -np.sqrt(np.maximum(sq_dist, 0.0))
            neg_dist -= neg_dist.max(axis=1, keepdims=True)
            exp_neg = np.exp(neg_dist)
            probs = exp_neg / exp_neg.sum(axis=1, keepdims=True)
            label_names = [
                self._code_to_label[i]
                for i in range(self.n_labels)
                if self._code_to_label[i] != self.unlabeled_category_
            ]
            obs_names = adata.obs_names if indices is None else adata.obs_names[indices]
            return pd.DataFrame(probs, index=obs_names, columns=label_names)

        pred_indices = np.argmin(sq_dist, axis=1)
        return np.array([self._code_to_label[i] for i in pred_indices])

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        labels_key: str,
        unlabeled_category: str,
        layer: str | None = None,
        batch_key: str | None = None,
        size_factor_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_labels_key)s
        %(param_unlabeled_category)s
        %(param_layer)s
        %(param_batch_key)s
        %(param_size_factor_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            LabelsWithUnlabeledObsField(REGISTRY_KEYS.LABELS_KEY, labels_key, unlabeled_category),
            NumericalObsField(REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
        ]
        adata_minify_type = _get_adata_minify_type(adata)
        if adata_minify_type is not None:
            anndata_fields += cls._get_fields_for_adata_minification(adata_minify_type)
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
