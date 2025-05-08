from __future__ import annotations

import logging
from functools import partial
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pyro
from pyro.infer import Trace_ELBO

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data._utils import get_anndata_attribute
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LabelsWithUnlabeledObsField,
    LayerField,
    ObsmField,
)
from scvi.dataloaders import AnnTorchDataset
from scvi.model._utils import (
    scrna_raw_counts_properties,
)
from scvi.model.base import ArchesMixin, BaseModelClass, PyroSampleMixin, PyroSviTrainMixin
from scvi.model.base._de_core import _de_core
from scvi.utils import de_dsp, setup_anndata_dsp

from ._module import RESOLVAE
from ._utils import ResolVIPredictiveMixin

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from typing import Literal

    from anndata import AnnData

logger = logging.getLogger(__name__)


class RESOLVI(
    PyroSviTrainMixin, PyroSampleMixin, ResolVIPredictiveMixin, BaseModelClass, ArchesMixin
):
    """
    ResolVI addresses noise and bias in single-cell resolved spatial transcriptomics data.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.SCVI.setup_anndata`.
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
    **model_kwargs
        Keyword args for :class:`~scvi.module.VAE`

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

    1. :doc:`/tutorials/notebooks/api_overview`
    2. :doc:`/tutorials/notebooks/harmonization`
    3. :doc:`/tutorials/notebooks/scarches_scvi_tools`
    4. :doc:`/tutorials/notebooks/scvi_in_R`
    """

    _module_cls = RESOLVAE
    _block_parameters = []

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 32,
        n_hidden_encoder: int = 128,
        n_latent: int = 10,
        n_layers: int = 2,
        dropout_rate: float = 0.05,
        dispersion: Literal["gene", "gene-batch"] = "gene",
        gene_likelihood: Literal["nb", "poisson"] = "nb",
        background_ratio=None,
        median_distance=None,
        semisupervised=False,
        mixture_k=50,
        downsample_counts=True,
        **model_kwargs,
    ):
        pyro.clear_param_store()

        super().__init__(adata)
        if semisupervised:
            self._set_indices_and_labels()

        results = self.compute_dataset_dependent_priors()

        n_cats_per_cov = (
            self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )
        n_labels = self.summary_stats.n_labels - 1

        if background_ratio is None:
            background_ratio = results["background_ratio"]
        if median_distance is None:
            median_distance = results["median_distance"]
        if downsample_counts:
            downsample_counts_mean = results["mean_log_counts"]
            downsample_counts_std = results["std_log_counts"]
        else:
            downsample_counts_mean = None
            downsample_counts_std = 1.0

        expression_anntorchdata = AnnTorchDataset(
            self.adata_manager,
            getitem_tensors=["X"],
            load_sparse_tensor=True,
        )
        self.module = self._module_cls(
            n_input=self.summary_stats.n_vars,
            n_batch=self.summary_stats.n_batch,
            n_cats_per_cov=n_cats_per_cov,
            n_labels=n_labels,
            mixture_k=mixture_k,
            expression_anntorchdata=expression_anntorchdata,
            n_neighbors=self.summary_stats.n_distance_neighbor,
            n_obs=self.summary_stats["n_ind_x"],
            n_hidden=n_hidden,
            n_hidden_encoder=n_hidden_encoder,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            background_ratio=background_ratio,
            median_distance=median_distance,
            semisupervised=semisupervised,
            downsample_counts_mean=downsample_counts_mean,
            downsample_counts_std=downsample_counts_std,
            **model_kwargs,
        )
        self._model_summary_string = (
            f"RESOLVI Model with the following params: \nn_hidden: {n_hidden} "
            f"n_latent: {n_latent}, n_layers: {n_layers}, dropout_rate: "
            f"{dropout_rate}, dispersion: {dispersion}, gene_likelihood: {gene_likelihood} "
            f"n_neighbors: {self.summary_stats.n_distance_neighbor}"
        )
        self.init_params_ = self._get_init_params(locals())

    def train(
        self,
        max_epochs: int = 50,
        lr: float = 3e-3,
        lr_extra: float = 1e-2,
        extra_lr_parameters: tuple = ("per_neighbor_diffusion_map", "u_prior_means"),
        batch_size: int = 512,
        weight_decay: float = 0.0,
        eps: float = 1e-4,
        n_steps_kl_warmup: int | None = None,
        n_epochs_kl_warmup: int | None = 20,
        plan_kwargs: dict | None = None,
        expose_params: list = (),
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
        lr_extra
            Learning rate for parameters (non-amortized and custom ones)
        extra_lr_parameters
            List of parameters to train with `lr_extra` learning rate.
        batch_size
            Minibatch size to use during training.
        weight_decay
            weight decay regularization term for optimization
        eps
            Optimizer eps
        n_steps_kl_warmup
            Number of training steps (minibatches) to scale weight on KL divergences from 0 to 1.
            Only activated when `n_epochs_kl_warmup` is set to None. If `None`, defaults
            to `floor(0.75 * adata.n_obs)`.
        n_epochs_kl_warmup
            Number of epochs to scale weight on KL divergences from 0 to 1.
            Overrides `n_steps_kl_warmup` when both are not `None`.
        plan_kwargs
            Keyword args for :class:`~resolvi.train.PyroTrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        expose_params
            List of parameters to train if running model in Arches mode.
        **kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        blocked = self._block_parameters.copy()
        for name, param in self.module.named_parameters():
            if not param.requires_grad:
                blocked.append(name)
                param.requires_grad = True
        blocked = set(blocked) - set(expose_params)

        if blocked:
            print("Running scArches. Set lr to 0 and blocking variables.")

        def per_param_callable(module_name, param_name):
            store_name = f"{module_name}$$${param_name}" if "." in param_name else param_name
            if store_name in blocked:
                return {"lr": 0.0, "weight_decay": 0, "eps": eps}
            if store_name in extra_lr_parameters:
                return {"lr": lr_extra, "weight_decay": weight_decay, "eps": eps}
            else:
                return {"lr": lr, "weight_decay": weight_decay, "eps": eps}

        optim = pyro.optim.Adam(per_param_callable)

        if plan_kwargs is None:
            plan_kwargs = {}
        plan_kwargs.update(
            {
                "optim_kwargs": {"lr": lr, "weight_decay": weight_decay, "eps": eps},
                "optim": optim,
                "blocked": blocked,
                "n_epochs_kl_warmup": n_epochs_kl_warmup
                if n_epochs_kl_warmup is not None
                else max_epochs,
                "n_steps_kl_warmup": n_steps_kl_warmup,
                "loss_fn": Trace_ELBO(
                    num_particles=5, vectorize_particles=True, retain_graph=True
                ),
            }
        )

        super().train(
            max_epochs=max_epochs,
            train_size=1.0,
            plan_kwargs=plan_kwargs,
            batch_size=batch_size,
            **kwargs,
        )

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        batch_key: str | None = None,
        labels_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        prepare_data: bool | None = True,
        prepare_data_kwargs: dict = None,
        unlabeled_category: str = "unknown",
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        %(param_batch_key)s
        %(param_labels_key)s
        %(param_cat_cov_keys)s
        prepare_data
            If True, prepares AnnData for training. Computes spatial neighbors and distances.
        prepare_data_kwargs
            Keyword args for :meth:`scvi.external.RESOLVI._prepare_data`
        %(param_unlabeled_category)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        if layer is None:
            x = adata.X
        else:
            x = adata.layers[layer]
        assert np.min(x.sum(axis=1)) > 0, (
            "Please filter cells with less than 5 counts prior to running resolVI."
        )
        if batch_key is not None:
            adata.obs["_indices"] = (
                adata.obs[batch_key].astype(str) + "_" + adata.obs_names.astype(str)
            )
        else:
            adata.obs["_indices"] = adata.obs_names.astype(str)
        adata.obs["_indices"] = adata.obs["_indices"].astype("category")
        assert not adata.obs["_indices"].duplicated(keep="first").any(), (
            "obs_names need to be unique prior to running resolVI."
        )
        if labels_key is None:
            label_field = CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key)
        else:
            label_field = LabelsWithUnlabeledObsField(
                REGISTRY_KEYS.LABELS_KEY, labels_key, unlabeled_category
            )

        if prepare_data:
            if prepare_data_kwargs is None:
                prepare_data_kwargs = {}
            RESOLVI._prepare_data(adata, batch_key=batch_key, **prepare_data_kwargs)

        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            ObsmField("index_neighbor", "index_neighbor"),
            ObsmField("distance_neighbor", "distance_neighbor"),
            CategoricalObsField(REGISTRY_KEYS.INDICES_KEY, "_indices"),
            label_field,
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    @staticmethod
    def _prepare_data(
        adata, n_neighbors=10, spatial_rep="X_spatial", batch_key=None, slice_key=None, **kwargs
    ):
        if slice_key is not None:
            batch_key = slice_key
        try:
            import scanpy
            from sklearn.neighbors._base import _kneighbors_from_graph
        except ImportError as err:
            raise ImportError(
                "Please install scanpy and scikit-learn -- `pip install scanpy`"
            ) from err

        # Spatial neighbors are batch dependent.
        if batch_key is None:
            indices = [np.arange(adata.n_obs)]
        else:
            indices = [
                np.where(adata.obs[batch_key] == i)[0] for i in adata.obs[batch_key].unique()
            ]

        distance_neighbor = 1e6 * np.ones([adata.n_obs, n_neighbors])
        index_neighbor = np.zeros([adata.n_obs, n_neighbors], dtype=int)

        for index in indices:
            sub_data = adata[index].copy()
            try:
                import rapids_singlecell

                print("RAPIDS SingleCell is installed and can be imported")
                rapids_singlecell.pp.neighbors(
                    sub_data, n_neighbors=n_neighbors + 5, use_rep=spatial_rep
                )
            except ImportError:
                scanpy.pp.neighbors(sub_data, n_neighbors=n_neighbors + 5, use_rep=spatial_rep)
            distances = sub_data.obsp["distances"] ** 2

            distance_neighbor[index, :], index_neighbor_batch = _kneighbors_from_graph(
                distances, n_neighbors, return_distance=True
            )
            index_neighbor[index, :] = index[index_neighbor_batch]

        adata.obsm["X_spatial"] = adata.obsm[spatial_rep]
        adata.obsm["index_neighbor"] = index_neighbor
        adata.obsm["distance_neighbor"] = distance_neighbor

    def compute_dataset_dependent_priors(self, n_small_genes=None):
        x = self.adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)
        n_small_genes = x.shape[1] // 50 if n_small_genes is None else int(n_small_genes)
        # Computing library size over low expressed genes (expectation for background).
        # Handles sparse and dense counts.
        smallest_means = x[:, np.array(x.sum(0)).flatten().argsort()[:n_small_genes]].mean(
            1
        ) / np.array(x.mean(1))
        background_ratio = np.mean(np.array(smallest_means))

        # Median distance for empiric expectation of kernel size in diffusion
        distance = self.adata_manager.get_from_registry("distance_neighbor")
        median_distance = np.median(np.partition(distance, 5)[:, 5])
        log_library_size = np.log1p(np.array(x.sum(1)))
        mean_log_counts = np.median(log_library_size)
        std_log_counts = np.std(log_library_size)

        return {
            "background_ratio": background_ratio,
            "median_distance": median_distance,
            "mean_log_counts": mean_log_counts,
            "std_log_counts": std_log_counts,
        }

    @de_dsp.dedent
    def differential_expression(
        self,
        adata: AnnData | None = None,
        groupby: str | None = None,
        group1: Iterable[str] | None = None,
        group2: str | None = None,
        idx1: Sequence[int] | Sequence[bool] | None = None,
        idx2: Sequence[int] | Sequence[bool] | None = None,
        subset_idx: Sequence[int] | None = None,
        mode: Literal["vanilla", "change"] = "change",
        delta: float = 0.25,
        batch_size: int | None = None,
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: Iterable[str] | None = None,
        batchid2: Iterable[str] | None = None,
        fdr_target: float = 0.05,
        silent: bool = False,
        weights: Literal["uniform", "importance"] | None = "uniform",
        filter_outlier_cells: bool = False,
        **kwargs,
    ) -> pd.DataFrame:
        r"""A unified method for differential expression analysis.

        Implements `"vanilla"` DE :cite:p:`Lopez18` and `"change"` mode DE :cite:p:`Boyeau19`.

        Parameters
        ----------
        %(de_adata)s
        %(de_groupby)s
        %(de_group1)s
        %(de_group2)s
        %(de_idx1)s
        %(de_idx2)s
        %(de_subset_idx)s
        %(de_mode)s
        %(de_delta)s
        %(de_batch_size)s
        %(de_all_stats)s
        %(de_batch_correction)s
        %(de_batchid1)s
        %(de_batchid2)s
        %(de_fdr_target)s
        %(de_silent)s
        weights
        filter_outlier_cells
            Whether to filter outlier cells with
            :meth:`~scvi.model.base.DifferentialComputation.filter_outlier_cells`
        **kwargs
            Keyword args for :meth:`scvi.model.base.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential expression DataFrame.
        """
        adata = self._validate_anndata(adata)

        model_fn = partial(
            self.get_normalized_expression_importance,
            return_numpy=True,
            n_samples=5,
            batch_size=batch_size,
            weights=weights,
            return_mean=False,
        )

        representation_fn = self.get_latent_representation if filter_outlier_cells else None

        result = _de_core(
            adata_manager=self.get_anndata_manager(adata, required=True),
            model_fn=model_fn,
            representation_fn=representation_fn,
            groupby=groupby,
            group1=group1,
            group2=group2,
            idx1=idx1,
            idx2=idx2,
            subset_idx=subset_idx,
            all_stats=all_stats,
            all_stats_fn=scrna_raw_counts_properties,
            col_names=adata.var_names,
            mode=mode,
            batchid1=batchid1,
            batchid2=batchid2,
            delta=delta,
            batch_correction=batch_correction,
            fdr=fdr_target,
            silent=silent,
            **kwargs,
        )

        return result

    @de_dsp.dedent
    def differential_niche_abundance(
        self,
        adata: AnnData | None = None,
        groupby: str | None = None,
        group1: Iterable[str] | None = None,
        group2: str | None = None,
        neighbor_key: str | None = None,
        idx1: Sequence[int] | Sequence[bool] | None = None,
        idx2: Sequence[int] | Sequence[bool] | None = None,
        subset_idx: Sequence[int] | None = None,
        mode: Literal["vanilla", "change"] = "change",
        delta: float = 0.25,
        batch_size: int | None = None,
        fdr_target: float = 0.05,
        silent: bool = False,
        filter_outlier_cells: bool = False,
        pseudocounts: float = 1e-3,
        **kwargs,
    ) -> pd.DataFrame:
        r"""A unified method for niche differential abundance analysis.

        Implements `"vanilla"` DE :cite:p:`Lopez18` and `"change"` mode DE :cite:p:`Boyeau19`.

        Parameters
        ----------
        %(de_adata)s
        %(de_groupby)s
        %(de_group1)s
        %(de_group2)s
        neighbor_key
            Obsm key containing the spatial neighbors of each cell.
        %(de_idx1)s
        %(de_idx2)s
        %(de_subset_idx)s
        %(de_mode)s
        %(de_delta)s
        %(de_batch_size)s
        %(de_fdr_target)s
        %(de_silent)s
        filter_outlier_cells
            Whether to filter outlier cells with
            :meth:`~scvi.model.base.DifferentialComputation.filter_outlier_cells`
        pseudocounts
            pseudocount offset used for the mode `change`.
            When None, observations from non-expressed genes are used to estimate its value.
        **kwargs
            Keyword args for :meth:`scvi.model.base.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential expression DataFrame.
        """
        adata = self._validate_anndata(adata)

        model_fn = partial(
            self.get_neighbor_abundance,
            return_numpy=True,
            n_samples=5,
            batch_size=batch_size,
            return_mean=False,
            neighbor_key=neighbor_key,
        )

        representation_fn = self.get_latent_representation if filter_outlier_cells else None

        result = _de_core(
            adata_manager=self.get_anndata_manager(adata, required=True),
            model_fn=model_fn,
            representation_fn=representation_fn,
            groupby=groupby,
            group1=group1,
            group2=group2,
            idx1=idx1,
            idx2=idx2,
            subset_idx=subset_idx,
            all_stats=False,
            all_stats_fn=scrna_raw_counts_properties,
            col_names=self._label_mapping[:-1],
            mode=mode,
            batchid1=None,
            batchid2=None,
            delta=delta,
            batch_correction=False,
            fdr=fdr_target,
            silent=silent,
            pseudocounts=pseudocounts,
            **kwargs,
        )

        return result

    def predict(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        soft: bool = False,
        batch_size: int | None = 500,
        num_samples: int | None = 30,
    ) -> np.ndarray | pd.DataFrame:
        """
        Return cell label predictions.

        Parameters
        ----------
        adata
            AnnData object that has been registered via :meth:`~scvi.model.SCANVI.setup_anndata`.
        indices
            Subsample AnnData to these indices.
        soft
            If True, returns per class probabilities
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        num_samples
            Samples to draw from the posterior for cell-type prediction.
        """
        adata = self._validate_anndata(adata)

        if indices is None:
            indices = np.arange(adata.n_obs)

        sampled_prediction = self.sample_posterior(
            adata=adata,
            indices=indices,
            model=self.module.model_corrected,
            return_sites=["probs_prediction"],
            num_samples=num_samples,
            return_samples=False,
            batch_size=batch_size,
            summary_frequency=10,
            return_observed=True,
        )
        y_pred = sampled_prediction["post_sample_means"]["probs_prediction"]

        if not soft:
            y_pred = y_pred.argmax(axis=1)
            predictions = [self._code_to_label[p] for p in y_pred]
            return np.array(predictions)
        else:
            n_labels = len(y_pred[0])
            predictions = pd.DataFrame(
                y_pred,
                columns=self._label_mapping[:n_labels],
                index=adata.obs_names[indices],
            )
            return predictions

    def _set_indices_and_labels(self):
        """Set indices for labeled and unlabeled cells."""
        labels_state_registry = self.adata_manager.get_state_registry(REGISTRY_KEYS.LABELS_KEY)
        self.original_label_key = labels_state_registry.original_key
        self.unlabeled_category_ = labels_state_registry.unlabeled_category

        labels = get_anndata_attribute(
            self.adata,
            self.adata_manager.data_registry.labels.attr_name,
            self.original_label_key,
        ).ravel()
        self._label_mapping = labels_state_registry.categorical_mapping

        # set unlabeled and labeled indices
        self._unlabeled_indices = np.argwhere(labels == self.unlabeled_category_).ravel()
        self._labeled_indices = np.argwhere(labels != self.unlabeled_category_).ravel()
        self._code_to_label = dict(enumerate(self._label_mapping))
