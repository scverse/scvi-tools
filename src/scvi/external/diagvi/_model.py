"""DIAGVI model for multi-modal integration with guidance graphs."""

from __future__ import annotations

import logging
import os
import warnings
from collections.abc import Sequence
from itertools import cycle
from typing import TYPE_CHECKING, Literal

import numpy as np
import scipy.sparse
import torch
from mudata import MuData
from torch.utils.data import DataLoader

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data._constants import _MODEL_NAME_KEY, _SETUP_ARGS_KEY
from scvi.data.fields import CategoricalObsField, LabelsWithUnlabeledObsField, LayerField
from scvi.dataloaders import DataSplitter
from scvi.external.diagvi._utils import (
    _check_guidance_graph_consistency,
    _construct_guidance_graph,
    _load_saved_diagvi_files,
)
from scvi.model._utils import get_max_epochs_heuristic, parse_device_args, use_distributed_sampler
from scvi.model.base import BaseModelClass, VAEMixin
from scvi.module._constants import MODULE_KEYS
from scvi.train import Trainer
from scvi.utils import dependencies, setup_anndata_dsp
from scvi.utils._docstrings import devices_dsp

from ._module import DIAGVAE
from ._task import DiagTrainingPlan

if TYPE_CHECKING:
    from typing import Literal

    import pandas as pd
    from anndata import AnnData
    from torch_geometric.data import Data

    from scvi.dataloaders import AnnDataLoader

logger = logging.getLogger(__name__)


class DIAGVI(BaseModelClass, VAEMixin):
    """Diagonal Multi-Modal Integration Variational Inference (DIAGVI) model.

    Integrates multi-modal single-cell data using a guidance graph and supports
    semi-supervised learning and GMM priors.

    The model architecture is inspired by GLUE (Cao & Gao, 2022).
    This scvi-tools implementation is based on GimVI.
    Handling of continuous data in decoder is inspired by CytoVI.

    Parameters
    ----------
    adatas
        Dictionary mapping input names to AnnData objects.
    guidance_graph
        Precomputed guidance graph. If None, it will be constructed from the data
        by using overlapping feature names.
    mapping_df
        DataFrame specifying feature correspondences between modalities
        (used to compute the guidance graph).
    n_latent
        Dimensionality of the latent space.
    n_hidden
        Number of nodes per hidden layer.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    dropout_rate
        Dropout rate for neural networks.
    **model_kwargs
        Additional keyword arguments for :class:`~scvi.external.diagvi._module.DIAGVAE`.

    Examples
    --------
    >>> adatas = {"rna_data": adata_rna, "protein_data": adata_protein}
    >>> model = DIAGVI(adatas)
    >>> model.train()
    """

    _module_cls = DIAGVAE
    _data_splitter_cls = DataSplitter
    _training_plan_cls = DiagTrainingPlan

    def __init__(
        self,
        adatas: dict[str, AnnData] | MuData,
        guidance_graph: Data | None = None,
        mapping_df: pd.DataFrame | None = None,
        n_latent: int = 50,
        n_hidden: int = 256,
        n_layers: int = 2,
        dropout_rate: float = 0.1,
        **model_kwargs,
    ):
        super().__init__()
        # Handle MuData input by extracting modalities as dict
        if isinstance(adatas, MuData):
            self.adatas = {mod_key: adatas.mod[mod_key] for mod_key in adatas.mod.keys()}
        else:
            self.adatas = adatas
        self.input_names = list(self.adatas.keys())
        self.adata_managers = {
            name: self._get_most_recent_anndata_manager(adata, required=True)
            for name, adata in self.adatas.items()
        }

        self.registries_ = []
        for adm in self.adata_managers.values():
            self._register_manager_for_instance(adm)
            self.registries_.append(adm.registry)
        
        sum_stats = {name: adm.summary_stats for name, adm in self.adata_managers.items()}
        self.summary_stats = sum_stats
        n_inputs = {name: s["n_vars"] for name, s in sum_stats.items()}
        n_batches = {name: s["n_batch"] for name, s in sum_stats.items()}
        n_labels = {name: s["n_labels"] for name, s in sum_stats.items()}

        generative_distributions = {
            name: adata.uns["diagvi_likelihood"] for name, adata in self.adatas.items()
        }
        normalize_lib = {
            name: adata.uns["diagvi_normalize_lib"] for name, adata in self.adatas.items()
        }
        gmm_priors = {name: adata.uns["diagvi_gmm_prior"] for name, adata in self.adatas.items()}
        semi_supervised = {
            name: adata.uns["diagvi_semi_supervised"] for name, adata in self.adatas.items()
        }
        n_mixture_components = {
            name: adata.uns["diagvi_n_mixture_components"] for name, adata in self.adatas.items()
        }

        # Compute guidance graph if not provided, do a sanity check
        if guidance_graph is not None:
            self.guidance_graph = guidance_graph
        else:
            self.guidance_graph = _construct_guidance_graph(self.adatas, mapping_df)
        _check_guidance_graph_consistency(self.guidance_graph, self.adatas)
        self.module = self._module_cls(
            n_inputs=n_inputs,
            n_batches=n_batches,
            n_labels=n_labels,
            modality_likelihoods=generative_distributions,
            normalize_lib=normalize_lib,
            guidance_graph=self.guidance_graph,
            use_gmm_prior=gmm_priors,
            semi_supervised=semi_supervised,
            n_mixture_components=n_mixture_components,
            n_latent=n_latent,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            **model_kwargs,
        )

        self._model_summary_string = (
            f"DiagVI Model with the following params: input names: {self.input_names}, "
            f"n_inputs: {n_inputs}, n_batches: {n_batches}, "
            f"generative distributions: {generative_distributions}, "
            f" n_latent: {n_latent}, n_hidden: {n_hidden}, n_layers: {n_layers}."
        )

        self.init_params_ = self._get_init_params(locals())
        logger.info(self._model_summary_string)

    @devices_dsp.dedent
    def train(
        self,
        max_epochs: int | None = None,
        batch_size: int = 1024,
        train_size: float = 0.9,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        shuffle_set_split: bool = True,
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        **kwargs,
    ):
        """Train the DIAGVI model.

        Parameters
        ----------
        max_epochs
            Maximum number of training epochs. If None, a heuristic is used.
        %(param_accelerator)s
        %(param_devices)s
        train_size
            Proportion of data to use for training (rest for validation).
        shuffle_set_split
            Whether to shuffle data before splitting into train/validation.
        datasplitter_kwargs
            Additional keyword arguments for the DataSplitter and DataLoaders.
        plan_kwargs
            Additional keyword arguments for the training plan.
        **kwargs
            Additional keyword arguments for the Trainer.
        """
        # Compute max_epochs if not provided
        if max_epochs is None:
            min_obs = np.min(
                [
                    self.summary_stats[self.input_names[0]]["n_cells"],
                    self.summary_stats[self.input_names[1]]["n_cells"],
                ]
            )
            max_epochs = get_max_epochs_heuristic(min_obs)
            logger.info(f"max_epochs was approximated to {max_epochs}")
        
        accelerator, devices, device = parse_device_args(
            accelerator=accelerator,
            devices=devices,
            return_device="torch",
        )

        # Extract early stopping parameters from kwargs with defaults
        early_stopping = kwargs.pop("early_stopping", True)
        early_stopping_monitor = kwargs.pop("early_stopping_monitor", "validation_loss")
        early_stopping_min_delta = kwargs.pop("early_stopping_min_delta", 1)
        early_stopping_patience = kwargs.pop("early_stopping_patience", 10)
        early_stopping_mode = kwargs.pop("early_stopping_mode", "min")

        # Initialize Trainer
        self.trainer = Trainer(
            max_epochs=max_epochs,
            enable_progress_bar=True,
            progress_bar_refresh_rate=1,
            early_stopping=early_stopping,
            early_stopping_monitor=early_stopping_monitor,
            early_stopping_min_delta=early_stopping_min_delta,
            early_stopping_patience=early_stopping_patience,
            early_stopping_mode=early_stopping_mode,
            accelerator=accelerator,
            devices=devices,
            **kwargs,
        )

        # Set up DataLoaders with DataSplitter
        validation_size = 1 - train_size
        self.train_indices_, self.test_indices_, self.validation_indices_ = [], [], []
        train_dls, test_dls, val_dls = {}, {}, {}
        datasplitter_kwargs = datasplitter_kwargs if isinstance(datasplitter_kwargs, dict) else {}
        for name, adm in self.adata_managers.items():
            ds = self._data_splitter_cls(
                adm,
                train_size=train_size,
                validation_size=validation_size,
                batch_size=batch_size,
                shuffle_set_split=shuffle_set_split,
                distributed_sampler=use_distributed_sampler(kwargs.get("strategy", None)),
                **datasplitter_kwargs,
            )
            ds.setup()
            train_dls[name] = ds.train_dataloader()
            test_dls[name] = ds.test_dataloader()
            val_dls[name] = ds.val_dataloader()
            self.train_indices_.append(ds.train_idx)
            self.test_indices_.append(ds.test_idx)
            self.validation_indices_.append(ds.val_idx)

        train_dl = TrainDL(train_dls)
        val_dl = TrainDL(val_dls)
        
        # Initialize and run training plan
        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else {}
        self._training_plan = self._training_plan_cls(
            self.module,
            **plan_kwargs,
        )
        if train_size == 1.0:
            self.trainer.fit(self._training_plan, train_dl)
        else:
            self.trainer.fit(
                self._training_plan, train_dataloaders=train_dl, val_dataloaders=val_dl
            )
        
        try:
            self.history_ = self.trainer.logger.history
        except AttributeError:
            self.history_ = None
        
        self.module.eval()
        self.to_device(device)
        self.is_trained_ = True

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        batch_key: str | None = None,
        labels_key: str | None = None,
        likelihood: Literal["nb", "zinb", "nbmixture", "normal", "log1pnormal", "ziln", "zig"] = "nb",
        normalize_lib: bool = True,
        gmm_prior: bool = False,
        semi_supervised: bool = False,
        n_mixture_components: int = 10,
        unlabeled_category: str = "unknown",
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        adata
            AnnData object to register.
        %(param_batch_key)s
        %(param_labels_key)s
        %(param_layer)s
        likelihood
            Likelihood model for this modality (default: 'nb').
            One of:
            - 'nb' : Negative Binomial
            - 'zinb' : Zero-Inflated Negative Binomial
            - 'nbmixture' : Negative Binomial Mixture (for protein data with background/foreground)
            - 'normal' : Normal distribution
            - 'log1pnormal' : Log1p Normal distribution
            - 'ziln' : Zero-Inflated Log Normal distribution
            - 'zig' : Zero-Inflated Gamma distribution
        normalize_lib
            Whether to normalize counts with library size in the model.
        gmm_prior
            Whether to use a GMM prior for this modality.
        semi_supervised
            Whether to use semi-supervised classification for this modality.
        n_mixture_components
            Number of mixture components for the GMM prior. If semi_supervised is True,
            this parameter is ignored and set to the number of unique labels in labels_key.
        unlabeled_category
            Category for unlabeled cells in labels_key.
        **kwargs
            Additional keyword arguments.
        """
        # Checks for valid likelihood and data compatibility
        continuous_likelihoods = {"ziln", "zig", "lognormal", "log1pnormal", "gamma", "normal"}
        if likelihood in continuous_likelihoods:
            # Check for negative values (invalid for most continuous likelihoods)
            data = adata.X if layer is None else adata.layers[layer]
            if scipy.sparse.issparse(data):
                data_min = data.min()
            else:
                data_min = np.min(data)
            
            if data_min < 0 and likelihood not in {"normal"}:
                raise ValueError(
                    f"Likelihood '{likelihood}' requires non-negative data, "
                    f"but found minimum value {data_min}. Consider using 'normal' instead."
                )
            # Warn about zeros for non-zero-inflated likelihoods
            if likelihood in {"lognormal", "gamma"}:
                if scipy.sparse.issparse(data):
                    has_zeros = (data.data == 0).any() or data.nnz < np.prod(data.shape)
                else:
                    has_zeros = (data == 0).any()
                
                if has_zeros:
                    warnings.warn(
                        f"Data contains zeros but likelihood '{likelihood}' does not model zeros. "
                        f"Consider using '{{'ziln' if likelihood == 'lognormal' else 'zig'}}' instead.",
                        UserWarning,
                        stacklevel=settings.warnings_stacklevel,
                    )

        if scipy.sparse.issparse(adata.X) and not isinstance(adata.X, scipy.sparse.csr_matrix):
            adata.X = adata.X.tocsr()
        if layer in adata.layers and not isinstance(adata.layers[layer], scipy.sparse.csr_matrix):
            adata.layers[layer] = adata.layers[layer].tocsr()
        
        adata.uns["diagvi_likelihood"] = likelihood
        adata.uns["diagvi_normalize_lib"] = normalize_lib
        adata.uns["diagvi_gmm_prior"] = gmm_prior
        adata.uns["diagvi_semi_supervised"] = semi_supervised

        # If semi-supervised, set n_mixture_components to number of unique labels
        if semi_supervised and labels_key is not None:
            n_mixture_components = adata.obs[labels_key].nunique()
        adata.uns["diagvi_n_mixture_components"] = n_mixture_components

        setup_method_args = cls._get_setup_method_args(**locals())
        if labels_key is None:
            label_field = CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key)
        else:
            label_field = LabelsWithUnlabeledObsField(
                REGISTRY_KEYS.LABELS_KEY, labels_key, unlabeled_category
            )
        
        is_count_data = likelihood in {"nb", "zinb", "nbmixture"}
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=is_count_data),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            label_field,
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    @classmethod
    def setup_mudata(
        cls,
        mdata: MuData,
        modalities: list[str],
        batch_key: dict[str, str] | str | None = None,
        labels_key: dict[str, str] | str | None = None,
        layer: dict[str, str | None] | str | None = None,
        likelihood: dict[str, Literal["nb", "zinb", "nbmixture", "normal"]] | str = "nb",
        normalize_lib: dict[str, bool] | bool = True,
        gmm_prior: dict[str, bool] | bool = False,
        semi_supervised: dict[str, bool] | bool = False,
        n_mixture_components: dict[str, int] | int = 10,
        unlabeled_category: dict[str, str] | str = "unknown",
        **kwargs,
    ):
        """Register a MuData object for use with DIAGVI.

        Parameters
        ----------
        mdata
            MuData object containing multiple modalities.
        modalities
            List of two modality names from mdata.mod to use.
        batch_key
            Key(s) in adata.obs for batch annotation. Can be a single string
            (applied to all modalities) or a dict mapping modality names to keys.
        labels_key
            Key(s) in adata.obs for cell type labels.
        layer
            Layer(s) in adata to use as input.
        likelihood
            Likelihood model(s) for each modality. Options: 'nb', 'zinb',
            'nbmixture', 'normal'. Default is 'nb'.
        normalize_lib
            Whether to normalize counts with library size in the model for each modality.
        gmm_prior
            Whether to use GMM prior for each modality. Default is False.
        semi_supervised
            Whether to use semi-supervised learning for each modality.
            Default is False.
        n_mixture_components
            Number of GMM mixture components for each modality. Default is 10.
        unlabeled_category
            Category name for unlabeled cells. Default is 'unknown'.
        **kwargs
            Additional keyword arguments passed to setup_anndata.
        """
        for mod_key in modalities:
            adata = mdata.mod[mod_key]

            batch_key_mod = batch_key[mod_key] if isinstance(batch_key, dict) else batch_key
            labels_key_mod = labels_key[mod_key] if isinstance(labels_key, dict) else labels_key
            layer_mod = layer[mod_key] if isinstance(layer, dict) else layer
            likelihood_mod = likelihood[mod_key] if isinstance(likelihood, dict) else likelihood
            normalize_lib_mod = normalize_lib[mod_key] if isinstance(normalize_lib, dict) else normalize_lib
            gmm_prior_mod = gmm_prior[mod_key] if isinstance(gmm_prior, dict) else gmm_prior
            semi_supervised_mod = (
                semi_supervised[mod_key] if isinstance(semi_supervised, dict) else semi_supervised
            )
            n_mixture_components_mod = (
                n_mixture_components[mod_key]
                if isinstance(n_mixture_components, dict)
                else n_mixture_components
            )
            unlabeled_category_mod = (
                unlabeled_category[mod_key]
                if isinstance(unlabeled_category, dict)
                else unlabeled_category
            )

            cls.setup_anndata(
                adata=adata,
                batch_key=batch_key_mod,
                labels_key=labels_key_mod,
                layer=layer_mod,
                likelihood=likelihood_mod,
                normalize_lib=normalize_lib_mod,
                gmm_prior=gmm_prior_mod,
                semi_supervised=semi_supervised_mod,
                n_mixture_components=n_mixture_components_mod,
                unlabeled_category=unlabeled_category_mod,
                **kwargs,
            )

    @staticmethod
    @dependencies("torch_geometric")
    def construct_custom_guidance_graph(
        input_dict: dict[str, AnnData],
        mapping_df: pd.DataFrame,
        weight: float = 1.0,
        sign: float = 1.0
    ) -> Data:
        """Constructs a custom guidance graph defined by the mapping_df.

        Parameters
        ----------
        input_dict
            Dictionary mapping modality names
            (e.g., "scRNAseq", "seqFISH") to their respective AnnData objects.
            The keys in this dictionary must match the column names in `mapping_df`.
        mapping_df
            DataFrame specifying feature correspondences between modalities.
            Each column should correspond to a modality name in input_dict,
            and each row defines a feature pair.
        weight
            Edge weight assigned to cross-modality edges.
        sign
            Edge sign assigned to cross-modality edges.

        Returns
        -------
        A PyTorch Geometric Data object representing the guidance graph,
        including node features, edge indices, edge weights, edge signs,
        including node features, edge indices, edge weights, edge signs,
        and modality-specific feature indices.

        Notes
        -----
        - Features are renamed to include the modality name as a suffix.
        - Self-loops are added for all features with weight=1.0 and sign=1.0.
        - Missing features in the mapping are silently skipped.
        """
        from torch_geometric.data import Data

        modality_names = list(input_dict.keys())
        adata1, adata2 = list(input_dict.values())

        ad_1_ft = [f"{f}_{modality_names[0]}" for f in adata1.var_names]
        ad_2_ft = [f"{f}_{modality_names[1]}" for f in adata2.var_names]

        all_features = ad_1_ft + ad_2_ft
        feature_to_index = {f: i for i, f in enumerate(all_features)}

        edge_index = []
        edge_weight = []
        edge_sign = []

        adata_adt_vars = []
        adata_gex_vars = []

        for ft_pair in range(mapping_df.shape[0]):
            pair = mapping_df.iloc[ft_pair, :]
            diss_ft = f"{pair[modality_names[0]]}_{modality_names[0]}"
            sp_ft = f"{pair[modality_names[1]]}_{modality_names[1]}"

            if diss_ft not in feature_to_index or sp_ft not in feature_to_index:
                continue

            adata_adt_vars.append(sp_ft)
            adata_gex_vars.append(diss_ft)

            i = feature_to_index[diss_ft]
            j = feature_to_index[sp_ft]

            edge_index += [[i, j], [j, i]]
            edge_weight += [weight, weight]
            edge_sign += [sign, sign]

        for feature in all_features:
            i = feature_to_index[feature]
            edge_index.append([i, i])
            edge_weight.append(1.0)
            edge_sign.append(1.0)

        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
        edge_weight = torch.tensor(edge_weight, dtype=torch.float)
        edge_sign = torch.tensor(edge_sign, dtype=torch.float)

        x = torch.eye(len(all_features))

        indices_1 = torch.tensor([feature_to_index[f] for f in ad_1_ft], dtype=torch.long)
        indices_2 = torch.tensor([feature_to_index[f] for f in ad_2_ft], dtype=torch.long)

        guidance_graph = Data(
            x=x,
            edge_index=edge_index,
            edge_weight=edge_weight,
            edge_sign=edge_sign,
            **{
                f"{modality_names[0]}_indices": indices_1,
                f"{modality_names[1]}_indices": indices_2,
            },
        )

        return guidance_graph

    def _make_scvi_dls(
        self, adatas: list[AnnData] = None, batch_size: int = 1024
    ) -> list[AnnDataLoader]:
        """Create AnnDataLoaders for each modality.

        Parameters
        ----------
        adatas
            List of AnnData objects. If None, uses self.adatas.
        batch_size
            Minibatch size for data loading.

        Returns
        -------
        List of data loaders, one per modality.
        """
        if adatas is None:
            adatas = self.adatas
        post_list = [self._make_data_loader(ad, batch_size=batch_size) for ad in adatas]
        for i, dl in enumerate(post_list):
            dl.mode = i
        return post_list

    @torch.inference_mode()
    def get_latent_representation(
        self,
        adatas: dict[str, AnnData] | list[AnnData] | None = None,
        batch_size: int = 1024,
    ) -> dict[str, np.ndarray]:
        """Return the latent space embedding for each dataset.

        Parameters
        ----------
        adatas
            Modalities to compute embeddings for. If None, uses all modalities.
        batch_size
            Minibatch size for data loading.

        Returns
        -------
        Dictionary mapping modality names to latent embeddings.
        """
        if adatas is None:
            adatas = self.adatas
        if isinstance(adatas, dict):
            input_names = list(adatas.keys())
            adata_list = list(adatas.values())
        else:
            input_names = self.input_names
            adata_list = adatas
        
        scdls = self._make_scvi_dls(adata_list, batch_size=batch_size)
        self.module.eval()

        latents = {}
        for input_name, scdl in zip(input_names, scdls, strict=False):
            latent = []
            for tensors in scdl:
                latent.append(
                    self.module.inference(
                        **self.module._get_inference_input(tensors), mode=input_name
                    )["z"]
                    .cpu()
                    .detach()
                )
            latent = torch.cat(latent).numpy()
            latents[input_name] = latent
        
        return latents

    @torch.inference_mode()
    def posterior_predictive_sample(
        self,
        adatas: dict[str, AnnData] | None = None,
        indices: dict[str, Sequence[int]] | None = None,
        n_samples: int = 1,
        batch_size: int | None = None,
    ) -> dict[str, np.ndarray]:
        """Generate posterior predictive samples for each modality.

        This method generates samples from the posterior predictive distribution
        p(x̂|x) for each modality, which can be used for model criticism via
        :class:`~scvi.criticism.PosteriorPredictiveCheck`.

        Parameters
        ----------
        adatas
            Dictionary mapping modality names to AnnData objects. If ``None``,
            defaults to the AnnData objects used to initialize the model.
        indices
            Dictionary mapping modality names to indices of cells to use.
            If ``None``, all cells are used for each modality.
        n_samples
            Number of posterior predictive samples to generate per cell.
        batch_size
            Minibatch size for data loading. Defaults to ``scvi.settings.batch_size``.

        Returns
        -------
        Dictionary mapping modality names to arrays of shape
        ``(n_cells, n_features, n_samples)`` if ``n_samples > 1``, else
        ``(n_cells, n_features)``.
        """
        if adatas is None:
            adatas = self.adatas

        if indices is None:
            indices = {name: None for name in adatas.keys()}

        # Validate that all provided modalities are valid
        if not set(adatas.keys()).issubset(set(self.input_names)):
            invalid_keys = set(adatas.keys()) - set(self.input_names)
            raise ValueError(
                f"Invalid modality names: {invalid_keys}. "
                f"Must be subset of model input names: {set(self.input_names)}"
            )

        self.module.eval()

        # Create data loaders for each modality
        data_loaders = {}
        for name, adata in adatas.items():
            mode_indices = indices.get(name)
            dl = self._make_data_loader(
                adata,
                indices=mode_indices,
                batch_size=batch_size or settings.batch_size,
            )
            data_loaders[name] = dl

        # Collect samples for each modality
        samples_per_modality = {name: [] for name in adatas.keys()}

        # Process each modality separately to handle different numbers of cells
        for mode in adatas.keys():
            mode_samples = []
            for tensors in data_loaders[mode]:
                # Create a dict with only this modality's tensors
                mode_tensors_dict = {mode: tensors}

                # Run module.sample which handles inference + generative + sampling
                mode_result = self.module.sample(mode_tensors_dict, n_samples=n_samples)

                mode_sample = mode_result[mode]

                if n_samples > 1:
                    # mode_sample shape: (n_samples, n_cells, n_features)
                    # Transpose to (n_cells, n_features, n_samples)
                    mode_sample = mode_sample.permute(1, 2, 0)

                mode_samples.append(mode_sample.numpy())

            # Concatenate along cell dimension (axis 0)
            samples_per_modality[mode] = np.concatenate(mode_samples, axis=0)

        return samples_per_modality

    @torch.inference_mode()
    def get_imputed_values(
        self,
        source_name: str,
        source_adata: AnnData | None = None,
        deterministic: bool = True,
        batch_size: int = 1024,
        target_batch: int | str | np.ndarray | None = None,
        target_libsize: float | np.ndarray | None = None,
    ) -> np.ndarray:
        """Return imputed values for a given source modality.

        Parameters
        ----------
        source_name
            Name of the source modality.
        source_adata
            AnnData object for the source modality. If None, uses the registered AnnData.
        batch_size
            Minibatch size for data loading.
        target_batch
            Target batch index or array for imputation.
        target_libsize
            Target library size(s) for imputation.

        Returns
        -------
        Imputed values.
        """
        # Choose source adata according to mode
        if source_name not in self.adatas:
            raise ValueError(f"`source_name` must be one of {list(self.adatas.keys())}!")
        if source_adata is None:
            source_adata = self.adatas[source_name]
        
        # Get batch categories for this modality
        batch_manager = self.adata_managers[source_name]
        batch_categories = batch_manager.get_state_registry(REGISTRY_KEYS.BATCH_KEY).get(
            "categorical_mapping", None
        )

        self.module.eval()
        dl = self._make_data_loader(source_adata, batch_size=batch_size)
        reconstructed_values = []
        for tensor in dl:
            inference_output = self.module.inference(
                **self.module._get_inference_input(tensor),
                mode=source_name,
                deterministic=deterministic,
            )
            # Use the feature embedding of the other modality for reconstruction
            inference_output["v"] = inference_output["v_other"]
            
            generative_input = self.module._get_generative_input(
                tensor,
                inference_output,
            )

            # Handle target batch
            batch_size_ = generative_input[MODULE_KEYS.LIBRARY_KEY].shape[0]
            if target_batch is not None:
                b = np.asarray(target_batch)
                if b.size == 1:
                    b = np.full(batch_size_, b.item())
                elif b.size != batch_size_:  # raise error if wrong size
                    raise ValueError("`target_batch` must have the same size as adata!")
                if batch_categories is not None and not np.issubdtype(b.dtype, np.integer):
                    b = np.array([np.where(batch_categories == lbl)[0][0] for lbl in b])
            # Use batch index zero if no target batch is provided
            else:
                b = np.zeros(batch_size_, dtype=int)
            b = np.asarray(b, dtype=int).reshape(-1, 1)
            generative_input[MODULE_KEYS.BATCH_INDEX_KEY] = torch.tensor(
                b,
                dtype=torch.long,
                device=generative_input[MODULE_KEYS.LIBRARY_KEY].device,
            )
            
            # Handle target libsize
            if target_libsize is not None:
                l = target_libsize
                if not isinstance(l, np.ndarray):
                    l = np.asarray(l)
                l = l.squeeze()
                if l.ndim == 0:
                    l = l[np.newaxis]
                elif l.ndim > 1:
                    raise ValueError("`target_libsize` cannot be >1 dimensional")
                if l.size == 1:
                    l = np.repeat(l, batch_size_)
                if l.size != batch_size_:
                    raise ValueError("`target_libsize` must have the same size as the batch!")
                l = l.reshape((-1, 1))
                generative_input[MODULE_KEYS.LIBRARY_KEY] = torch.tensor(
                    l,
                    dtype=generative_input[MODULE_KEYS.LIBRARY_KEY].dtype,
                    device=generative_input[MODULE_KEYS.LIBRARY_KEY].device,
                )
            
            # Determine target modality
            target_names = [m for m in self.input_names if m != source_name]
            if len(target_names) != 1:
                raise ValueError("There must be exactly two modalities defined.")
            target_name = target_names[0]
            
            generative_output = self.module.generative(**generative_input, mode=target_name)
            # Use distribution mean for correct expected value across all likelihoods
            px_dist = generative_output[MODULE_KEYS.PX_KEY]
            reconstructed_values.append(px_dist.mean.cpu().detach())
        
        reconstructed_value = torch.cat(reconstructed_values).numpy()
        return reconstructed_value

    def save(
        self,
        dir_path: str,
        prefix: str | None = None,
        overwrite: bool = False,
        save_anndata: bool = False,
        save_kwargs: dict | None = None,
    ):
        """Save the DIAGVI model and optionally the AnnData objects.

        Parameters
        ----------
        dir_path
            Directory path to save the model and AnnData objects.
        prefix
            Prefix for the saved model file name.
        overwrite
            Whether to overwrite the directory if it already exists.
        save_anndata
            Whether to save the AnnData objects used to train the model.
        save_kwargs
            Additional keyword arguments for torch.save.

        Raises
        ------
        ValueError
            If dir_path exists and overwrite is False.
        """
        if not os.path.exists(dir_path) or overwrite:
            os.makedirs(dir_path, exist_ok=overwrite)
        else:
            raise ValueError(
                f"{dir_path} already exists. Please provide an unexisting directory for saving."
            )

        file_name_prefix = prefix or ""
        save_kwargs = save_kwargs or {}

        adatas = self.adatas
        if save_anndata:
            for key in self.input_names:
                ad = adatas[key]
                save_path = os.path.join(dir_path, f"adata_{key}.h5ad")
                ad.write(save_path)

        model_state_dict = self.module.state_dict()

        var_names = {}
        for key in self.input_names:
            var_names[key] = adatas[key].var_names

        # Get all the user attributes
        user_attributes = self._get_user_attributes()

        # Only save the public attributes with _ at the very end
        user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}

        model_save_path = os.path.join(dir_path, f"{file_name_prefix}model.pt")

        torch.save(
            {
                "model_state_dict": model_state_dict,
                "names": self.input_names,
                **{f"var_names_{mod}": var_names[mod] for mod in self.input_names},
                "attr_dict": user_attributes,
            },
            model_save_path,
            **save_kwargs,
        )

    @classmethod
    @devices_dsp.dedent
    def load(
        cls,
        dir_path: str,
        accelerator: str = "auto",
        device: int | str = "auto",
        prefix: str | None = None,
        backup_url: str | None = None,
    ) -> DIAGVI:
        """Load a saved DIAGVI model from disk.

        Parameters
        ----------
        dir_path
            Directory path where the model was saved.
        %(param_accelerator)s
        %(param_device)s
        prefix
            Prefix used when saving the model.
        backup_url
            URL to download the model from if not found locally.

        Returns
        -------
        Loaded DIAGVI model instance.

        Raises
        ------
        ValueError
            If the saved model is from a different class or missing setup inputs.
        """
        _, _, device = parse_device_args(
            accelerator=accelerator,
            devices=device,
            return_device="torch",
        )

        (
            attr_dict,
            var_names_dict,
            model_state_dict,
            adatas,
        ) = _load_saved_diagvi_files(
            dir_path,
            prefix=prefix,
            map_location=device,
            backup_url=backup_url,
        )

        for mod in adatas.keys():
            saved_var_names = var_names_dict[mod]
            user_var_names = adatas[mod].var_names.astype(str)

            if not np.array_equal(saved_var_names, user_var_names):
                warnings.warn(
                    "var_names for adata passed in does not match var_names of adata "
                    "used to train the model. For valid results, the vars need to be the"
                    "same and in the same order as the adata used to train the model.",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )

        registries = attr_dict.pop("registries_")
        for adata, registry in zip(adatas, registries, strict=True):
            if _MODEL_NAME_KEY in registry and registry[_MODEL_NAME_KEY] != cls.__name__:
                raise ValueError("It appears you are loading a model from a different class.")

            if _SETUP_ARGS_KEY not in registry:
                raise ValueError(
                    "Saved model does not contain original setup inputs. "
                    "Cannot load the original setup."
                )

            cls.setup_anndata(adatas[adata], source_registry=registry, **registry[_SETUP_ARGS_KEY])

        # Get the parameters for the class init signature
        init_params = attr_dict.pop("init_params_")

        # New saving and loading, enable backwards compatibility
        if "non_kwargs" in init_params.keys():
            # Grab all the parameters except for kwargs (is a dict)
            non_kwargs = init_params["non_kwargs"]
            kwargs = init_params["kwargs"]
            kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        else:
            # Grab all the parameters except for kwargs (is a dict)
            non_kwargs = {k: v for k, v in init_params.items() if not isinstance(v, dict)}
            kwargs = {k: v for k, v in init_params.items() if isinstance(v, dict)}
            kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}

        non_kwargs.pop("adatas", None)
        kwargs.pop("adatas", None)

        model = cls(adatas, **non_kwargs, **kwargs)

        for attr, val in attr_dict.items():
            setattr(model, attr, val)

        model.module.load_state_dict(model_state_dict)
        model.module.eval()
        model.to_device(device)
        return model


class TrainDL(DataLoader):
    """DataLoader that creates batch structure for training process by combining multiple data loaders.

    The data loaders for each modality are zipped together, and smaller datasets are cycled
    to match the size of the largest dataset. Each yielded batch is a dictionary mapping
    modality names to their respective batches.

    Parameters
    ----------
    data_loader_dict
        Dictionary mapping modality names to their respective DataLoaders.
    """

    def __init__(self, data_loader_dict: dict[str, DataLoader]):
        self.data_loader_dict = data_loader_dict
        self.input_names = list(data_loader_dict.keys())
        self.data_loader_list = list(data_loader_dict.values())

        # Identify the largest dataset
        self.largest_train_dl_idx = np.argmax(
            [len(dl.indices) for dl in self.data_loader_list]
        )

        # DataLoader corresponding to the largest dataset
        self.largest_dl = self.data_loader_list[
            self.largest_train_dl_idx
        ]

        super().__init__(
            self.largest_dl,
            num_workers=settings.dl_num_workers,
            persistent_workers=getattr(settings, "dl_persistent_workers", False),
            pin_memory=getattr(settings, "dl_pin_memory_gpu_training", False),
        )

    # Number of batches per epoch is determined by the larger dataset
    def __len__(self):
        return len(self.largest_dl)

    # Cyclic iteration
    def __iter__(self):
        # Produces batches by zipping together
        # Ensures that smaller datasets are repeated indefinitely
        train_dls = [
            (
                dl if i == self.largest_train_dl_idx else cycle(dl)
            )  # Repeat smaller dataset indefinitely
            for i, dl in enumerate(self.data_loader_list)
        ]

        # Return a dict of modality_name -> batch for each step
        for batches in zip(*train_dls, strict=False):
            yield dict(zip(self.input_names, batches, strict=False))
