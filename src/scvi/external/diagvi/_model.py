import logging
import os
import warnings
from itertools import cycle
from typing import Literal

import numpy as np
import pandas as pd
import scipy.sparse
import torch
from anndata import AnnData
from sklearn.neighbors import NearestNeighbors
from torch.utils.data import DataLoader

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data._constants import _MODEL_NAME_KEY, _SETUP_ARGS_KEY
from scvi.data.fields import CategoricalObsField, LabelsWithUnlabeledObsField, LayerField
from scvi.dataloaders import AnnDataLoader, DataSplitter
from scvi.external.diagvi._utils import (
    _check_guidance_graph_consisteny,
    _construct_guidance_graph,
    _load_saved_diagvi_files,
)
from scvi.model._utils import get_max_epochs_heuristic, parse_device_args, use_distributed_sampler
from scvi.model.base import BaseModelClass, VAEMixin
from scvi.module._constants import MODULE_KEYS
from scvi.train import Trainer
from scvi.utils import dependencies

from ._module import DIAGVAE
from ._task import DiagTrainingPlan

logger = logging.getLogger(__name__)


class DIAGVI(BaseModelClass, VAEMixin):
    """
    Diagonal Multi-Modal Integration Variational Inference (DIAGVI) model.

    Integrates multi-modal single-cell data using a guidance graph and supports
    semi-supervised learning and GMM priors.

    Parameters
    ----------
    adatas : dict[str, AnnData]
        Dictionary mapping input names to AnnData objects.
    guidance_graph : Data, optional
        Precomputed guidance graph. If None, it will be constructed from the data
        by using overlapping feature names.
    mapping_df : pd.DataFrame, optional
        DataFrame specifying feature correspondences between modalities
        (used to compute the guidance graph).
    **model_kwargs : dict
        Additional keyword arguments for the DIAGVAE module.

    Examples
    --------
    >>> adatas = {"rna_data": adata_rna, "protein_data": adata_protein}
    >>> model = DIAGVI(adatas)
    >>> model.train()
    """

    # from torch_geometric.data import Data
    def __init__(
        self,
        adatas: dict[str, AnnData],
        # guidance_graph: Data | None = None,
        guidance_graph=None,
        mapping_df: pd.DataFrame | None = None,
        common_scale: bool | bool = True,
        **model_kwargs: dict,
    ) -> None:
        super().__init__()
        self.adatas = adatas
        self.input_names = list(adatas.keys())
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
        modalities = {name: adata.uns["diagvi_modality"] for name, adata in self.adatas.items()}
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
        _check_guidance_graph_consisteny(self.guidance_graph, adatas)
        self.module = DIAGVAE(
            n_inputs=n_inputs,
            n_batches=n_batches,
            n_labels=n_labels,
            modality_likelihoods=generative_distributions,
            modalities=modalities,
            guidance_graph=self.guidance_graph,
            use_gmm_prior=gmm_priors,
            semi_supervised=semi_supervised,
            n_mixture_components=n_mixture_components,
            common_scale=common_scale,
            **model_kwargs,
        )
        self._model_summary_string = (
            f"DiagVI Model with the following params: input names: {self.input_names}, "
            f"n_inputs: {n_inputs}, n_batches: {n_batches}, "
            f"generative distributions: {generative_distributions}"
        )
        self.init_params_ = self._get_init_params(locals())
        logger.info(self._model_summary_string)

    def train(
        self,
        max_epochs: int | None = None,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float = 0.9,
        shuffle_set_split: bool = True,
        batch_size: int = 1024,
        # datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,  # kwargs passed to trainingplan
        **kwargs,
    ) -> None:
        """
        Train the DIAGVI model.

        Parameters
        ----------
        max_epochs : int, optional
            Maximum number of training epochs. If None, a heuristic is used.
        accelerator : str, optional
            Accelerator to use for training (e.g., 'cpu', 'gpu', 'auto').
        devices : int, list[int], or str, optional
            Devices to use for training.
        train_size : float, optional
            Proportion of data to use for training (rest for validation).
        shuffle_set_split : bool, optional
            Whether to shuffle data before splitting into train/validation.
        batch_size : int, optional
            Minibatch size for training.
        plan_kwargs : dict, optional
            Additional keyword arguments for the training plan.
        **kwargs
            Additional keyword arguments for the Trainer.
        """
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
        self.trainer = Trainer(
            max_epochs=max_epochs,
            enable_progress_bar=True,
            progress_bar_refresh_rate=1,
            early_stopping=True,
            early_stopping_monitor="validation_loss",
            early_stopping_min_delta=1,
            early_stopping_patience=10,
            early_stopping_mode="min",
            accelerator=accelerator,
            devices=devices,
            **kwargs,
        )
        validation_size = 1 - train_size
        self.train_indices_, self.test_indices_, self.validation_indices_ = [], [], []
        train_dls, test_dls, val_dls = {}, {}, {}
        for name, adm in self.adata_managers.items():
            ds = DataSplitter(
                adm,
                train_size=train_size,
                validation_size=validation_size,
                batch_size=batch_size,
                shuffle_set_split=shuffle_set_split,
                distributed_sampler=use_distributed_sampler(kwargs.get("strategy", None)),
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
        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else {}
        self._training_plan = DiagTrainingPlan(
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
    def setup_anndata(
        cls,
        adata: AnnData,
        batch_key: str | None = None,
        labels_key: str | None = None,
        layer: str | None = None,
        likelihood: Literal[
            "nb", "zinb", "nbmixture", "normal", "log1pnormal", "ziln", "zig"
        ] = "nb",
        modality: Literal["rna", "protein"] = "rna",
        gmm_prior: bool = False,
        semi_supervised: bool = False,
        n_mixture_components: int = 10,
        unlabeled_category: str = "unknown",
        **kwargs: dict,
    ) -> None:
        """
        Register an AnnData object for use with DIAGVI.

        Parameters
        ----------
        adata : AnnData
            AnnData object to register.
        batch_key : str, optional
            Key in adata.obs for batch annotation.
        labels_key : str, optional
            Key in adata.obs for cell type labels.
        layer : str, optional
            Layer in adata to use as input.
        likelihood : {'nb', 'zinb', 'nbmixture', 'normal', 'log1pnormal', 'ziln', 'zig'}, optional
            Likelihood model for this modality (default: 'nb').
        modality : {'rna', 'protein'}, optional
            Name of the modality (default: 'rna').
        gmm_prior : bool, optional
            Whether to use a GMM prior for this modality.
        semi_supervised : bool, optional
            Whether to use semi-supervised classification for this modality.
        n_mixture_components : int, optional
            Number of mixture components for the GMM prior.
        unlabeled_category : str, optional
            Category for unlabeled cells in labels_key.
        **kwargs : dict
            Additional keyword arguments.
        """
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
        adata.uns["diagvi_modality"] = modality
        adata.uns["diagvi_gmm_prior"] = gmm_prior
        adata.uns["diagvi_semi_supervised"] = semi_supervised
        adata.uns["diagvi_n_mixture_components"] = n_mixture_components
        # Set up the anndata object for the model
        setup_method_args = cls._get_setup_method_args(
            **locals()
        )  # returns dict organizing the args used to call setup anndata
        if labels_key is None:
            label_field = CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key)
        else:
            label_field = LabelsWithUnlabeledObsField(
                REGISTRY_KEYS.LABELS_KEY, labels_key, unlabeled_category
            )
        # TODO: think about how to handle `is_count_data` below
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            label_field,
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata)
        cls.register_manager(adata_manager)

    @classmethod
    @dependencies("torch_geometric")
    def construct_custom_guidance_graph(self, input_dict, mapping_df, weight=1.0, sign=1.0):
        """
        Constructs a custom guidance graph defined by the mapping_df

        Parameters
        ----------
        input_dict : dict
            Dictionary mapping modality names
            (e.g., "scRNAseq", "seqFISH") to their respective AnnData objects.
            The keys in this dictionary must match the column names in `mapping_df`.
        mapping_df : pandas.DataFrame
            DataFrame specifying feature correspondences between modalities.
            Each column should correspond to a modality name in `input_dict`,
            and each row defines a feature pair.
        weight : float, optional
            Edge weight assigned to all edges in the graph (default: 1.0).
        sign : float, optional
            Edge sign assigned to all edges in the graph (default: 1.0).

        Returns
        -------
        guidance_graph : torch_geometric.data.Data
            A PyTorch Geometric Data object representing the guidance graph,
            including node features,
            edge indices, edge weights, edge signs,
            and modality-specific feature indices.

        Notes
        -----
        - The function renames features in each AnnData object
        to include the modality name as a suffix.
        - Self-loops are added for all features.
        - The mapping DataFrame must have column names that exactly match
        the keys in `input_dict`.
        """
        from torch_geometric.data import Data

        modality_names = list(input_dict.keys())
        adata1, adata2 = list(input_dict.values())

        ad_1_ft = [f"{f}_{modality_names[0]}" for f in adata1.var_names]
        ad_2_ft = [f"{f}_{modality_names[1]}" for f in adata2.var_names]

        adata1.var_names = ad_1_ft
        adata2.var_names = ad_2_ft

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

        # Add self-loops
        for feature in all_features:
            i = feature_to_index[feature]
            edge_index.append([i, i])
            edge_weight.append(weight)
            edge_sign.append(sign)

        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
        edge_weight = torch.tensor(edge_weight, dtype=torch.float)
        edge_sign = torch.tensor(edge_sign, dtype=torch.float)

        x = torch.eye(len(all_features))  # node features as identity for simplicity

        # Extract seq/spa indices
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
        if adatas is None:
            adatas = self.adatas
        post_list = [self._make_data_loader(ad, batch_size=batch_size) for ad in adatas]
        for i, dl in enumerate(post_list):
            dl.mode = i
        return post_list

    @torch.inference_mode()
    def get_latent_representation(
        self,
        adatas: dict[str, AnnData] | list[AnnData] = None,
        batch_size: int = 1024,
    ) -> dict[np.ndarray]:
        """
        Return the latent space embedding for each dataset.

        Parameters
        ----------
        adatas : dict[str, AnnData] or list[AnnData], optional
            Modalities to compute embeddings for. If None, uses all modalities.
        batch_size : int, optional
            Minibatch size for data loading.

        Returns
        -------
        dict[np.ndarray]
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

    def compute_per_feature_confidence(self, feature_embedding, conf_method):
        knn = NearestNeighbors(n_neighbors=11, metric="cosine")
        knn.fit(feature_embedding)
        distances, indices = knn.kneighbors(feature_embedding)
        # Remove self (first neighbor is always the point itself)
        distances = distances[:, 1:]
        indices = indices[:, 1:]
        score = []
        # compute statistics for every node
        if conf_method == "min":
            score = distances.min(axis=1)
        elif conf_method == "mean":
            score = distances.mean(axis=1)
        elif conf_method == "max":
            score = distances.max(axis=1)
        elif conf_method == "median":
            score = np.median(distances, axis=1)
        return score

    @torch.inference_mode()
    def get_imputed_values(
        self,
        source_name: int,
        source_adata: AnnData | None = None,
        batch_size: int = 1024,
        target_batch: int | None = None,
        target_libsize: float | None = None,
        conf_method: Literal["min", "mean", "max", "median"] = "min",
        min_max_scale: bool = True,
    ) -> list[np.ndarray]:
        """
        Return imputed values and feature confidence scores for a given source modality.

        Parameters
        ----------
        source_name : str
            Name of the source modality.
        source_adata : AnnData, optional
            AnnData object for the source modality. If None, uses the registered AnnData.
        batch_size : int, optional
            Minibatch size for data loading.
        target_batch : int or array-like, optional
            Target batch index or array for imputation.
        target_libsize : float or array-like, optional
            Target library size(s) for imputation.
        conf_method : {'min', 'mean', 'max', 'median'}, optional
            Method for aggregating kNN distances for confidence scoring.
        min_max_scale : bool, optional
            Whether to min-max scale the confidence scores.

        Returns
        -------
        tuple (np.ndarray, np.ndarray)
            Imputed values and per-feature confidence scores.
        """
        # choose source adata according to mode
        if source_name not in self.adatas:
            raise ValueError(f"`source_name` must be one of {list(self.adatas.keys())}!")
        # Choose source AnnData
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
                **self.module._get_inference_input(tensor), mode=source_name
            )
            # use the feature embedding of the other modality for reconstruction
            inference_output["v"] = inference_output["v_other"]
            generative_input = self.module._get_generative_input(
                tensor,
                inference_output,
            )
            # handle target batch
            batch_size_ = generative_input[MODULE_KEYS.LIBRARY_KEY].shape[0]
            if target_batch is not None:
                b = np.asarray(target_batch)
                if b.size == 1:  # if size 1: broadcast to array of length Batch size
                    b = np.full(batch_size_, b.item())
                elif b.size != batch_size_:  # raise error if wrong size
                    raise ValueError("`target_batch` must have the same size as adata!")
                # Map categorical names to indices if needed
                if batch_categories is not None and not np.issubdtype(b.dtype, np.integer):
                    b = np.array([np.where(batch_categories == lbl)[0][0] for lbl in b])
            # use batch index zero if no target batch is provided
            else:
                b = np.zeros(batch_size_, dtype=int)
            b = np.asarray(b, dtype=int).reshape(-1, 1)
            generative_input[MODULE_KEYS.BATCH_INDEX_KEY] = torch.tensor(
                b,
                dtype=torch.long,
                device=generative_input[MODULE_KEYS.LIBRARY_KEY].device,
            )
            # handle target libsuze
            if target_libsize is not None:
                l = target_libsize
                if not isinstance(l, np.ndarray):
                    l = np.asarray(l)
                l = l.squeeze()
                if l.ndim == 0:  # Scalar
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
            # determine target modality
            target_names = [m for m in self.input_names if m != source_name]
            if len(target_names) != 1:
                raise ValueError("There must be exactly two modalities defined.")
            target_name = target_names[0]
            generative_output = self.module.generative(**generative_input, mode=target_name)
            # use distribution mean for correct expected value across all likelihoods
            # TODO: check if this is correct for all likelihoods
            px_dist = generative_output[MODULE_KEYS.PX_KEY]
            reconstructed_values.append(px_dist.mean.cpu().detach())
            # extract the final feature embedding
            feature_embedding = inference_output["v"].cpu().detach().numpy()
        score = self.compute_per_feature_confidence(feature_embedding, conf_method)
        if min_max_scale:
            if len(score) > 0:
                score_norm = (score - score.min()) / (score.max() - score.min())
            else:
                score_norm = score.copy()
        reconstructed_value = torch.cat(reconstructed_values).numpy()
        return reconstructed_value, score_norm

    def save(
        self,
        dir_path: str,
        prefix: str | None = None,
        overwrite: bool = False,
        save_anndata: bool = False,
        save_kwargs: dict | None = None,
        **anndata_write_kwargs,
    ):
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

        # get all the user attributes
        user_attributes = self._get_user_attributes()

        # only save the public attributes with _ at the very end
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
    # @devices_dsp.dedent
    def load(
        cls,
        dir_path: str,
        # adata_seq: AnnData | None = None,
        # adata_spatial: AnnData | None = None,
        accelerator: str = "auto",
        device: int | str = "auto",
        prefix: str | None = None,
        backup_url: str | None = None,
    ):
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

        # get the parameters for the class init signature
        init_params = attr_dict.pop("init_params_")

        # new saving and loading, enable backwards compatibility
        if "non_kwargs" in init_params.keys():
            # grab all the parameters except for kwargs (is a dict)
            non_kwargs = init_params["non_kwargs"]
            kwargs = init_params["kwargs"]
            kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        else:
            # grab all the parameters except for kwargs (is a dict)
            non_kwargs = {k: v for k, v in init_params.items() if not isinstance(v, dict)}
            kwargs = {k: v for k, v in init_params.items() if isinstance(v, dict)}
            kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}

        # Remove 'adatas' from non_kwargs and kwargs if present
        non_kwargs.pop("adatas", None)
        kwargs.pop("adatas", None)

        model = cls(adatas, **non_kwargs, **kwargs)

        for attr, val in attr_dict.items():
            setattr(model, attr, val)

        model.module.load_state_dict(model_state_dict)
        model.module.eval()
        model.to_device(device)
        return model


class TrainDL(DataLoader):  # creates batch structure for training process
    def __init__(self, data_loader_dict):
        # def __init__(self, data_loader_list, num_workers=0):
        self.data_loader_dict = data_loader_dict
        self.input_names = list(data_loader_dict.keys())
        self.data_loader_list = list(data_loader_dict.values())

        self.largest_train_dl_idx = np.argmax(
            [len(dl.indices) for dl in self.data_loader_list]
        )  # index of dl with largest num of samples
        self.largest_dl = self.data_loader_list[
            self.largest_train_dl_idx
        ]  # dl corresponding to the largest dataset
        super().__init__(
            self.largest_dl,
            # shuffle=True,
            num_workers=settings.dl_num_workers,
            persistent_workers=getattr(settings, "dl_persistent_workers", False),
        )

    # number of batches per epoch is determined by the larger dataset
    def __len__(self):
        return len(self.largest_dl)

    # cyclic iteration
    def __iter__(self):
        # produces batches by zipping together
        # ensures that smaller datasets are repeated indefinitely
        train_dls = [
            (
                dl if i == self.largest_train_dl_idx else cycle(dl)
            )  # repeat smaller dataset indefinitely
            for i, dl in enumerate(self.data_loader_list)
        ]
        # Return a dict of modality_name -> batch for each step
        for batches in zip(*train_dls, strict=False):
            yield dict(zip(self.input_names, batches, strict=False))
