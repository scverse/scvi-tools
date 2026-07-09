"""Utility functions for the DIAGVI model."""

from __future__ import annotations

import logging
import os
from collections.abc import Mapping
from itertools import cycle
from typing import TYPE_CHECKING

import anndata as ad
import numpy as np
import scipy.spatial
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader

from scvi.data._download import _download
from scvi.utils import dependencies

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any, Literal

    import pandas as pd
    from anndata import AnnData
    from torch_geometric.data import Data

logger = logging.getLogger(__name__)


class CyclicMultiDataLoader(DataLoader):
    """Combine multiple data loaders by cycling shorter loaders to match the longest one."""

    def __init__(
        self,
        data_loaders: Mapping[str, DataLoader] | Sequence[DataLoader],
        **data_loader_kwargs: Any,
    ):
        self._returns_mapping = isinstance(data_loaders, Mapping)
        if self._returns_mapping:
            self.input_names = list(data_loaders.keys())
            self.data_loader_list = list(data_loaders.values())
        else:
            self.input_names = None
            self.data_loader_list = list(data_loaders)

        if not self.data_loader_list:
            raise ValueError("At least one data loader is required.")

        self.largest_train_dl_idx = max(
            range(len(self.data_loader_list)),
            key=lambda idx: len(self.data_loader_list[idx].indices),
        )
        self.largest_dl = self.data_loader_list[self.largest_train_dl_idx]
        super().__init__(self.largest_dl, **data_loader_kwargs)

    def __len__(self):
        return len(self.largest_dl)

    def __iter__(self):
        data_loaders = [
            dl if i == self.largest_train_dl_idx else cycle(dl)
            for i, dl in enumerate(self.data_loader_list)
        ]

        if self._returns_mapping:
            for batches in zip(*data_loaders, strict=False):
                yield dict(zip(self.input_names, batches, strict=False))
        else:
            yield from zip(*data_loaders, strict=False)


@dependencies("torch_geometric")
def _construct_guidance_graph(
    adatas: dict[str, AnnData],
    mapping_df: pd.DataFrame | None,
    weight: float = 1.0,
    sign: float = 1.0,
) -> Data:
    """Construct a guidance graph linking features across modalities.

    Creates a bipartite graph where nodes represent features from each modality
    and edges connect corresponding features based on the mapping DataFrame or
    shared feature names.

    Parameters
    ----------
    adatas
        Dictionary mapping modality names to AnnData objects.
    mapping_df
        DataFrame with columns matching modality names, containing feature
        mappings. If None, uses shared feature names.
    weight
        Edge weight for cross-modality connections.
    sign
        Edge sign for cross-modality connections.

    Returns
    -------
    PyTorch Geometric Data object with node features, edge indices,
    edge weights, edge signs, and modality index tensors.

    Raises
    ------
    ValueError
        If not exactly two modalities are provided or no overlapping features
        exist when mapping_df is None.
    """
    from torch_geometric.data import Data

    if len(adatas) != 2:
        raise ValueError("Exactly two modalities are required.")
    input_names = list(adatas.keys())
    adata1, adata2 = adatas[input_names[0]], adatas[input_names[1]]

    if mapping_df is not None:
        features1 = list(adata1.var_names)
        features2 = list(adata2.var_names)
    else:
        shared_features = set(adata1.var_names) & set(adata2.var_names)
        if not shared_features:
            raise ValueError("No overlapping features between the two modalities.")

        features1 = [f"{f}_{input_names[0]}" for f in adata1.var_names]
        features2 = [f"{f}_{input_names[1]}" for f in adata2.var_names]

    all_features = features1 + features2
    feature_to_index = {f: i for i, f in enumerate(all_features)}

    edge_index = []
    edge_weight = []
    edge_sign = []

    if mapping_df is not None:
        for ft_pair in range(mapping_df.shape[0]):
            pair = mapping_df.iloc[ft_pair, :]
            diss_ft = pair[input_names[0]]
            sp_ft = pair[input_names[1]]

            i = feature_to_index[diss_ft]
            j = feature_to_index[sp_ft]

            edge_index += [[i, j], [j, i]]
            edge_weight += [weight, weight]
            edge_sign += [sign, sign]

    else:
        for feature in shared_features:
            i = feature_to_index[f"{feature}_{input_names[0]}"]
            j = feature_to_index[f"{feature}_{input_names[1]}"]

            edge_index += [[i, j], [j, i]]
            edge_weight += [weight, weight]
            edge_sign += [sign, sign]

    for feature in all_features:
        i = feature_to_index[feature]
        edge_index.append([i, i])
        edge_weight.append(weight)
        edge_sign.append(sign)

    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    edge_weight = torch.tensor(edge_weight, dtype=torch.float)
    edge_sign = torch.tensor(edge_sign, dtype=torch.float)

    x = torch.eye(len(all_features))

    indices1 = torch.tensor([feature_to_index[f] for f in features1], dtype=torch.long)
    indices2 = torch.tensor([feature_to_index[f] for f in features2], dtype=torch.long)

    return Data(
        x=x,
        edge_index=edge_index,
        edge_weight=edge_weight,
        edge_sign=edge_sign,
        **{f"{input_names[0]}_indices": indices1, f"{input_names[1]}_indices": indices2},
    )


def _check_guidance_graph_consistency(graph: Data, adatas: dict[str, AnnData]):
    """Validate guidance graph structure and consistency with AnnData objects.

    Performs several consistency checks on the guidance graph:
    1. Node count matches total number of features across modalities
    2. Required edge attributes (edge_weight, edge_sign) are present
    3. Self-loops exist for all nodes
    4. Graph is symmetric (undirected)

    Parameters
    ----------
    graph
        PyTorch Geometric Data object representing the guidance graph.
    adatas
        Dictionary mapping modality names to AnnData objects.

    Raises
    ------
    ValueError
        If any consistency check fails.
    """
    n_expected = sum(adata.shape[1] for adata in adatas.values())

    # 1. Check variable coverage via counts
    if graph.num_nodes != n_expected:
        raise ValueError(
            f"Graph node count {graph.num_nodes} does not match expected {n_expected}."
        )

    # 2. Check edge attributes
    for attr in ["edge_weight", "edge_sign"]:
        if not hasattr(graph, attr):
            raise ValueError(f"Graph missing required edge attribute: {attr}")
        if getattr(graph, attr).shape[0] != graph.edge_index.shape[1]:
            raise ValueError(f"Edge attribute {attr} does not match number of edges.")

    # 3. Check self-loops
    src, tgt = graph.edge_index
    self_loops = src == tgt
    n_self_loops = self_loops.sum().item()
    if n_self_loops < graph.num_nodes:
        raise ValueError("Graph is missing self-loops for some nodes.")

    # 4. Check symmetry (for undirected graphs: For every edge (i, j), check that (j, i) exists)
    edge_set = {(i.item(), j.item()) for i, j in zip(src, tgt, strict=False)}
    for i, j in zip(src, tgt, strict=False):
        if (j.item(), i.item()) not in edge_set:
            raise ValueError(
                f"Graph is not symmetric: edge ({i.item()}, {j.item()}) has no counterpart."
            )

    # If all checks pass
    logger.info("Guidance graph consistency checks passed.")


def _load_saved_diagvi_files(
    dir_path: str,
    prefix: str | None = None,
    map_location: Literal["cpu", "cuda"] | None = None,
    backup_url: str | None = None,
) -> tuple[dict, dict[str, np.ndarray], dict, dict[str, AnnData | None]]:
    """Loads saved DiagVI model and AnnData files from a directory.

    Parameters
    ----------
    dir_path
        Directory path where the model and AnnData files are stored.
    prefix
        Optional prefix for the file names.
    map_location
        Device mapping for loading the model.
    backup_url
        Optional URL to download the model file if not found locally.

    Returns
    -------
    A tuple containing:
    - attr_dict: Dictionary of model attributes.
    - var_names: Dictionary of variable names for each modality.
    - model_state_dict: State dictionary of the model.
    - adatas: Dictionary of AnnData objects for each modality.

    Raises
    ------
    ValueError
        If model file cannot be loaded.
    """
    file_name_prefix = prefix or ""

    model_file_name = f"{file_name_prefix}model.pt"
    model_path = os.path.join(dir_path, model_file_name)

    try:
        _download(backup_url, dir_path, model_file_name)
        model = torch.load(model_path, map_location=map_location, weights_only=False)
    except FileNotFoundError as exc:
        raise ValueError(f"Failed to load model file at {model_path}. ") from exc

    names = model["names"]

    adatas = {}
    var_names = {}
    for name in names:
        adata_path = os.path.join(dir_path, f"{file_name_prefix}adata_{name}.h5ad")
        if os.path.exists(adata_path):
            adatas[name] = ad.read_h5ad(adata_path)
            var_names[name] = adatas[name].var_names
        else:
            adatas[name] = None
            key = f"var_names_{name}"
            if key in model:
                var_names[name] = model[key]

    model_state_dict = model["model_state_dict"]
    attr_dict = model["attr_dict"]

    return (
        attr_dict,
        var_names,
        model_state_dict,
        adatas,
    )


@dependencies("torch_geometric")
def compute_graph_loss(graph: Data, feature_embeddings: torch.Tensor) -> torch.Tensor:
    """Compute graph reconstruction loss using negative sampling.

    Uses structured negative sampling to compute a contrastive loss that
    encourages connected nodes to have similar embeddings and unconnected
    nodes to have dissimilar embeddings.

    Parameters
    ----------
    graph
        PyTorch Geometric Data object with edge_index.
    feature_embeddings
        Tensor of shape (n_features, embedding_dim) containing feature embeddings.

    Returns
    -------
    Scalar tensor containing the graph reconstruction loss.
    """
    import torch_geometric

    edge_index = graph.edge_index
    edge_index_neg = torch_geometric.utils.structured_negative_sampling(edge_index)

    pos_i, pos_j, neg_j = edge_index_neg[0], edge_index_neg[1], edge_index_neg[2]

    vi = feature_embeddings[pos_i]
    vj = feature_embeddings[pos_j]
    vj_neg = feature_embeddings[neg_j]

    pos_logits = (vi * vj).sum(dim=1)
    pos_loss = F.logsigmoid(pos_logits).mean()

    neg_logits = (vi * vj_neg).sum(dim=1)
    neg_loss = F.logsigmoid(-neg_logits).mean()

    total_loss = -(pos_loss + neg_loss) / 2

    return total_loss


def kl_divergence_graph(mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
    """Computes the KL divergence for graph latent variables.

    Parameters
    ----------
    mu
        Mean tensor of the latent variables.
    logvar
        Log-variance tensor of the latent variables.

    Returns
    -------
    The mean KL divergence as a tensor.
    """
    kl = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1)
    kl_mean = kl.mean()
    return kl_mean


def compute_foscttm(
    latents: dict[str, np.ndarray],
    indices: list[np.ndarray] | None = None,
    downsample: bool = False,
    n_obs: int = 10000,
) -> dict[str, float]:
    """Compute FOSCTTM (Fraction of Samples Closer Than True Match) metric.

    For paired multi-modal data, measures how well the model aligns corresponding
    cells in the latent space. Lower values are better (0 = perfect alignment).
    Assumes cells are aligned: mod1 cell i corresponds to mod2 cell i.

    Parameters
    ----------
    latents
        Dictionary mapping modality names to latent arrays.
    indices
        List of two arrays: indices for [mod1, mod2].
        Usually obtained from model.validation_indices after training.
    downsample
        Whether to downsample when the number of observations is large.
    n_obs
        Number of observations to keep if downsample is True.
    """
    # Extract latent representations for the two modalities
    latent_mod1 = latents[list(latents.keys())[0]]
    latent_mod2 = latents[list(latents.keys())[1]]

    # Subset if indices provided
    if indices is not None:
        latent_mod1 = latent_mod1[indices[0]]
        latent_mod2 = latent_mod2[indices[1]]

    # Validate shapes match
    if latent_mod1.shape[0] != latent_mod2.shape[0]:
        raise ValueError("Shapes do not match!")

    n_cells = latent_mod1.shape[0]

    # Downsample if requested and number of cells exceeds n_obs
    if n_cells > n_obs and downsample:
        np.random.seed(0)
        sample_indices = np.random.choice(n_cells, size=n_obs, replace=False)
        latent_mod1 = latent_mod1[sample_indices]
        latent_mod2 = latent_mod2[sample_indices]
        n_cells = n_obs

    # Compute pairwise distances
    distances = scipy.spatial.distance_matrix(latent_mod1, latent_mod2)

    # Compute FOSCTTM
    foscttm_mod1_to_mod2 = (distances < np.expand_dims(np.diag(distances), axis=1)).mean(axis=1)
    foscttm_mod2_to_mod1 = (distances < np.expand_dims(np.diag(distances), axis=0)).mean(axis=0)

    # Aggregate metrics
    mean_mod1_to_mod2 = foscttm_mod1_to_mod2.mean()
    mean_mod2_to_mod1 = foscttm_mod2_to_mod1.mean()

    mod1_name, mod2_name = list(latents.keys())
    foscttm_metrics = {
        f"foscttm/{mod1_name}_to_{mod2_name}": float(mean_mod1_to_mod2),
        f"foscttm/{mod2_name}_to_{mod1_name}": float(mean_mod2_to_mod1),
        "foscttm/mean": float((mean_mod1_to_mod2 + mean_mod2_to_mod1) / 2),
    }

    return foscttm_metrics
