import logging
import os
from typing import Literal

import anndata as ad
import numpy as np
import torch
import torch.nn.functional as F
from anndata import AnnData
from torch_geometric.data import Data
from torch_geometric.utils import structured_negative_sampling

from scvi.data._download import _download

logger = logging.getLogger(__name__)


def _construct_guidance_graph(adatas, mapping_df, weight=1.0, sign=1):
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

    # Build node list
    all_features = features1 + features2
    feature_to_index = {f: i for i, f in enumerate(all_features)}

    # Edges for matching features
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
    indices1 = torch.tensor([feature_to_index[f] for f in features1], dtype=torch.long)
    indices2 = torch.tensor([feature_to_index[f] for f in features2], dtype=torch.long)

    return Data(
        x=x,
        edge_index=edge_index,
        edge_weight=edge_weight,
        edge_sign=edge_sign,
        **{f"{input_names[0]}_indices": indices1, f"{input_names[1]}_indices": indices2},
    )


def _check_guidance_graph_consisteny(graph: Data, adatas: dict[AnnData]):
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
) -> tuple[dict, dict, np.ndarray, np.ndarray, dict, AnnData | None, AnnData | None]:
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

    model_state_dict = model["model_state_dict"]
    attr_dict = model["attr_dict"]

    return (
        attr_dict,
        var_names,
        model_state_dict,
        adatas,
    )


def compute_graph_loss(graph, feature_embeddings):
    edge_index = graph.edge_index
    edge_index_neg = structured_negative_sampling(edge_index)

    pos_i = edge_index_neg[0].cpu().numpy()
    pos_j = edge_index_neg[1].cpu().numpy()
    neg_j = edge_index_neg[2].cpu().numpy()

    vi = feature_embeddings[pos_i]
    vj = feature_embeddings[pos_j]
    vj_neg = feature_embeddings[neg_j]

    pos_logits = (vi * vj).sum(dim=1)
    pos_loss = F.logsigmoid(pos_logits).mean()

    neg_logits = (vi * vj_neg).sum(dim=1)
    neg_loss = F.logsigmoid(-neg_logits).mean()

    total_loss = -(pos_loss + neg_loss) / 2

    return total_loss


def compute_sinkhorn_lam(lam_sinkhorn, epoch_current, epoch_sinkhorn):
    lam_sinkhorn_curr = lam_sinkhorn
    if epoch_sinkhorn:
        if epoch_current < epoch_sinkhorn:
            lam_sinkhorn_curr = (epoch_current / epoch_sinkhorn) * lam_sinkhorn
    return lam_sinkhorn_curr


def kl_divergence_graph(mu, logvar):
    kl = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1)
    kl_mean = kl.mean()
    return kl_mean
