import torch
from torch_geometric.data import Data


def _construct_guidance_graph(adata_seq, adata_spatial, weight=1.0, sign=1):
    shared_features = set(adata_seq.var_names) & set(adata_spatial.var_names)
    if not shared_features:
        raise ValueError("No overlapping features between the two modalities.")

    features_seq = [f"{f}_seq" for f in adata_seq.var_names]
    features_spatial = [f"{f}_spatial" for f in adata_spatial.var_names]

    # Build node list
    all_features = features_seq + features_spatial
    feature_to_index = {f: i for i, f in enumerate(all_features)}

    # Edges for matching features
    edge_index = []
    edge_weight = []
    edge_sign = []

    for feature in shared_features:
        i = feature_to_index[feature + "_seq"]
        j = feature_to_index[feature + "_spatial"]

        edge_index += [[i, j], [j, i]]
        edge_weight += [weight, weight]
        edge_sign += [sign, sign]

    # Add self-loops
    for feature in all_features:
        i = feature_to_index[feature]
        edge_index.append([i, i])
        edge_weight.append(1.0)
        edge_sign.append(1)

    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    edge_weight = torch.tensor(edge_weight, dtype=torch.float)
    edge_sign = torch.tensor(edge_sign, dtype=torch.float)

    x = torch.eye(len(all_features))  # node features as identity for simplicity

    # Extract seq/spa indices
    seq_indices = torch.tensor([feature_to_index[f] for f in features_seq], dtype=torch.long)
    spa_indices = torch.tensor([feature_to_index[f] for f in features_spatial], dtype=torch.long)

    return Data(
        x=x,
        edge_index=edge_index,
        edge_weight=edge_weight,
        edge_sign=edge_sign,
        seq_indices=seq_indices,  # attach as attributes
        spa_indices=spa_indices,
    )
