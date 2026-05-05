"""Graph-aware dataloaders for spatial single-cell models."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import torch
from torch.utils.data import default_convert

from scvi import REGISTRY_KEYS
from scvi.data import AnnTorchDataset
from scvi.dataloaders._ann_dataloader import AnnDataLoader
from scvi.dataloaders._data_splitting import DataSplitter

if TYPE_CHECKING:
    from scvi.data import AnnDataManager


def _as_torch_tensor(array: np.ndarray | torch.Tensor) -> torch.Tensor:
    """Convert AnnTorchDataset output to a torch tensor without changing sparse layout."""
    if isinstance(array, np.ndarray):
        return torch.from_numpy(array)
    return array


class _GraphBatchConverter:
    """Convert an AnnTorchDataset batch into a PyG Data object."""

    def __init__(
        self,
        full_adata_manager: AnnDataManager,
        neighbor_indices_key: str,
        edge_obsm_keys: list[str],
        load_sparse_neighbor_tensor: bool,
        load_neighbor_expression: bool,
    ):
        self.neighbor_indices_key = neighbor_indices_key
        self.edge_obsm_keys = edge_obsm_keys
        self.load_neighbor_expression = load_neighbor_expression
        if load_neighbor_expression:
            self._full_dataset = AnnTorchDataset(
                full_adata_manager,
                getitem_tensors=[REGISTRY_KEYS.X_KEY],
                load_sparse_tensor=load_sparse_neighbor_tensor,
            )

    def __call__(self, batch: dict[str, np.ndarray | torch.Tensor]):
        try:
            from torch_geometric.data import Data
        except ImportError as error:
            raise ImportError(
                "torch_geometric is required for GraphDataLoader. "
                "Install it with: pip install torch_geometric"
            ) from error

        batch = default_convert(batch)
        ind_neighbors = batch[self.neighbor_indices_key].long()
        n_obs, n_neighbors = ind_neighbors.shape

        x = _as_torch_tensor(batch[REGISTRY_KEYS.X_KEY])

        center_idx = torch.arange(n_obs, dtype=torch.long).repeat_interleave(n_neighbors)
        neighbor_idx = torch.arange(n_obs * n_neighbors, dtype=torch.long)
        edge_index = torch.stack([center_idx, neighbor_idx], dim=0)

        edge_attrs = []
        for key in self.edge_obsm_keys:
            vals = batch[key].float()
            edge_attrs.append(vals.reshape(n_obs * n_neighbors, -1))
        edge_attr = torch.cat(edge_attrs, dim=1) if edge_attrs else None

        data_kwargs = dict(batch)
        data_kwargs.update(
            {
                "x": x,
                "edge_index": edge_index,
                "edge_attr": edge_attr,
                "distances_n": batch.get("distance_neighbor"),
            }
        )
        if self.load_neighbor_expression:
            flat_neighbors = ind_neighbors.cpu().numpy().ravel()
            data_kwargs["x_n"] = _as_torch_tensor(
                self._full_dataset[flat_neighbors][REGISTRY_KEYS.X_KEY]
            )
        return Data(**data_kwargs)


class GraphDataLoader(AnnDataLoader):
    """DataLoader that yields mini-batches as :class:`torch_geometric.data.Data` objects.

    Each batch contains center cells and their pre-fetched spatial neighbors. Neighbor
    expression is looked up from ``full_adata_manager`` so neighbors outside the current
    train/validation/test split are intentionally allowed, matching existing RESOLVI behavior.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` for the split being loaded.
    full_adata_manager
        :class:`~scvi.data.AnnDataManager` for all observations. Used for neighbor expression
        lookup, including cross-split neighbors.
    indices
        Observation indices to load from ``adata_manager``.
    neighbor_indices_key
        Registry key containing neighbor indices, shape ``[N, K]``.
    edge_obsm_keys
        Registry keys to flatten and concatenate into ``edge_attr``. Each key must have shape
        ``[N, K]`` or ``[N, K, D]``. Defaults to ``["distance_neighbor"]``.
    load_sparse_neighbor_tensor
        If ``True``, loads sparse neighbor expression as sparse torch tensors. This avoids
        densifying neighbor expression on the CPU before device transfer.
    load_neighbor_expression
        If ``False``, omits ``x_n`` and leaves neighbor expression gathering to the model. This is
        useful when a model keeps a device-resident expression cache.
    **kwargs
        Forwarded to :class:`~scvi.dataloaders.AnnDataLoader`.
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        full_adata_manager: AnnDataManager,
        indices: list[int] | list[bool] | None = None,
        neighbor_indices_key: str = "index_neighbor",
        edge_obsm_keys: list[str] | None = None,
        load_sparse_neighbor_tensor: bool = True,
        load_neighbor_expression: bool = True,
        **kwargs,
    ):
        if "collate_fn" in kwargs:
            raise ValueError("GraphDataLoader uses its own collate_fn to build graph batches.")
        if kwargs.pop("iter_ndarray", False):
            raise ValueError("GraphDataLoader does not support `iter_ndarray=True`.")

        super().__init__(adata_manager, indices=indices, **kwargs)
        self.neighbor_indices_key = neighbor_indices_key
        self.edge_obsm_keys = (
            list(edge_obsm_keys) if edge_obsm_keys is not None else ["distance_neighbor"]
        )
        self.load_sparse_neighbor_tensor = load_sparse_neighbor_tensor
        self.load_neighbor_expression = load_neighbor_expression
        self._graph_batch_converter = _GraphBatchConverter(
            full_adata_manager,
            neighbor_indices_key=self.neighbor_indices_key,
            edge_obsm_keys=self.edge_obsm_keys,
            load_sparse_neighbor_tensor=load_sparse_neighbor_tensor,
            load_neighbor_expression=load_neighbor_expression,
        )
        self.collate_fn = self._graph_batch_converter


class GraphDataSplitter(DataSplitter):
    """DataSplitter that creates :class:`GraphDataLoader` instances.

    Parameters
    ----------
    load_sparse_neighbor_tensor
        Forwarded to :class:`GraphDataLoader`.
    load_neighbor_expression
        Forwarded to :class:`GraphDataLoader`.
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        neighbor_indices_key: str = "index_neighbor",
        edge_obsm_keys: list[str] | None = None,
        load_sparse_neighbor_tensor: bool = True,
        load_neighbor_expression: bool = True,
        **kwargs,
    ):
        super().__init__(adata_manager, **kwargs)
        self.neighbor_indices_key = neighbor_indices_key
        self.edge_obsm_keys = (
            list(edge_obsm_keys) if edge_obsm_keys is not None else ["distance_neighbor"]
        )
        self.load_sparse_neighbor_tensor = load_sparse_neighbor_tensor
        self.load_neighbor_expression = load_neighbor_expression

    def _make_graph_dataloader(
        self,
        indices: np.ndarray,
        shuffle: bool,
        drop_last: bool,
    ) -> GraphDataLoader:
        return GraphDataLoader(
            self.adata_manager,
            full_adata_manager=self.adata_manager,
            indices=indices,
            shuffle=shuffle,
            drop_last=drop_last,
            load_sparse_tensor=self.load_sparse_tensor,
            pin_memory=self.pin_memory,
            neighbor_indices_key=self.neighbor_indices_key,
            edge_obsm_keys=self.edge_obsm_keys,
            load_sparse_neighbor_tensor=self.load_sparse_neighbor_tensor,
            load_neighbor_expression=self.load_neighbor_expression,
            **self.data_loader_kwargs,
        )

    def train_dataloader(self) -> GraphDataLoader:
        """Create graph train dataloader."""
        return self._make_graph_dataloader(
            self.train_idx,
            shuffle=True,
            drop_last=self.drop_last,
        )

    def val_dataloader(self) -> GraphDataLoader | None:
        """Create graph validation dataloader."""
        if len(self.val_idx) > 0:
            return self._make_graph_dataloader(
                self.val_idx,
                shuffle=False,
                drop_last=False,
            )

    def test_dataloader(self) -> GraphDataLoader | None:
        """Create graph test dataloader."""
        if len(self.test_idx) > 0:
            return self._make_graph_dataloader(
                self.test_idx,
                shuffle=False,
                drop_last=False,
            )
