"""RESOLVI-specific LightningDataModule for annbatch streaming."""

import numpy as np
import torch

from scvi import REGISTRY_KEYS
from scvi.dataloaders._custom_dataloaders import AnnbatchDataModule


class _TensorIndexDataset:
    """Wraps a (n_obs, n_vars) float tensor for neighbor expression lookup.

    Supports ``ds[array_of_indices, :]["X"]`` as expected by the RESOLVI module.
    """

    def __init__(self, x: torch.Tensor):
        self._x = x

    def __getitem__(self, idx):
        if isinstance(idx, tuple):
            idx = idx[0]
        if isinstance(idx, np.ndarray):
            idx = torch.from_numpy(idx.astype(np.int64))
        return {"X": self._x[idx]}


class RESOLVISparseDataModule(AnnbatchDataModule):
    """AnnbatchDataModule subclass for streaming RESOLVI training.

    Extends :class:`~scvi.dataloaders.AnnbatchDataModule` with spatial
    neighbour tensors and pre-computed dataset-level priors that RESOLVI
    needs when ``adata`` is not provided.

    Parameters
    ----------
    dataset
        ``annbatch.Loader`` configured with the prepared AnnData.
    x_all
        Full expression matrix as a float tensor — used for neighbour
        expression lookup during forward passes.
    n_neigh
        Number of spatial neighbours per cell.
    precomputed_priors
        Dict with keys ``background_ratio``, ``median_distance``,
        ``mean_log_counts``, ``std_log_counts``.
    batch_key, label_key, unlabeled_category, model_name,
    categorical_covariate_keys, var_names
        Forwarded to :class:`~scvi.dataloaders.AnnbatchDataModule`.
    """

    def __init__(
        self,
        dataset,
        x_all: torch.Tensor,
        n_neigh: int,
        precomputed_priors: dict,
        batch_key=None,
        label_key=None,
        unlabeled_category="Unknown",
        model_name="RESOLVI",
        categorical_covariate_keys=None,
        var_names=None,
    ):
        super().__init__(
            dataset,
            batch_key=batch_key,
            label_key=label_key,
            unlabeled_category=unlabeled_category,
            model_name=model_name,
            categorical_covariate_keys=categorical_covariate_keys,
            var_names=var_names,
        )
        self._n_neigh = n_neigh
        self.precomputed_priors = precomputed_priors
        self.expression_anntorchdata = _TensorIndexDataset(x_all)

    def on_before_batch_transfer(self, batch, dataloader_idx):
        """Extend base transfer with RESOLVI spatial-neighbour tensors."""
        result = super().on_before_batch_transfer(batch, dataloader_idx)
        obs = batch.get("obs")
        n_cells = result["X"].shape[0]

        # Global cell index used by Pyro plate subsampling
        if obs is not None and "_global_index" in obs.columns:
            ind_x = torch.from_numpy(obs["_global_index"].values.astype(np.int64)).unsqueeze(1)
        else:
            ind_x = torch.arange(n_cells, dtype=torch.long).unsqueeze(1)
        result[REGISTRY_KEYS.INDICES_KEY] = ind_x

        # Spatial neighbour indices/distances stored as flat obs columns
        if obs is not None and "_index_neighbor_0" in obs.columns:
            idx_cols = [
                obs[f"_index_neighbor_{k}"].values.astype(np.int64) for k in range(self._n_neigh)
            ]
            dist_cols = [
                obs[f"_distance_neighbor_{k}"].values.astype(np.float32)
                for k in range(self._n_neigh)
            ]
            result["index_neighbor"] = torch.from_numpy(np.stack(idx_cols, axis=1))
            result["distance_neighbor"] = torch.from_numpy(np.stack(dist_cols, axis=1))
        else:
            result["index_neighbor"] = torch.zeros((n_cells, self._n_neigh), dtype=torch.long)
            result["distance_neighbor"] = torch.zeros(
                (n_cells, self._n_neigh), dtype=torch.float32
            )

        # PyroJitGuideWarmup calls .to(device) on all values; remove None entries
        return {k: v for k, v in result.items() if v is not None}

    @property
    def registry(self) -> dict:
        base = super().registry
        fr = base["field_registries"]

        # n_distance_neighbor drives RESOLVI's n_neighbors hyperparameter
        fr["index_neighbor"] = {
            "data_registry": {"attr_name": "obsm", "attr_key": "index_neighbor"},
            "state_registry": {},
            "summary_stats": {"n_distance_neighbor": self._n_neigh},
        }
        fr["distance_neighbor"] = {
            "data_registry": {"attr_name": "obsm", "attr_key": "distance_neighbor"},
            "state_registry": {},
            "summary_stats": {},
        }
        # n_ind_x: total cell count for Pyro plate
        fr["ind_x"] = {
            "data_registry": {"attr_name": "obs", "attr_key": "_global_index"},
            "state_registry": {},
            "summary_stats": {"n_ind_x": self.n_obs},
        }
        return base
