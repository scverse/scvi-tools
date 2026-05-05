"""Unit tests for GraphDataLoader."""

import numpy as np
import pytest
import torch

from scvi.data import synthetic_iid
from scvi.external import RESOLVI


@pytest.fixture(scope="module")
def resolvi_adata():
    adata = synthetic_iid(generate_coordinates=True, n_regions=5)
    adata.obsm["X_spatial"] = adata.obsm["coordinates"]
    RESOLVI.setup_anndata(adata)
    return adata


@pytest.fixture(scope="module")
def adata_manager(resolvi_adata):
    return RESOLVI._get_most_recent_anndata_manager(resolvi_adata, required=True)


def test_graph_dataloader_yields_data_objects(adata_manager):
    """Each batch must be a torch_geometric Data object, not a plain dict."""
    from torch_geometric.data import Data

    from scvi.dataloaders import GraphDataLoader

    dl = GraphDataLoader(
        adata_manager,
        full_adata_manager=adata_manager,
        batch_size=32,
        shuffle=False,
    )
    batch = next(iter(dl))
    assert isinstance(batch, Data), f"Expected Data, got {type(batch)}"


def test_graph_dataloader_shapes(adata_manager):
    """x: [N, G], x_n: [N*K, G], edge_index: [2, N*K], edge_attr: [N*K, 1]."""
    from scvi.dataloaders import GraphDataLoader

    batch_size = 32
    dl = GraphDataLoader(
        adata_manager,
        full_adata_manager=adata_manager,
        batch_size=batch_size,
        shuffle=False,
    )
    batch = next(iter(dl))

    n_genes = adata_manager.adata.n_vars
    K = adata_manager.adata.obsm["index_neighbor"].shape[1]
    N = batch.x.shape[0]  # may be < batch_size on last batch

    assert batch.x.shape == (N, n_genes), f"x shape wrong: {batch.x.shape}"
    assert batch.x_n.shape == (N * K, n_genes), f"x_n shape wrong: {batch.x_n.shape}"
    assert batch.edge_index.shape == (2, N * K), (
        f"edge_index shape wrong: {batch.edge_index.shape}"
    )
    assert batch.edge_attr.shape[0] == N * K, f"edge_attr rows wrong: {batch.edge_attr.shape}"


def test_graph_dataloader_edge_index_correctness(adata_manager):
    """edge_index[0] = center indices 0..N-1 repeated K times."""
    from scvi.dataloaders import GraphDataLoader

    dl = GraphDataLoader(
        adata_manager,
        full_adata_manager=adata_manager,
        batch_size=16,
        shuffle=False,
    )
    batch = next(iter(dl))
    N = batch.x.shape[0]
    K = adata_manager.adata.obsm["index_neighbor"].shape[1]

    expected_src = torch.arange(N).repeat_interleave(K)
    torch.testing.assert_close(batch.edge_index[0], expected_src)

    expected_dst = torch.arange(N * K)
    torch.testing.assert_close(batch.edge_index[1], expected_dst)


def test_graph_dataloader_edge_attr_from_keys(adata_manager):
    """distance_neighbor in edge_obsm_keys appears in edge_attr as [N*K, 1]."""
    from scvi.dataloaders import GraphDataLoader

    dl = GraphDataLoader(
        adata_manager,
        full_adata_manager=adata_manager,
        batch_size=16,
        shuffle=False,
        edge_obsm_keys=["distance_neighbor"],
    )
    batch = next(iter(dl))
    N = batch.x.shape[0]
    K = adata_manager.adata.obsm["index_neighbor"].shape[1]

    # distance_neighbor is [N, K] flattened to [N*K, 1]
    assert batch.edge_attr.shape == (N * K, 1)
    assert batch.edge_attr.dtype == torch.float32


def test_graph_dataloader_cross_split_neighbors(adata_manager):
    """Neighbor indices outside train split must resolve without error."""
    from scvi.dataloaders import GraphDataLoader

    n_obs = adata_manager.adata.n_obs
    train_indices = np.arange(n_obs // 2)  # first half only

    dl = GraphDataLoader(
        adata_manager,
        full_adata_manager=adata_manager,  # full dataset — cross-split intentional
        indices=train_indices,
        batch_size=16,
        shuffle=False,
    )
    for batch in dl:
        assert batch.x_n is not None
        break


def test_graph_dataloader_missing_torch_geometric(adata_manager, monkeypatch):
    """ImportError with install hint when torch_geometric is absent."""
    import builtins
    import sys

    from scvi.dataloaders import GraphDataLoader

    real_import = builtins.__import__

    def mock_import(name, *args, **kwargs):
        if name == "torch_geometric" or name.startswith("torch_geometric."):
            raise ImportError("No module named 'torch_geometric'")
        return real_import(name, *args, **kwargs)

    dl = GraphDataLoader(
        adata_manager,
        full_adata_manager=adata_manager,
        batch_size=4,
        shuffle=False,
    )

    # Remove torch_geometric from the module cache so our __import__ mock takes effect.
    # Without this, Python returns the cached module and never calls __import__.
    tg_keys = [
        k for k in sys.modules if k == "torch_geometric" or k.startswith("torch_geometric.")
    ]
    for k in tg_keys:
        monkeypatch.delitem(sys.modules, k)

    monkeypatch.setattr(builtins, "__import__", mock_import)
    with pytest.raises(ImportError, match="torch_geometric"):
        next(iter(dl))


def test_graph_datasplitter_returns_graph_dataloaders(adata_manager):
    """train_dataloader and val_dataloader must return GraphDataLoader instances."""
    from scvi.dataloaders import GraphDataLoader, GraphDataSplitter

    splitter = GraphDataSplitter(
        adata_manager,
        batch_size=32,
        train_size=0.8,
        validation_size=0.1,
    )
    splitter.setup()

    assert isinstance(splitter.train_dataloader(), GraphDataLoader)
    assert isinstance(splitter.val_dataloader(), GraphDataLoader)


def test_graph_datasplitter_train_batches_are_data(adata_manager):
    """Iterating train dataloader from GraphDataSplitter yields Data objects."""
    from torch_geometric.data import Data

    from scvi.dataloaders import GraphDataSplitter

    splitter = GraphDataSplitter(adata_manager, batch_size=32, train_size=0.8)
    splitter.setup()
    batch = next(iter(splitter.train_dataloader()))

    assert isinstance(batch, Data)
    assert hasattr(batch, "x")
    assert hasattr(batch, "x_n")
    assert hasattr(batch, "edge_index")


def test_graph_datasplitter_custom_edge_keys(adata_manager):
    """edge_obsm_keys forwarded to the underlying GraphDataLoader."""
    from scvi.dataloaders import GraphDataSplitter

    splitter = GraphDataSplitter(
        adata_manager,
        batch_size=16,
        train_size=0.8,
        edge_obsm_keys=["distance_neighbor"],
    )
    splitter.setup()

    assert splitter.train_dataloader().edge_obsm_keys == ["distance_neighbor"]
