import os
from itertools import product

import anndata
import mudata
import numpy as np
import pytest
import torch

from scvi import REGISTRY_KEYS
from scvi.data import synthetic_iid
from scvi.dataloaders import AnnTorchDataset
from tests.dataset.utils import (
    generic_setup_adata_manager,
    generic_setup_mudata_manager,
)


@pytest.mark.parametrize(
    "batch_size, n_batches, n_genes, sparse",
    product([128], [1, 2], [50], [True, False]),
)
def test_basic_anntorchdataset(
    batch_size: int, n_batches: int, n_genes: int, sparse: bool
):
    adata = synthetic_iid(
        batch_size=batch_size, n_batches=n_batches, n_genes=n_genes, sparse=sparse
    )
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    # test init
    dataset = AnnTorchDataset(manager)
    assert hasattr(dataset, "data")  # data should be loaded on init
    assert len(dataset) == n_batches * batch_size
    attr_and_types = {k: np.float32 for k in manager.data_registry}
    assert dataset.attributes_and_types == attr_and_types

    # getitem_tensors default
    ## single index
    data = dataset[0]
    X, batch, labels = (
        data[REGISTRY_KEYS.X_KEY],
        data[REGISTRY_KEYS.BATCH_KEY],
        data[REGISTRY_KEYS.LABELS_KEY],
    )
    assert isinstance(X, np.ndarray)
    assert isinstance(batch, np.ndarray)
    assert isinstance(labels, np.ndarray)
    assert X.dtype == np.float32
    assert batch.dtype == np.float32
    assert labels.dtype == np.float32
    assert X.shape == (1, n_genes)  # single observation should be a row vector
    assert batch.shape == (1, 1)
    assert labels.shape == (1, 1)
    ## slice
    data = dataset[:10]
    X, batch, labels = (
        data[REGISTRY_KEYS.X_KEY],
        data[REGISTRY_KEYS.BATCH_KEY],
        data[REGISTRY_KEYS.LABELS_KEY],
    )
    assert X.shape == (10, n_genes)
    assert batch.shape == (10, 1)
    assert labels.shape == (10, 1)
    ## list of indices
    data = dataset[[1, 0, 99]]
    X, batch, labels = (
        data[REGISTRY_KEYS.X_KEY],
        data[REGISTRY_KEYS.BATCH_KEY],
        data[REGISTRY_KEYS.LABELS_KEY],
    )
    assert X.shape == (3, n_genes)
    assert batch.shape == (3, 1)
    assert labels.shape == (3, 1)

    # getitem_tensors with list
    dataset = AnnTorchDataset(
        manager, getitem_tensors=[REGISTRY_KEYS.X_KEY, REGISTRY_KEYS.BATCH_KEY]
    )
    assert dataset.attributes_and_types == {
        REGISTRY_KEYS.X_KEY: np.float32,
        REGISTRY_KEYS.BATCH_KEY: np.float32,
    }
    ## single index
    data = dataset[0]
    X, batch = data[REGISTRY_KEYS.X_KEY], data[REGISTRY_KEYS.BATCH_KEY]
    assert X.dtype == np.float32
    assert batch.dtype == np.float32
    assert X.shape == (1, n_genes)  # single observation should be a row vector
    assert batch.shape == (1, 1)
    ## slice
    data = dataset[:10]
    X, batch = data[REGISTRY_KEYS.X_KEY], data[REGISTRY_KEYS.BATCH_KEY]
    assert X.shape == (10, n_genes)
    assert batch.shape == (10, 1)
    ## list of indices
    data = dataset[[1, 0, 99]]
    X, batch = data[REGISTRY_KEYS.X_KEY], data[REGISTRY_KEYS.BATCH_KEY]
    assert X.shape == (3, n_genes)
    assert batch.shape == (3, 1)

    # getitem_tensors with dict
    getitem_tensors = {
        REGISTRY_KEYS.X_KEY: np.float64,
        REGISTRY_KEYS.BATCH_KEY: np.int32,
    }
    dataset = AnnTorchDataset(manager, getitem_tensors=getitem_tensors)
    assert dataset.attributes_and_types == getitem_tensors
    ## single index
    data = dataset[0]
    X, batch = data[REGISTRY_KEYS.X_KEY], data[REGISTRY_KEYS.BATCH_KEY]
    assert X.dtype == np.float64
    assert batch.dtype == np.int32
    assert X.shape == (1, n_genes)  # single observation should be a row vector
    assert batch.shape == (1, 1)
    ## slice
    data = dataset[:10]
    X, batch = data[REGISTRY_KEYS.X_KEY], data[REGISTRY_KEYS.BATCH_KEY]
    assert X.shape == (10, n_genes)
    assert batch.shape == (10, 1)
    ## list of indices
    data = dataset[[1, 0, 99]]
    X, batch = data[REGISTRY_KEYS.X_KEY], data[REGISTRY_KEYS.BATCH_KEY]
    assert X.shape == (3, n_genes)
    assert batch.shape == (3, 1)


@pytest.mark.parametrize(
    "batch_size, n_batches, n_genes, sparse",
    product([128], [1], [50], [True, False]),
)
def test_errors_anntorchdataset(
    batch_size: int, n_batches: int, n_genes: int, sparse: bool
):
    adata = synthetic_iid(
        batch_size=batch_size, n_batches=n_batches, n_genes=n_genes, sparse=sparse
    )
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    # basic error handling
    dataset = AnnTorchDataset(manager)
    with pytest.raises(IndexError):
        _ = dataset[1000]

    with pytest.raises(ValueError):
        _ = AnnTorchDataset(manager, getitem_tensors=["not_a_key"])

    with pytest.raises(ValueError):
        _ = AnnTorchDataset(manager, getitem_tensors={"not_a_key": np.float32})

    with pytest.raises(ValueError):
        _ = AnnTorchDataset(
            manager, getitem_tensors={REGISTRY_KEYS.X_KEY: "not_a_dtype"}
        )


@pytest.mark.parametrize(
    "batch_size, n_batches, n_genes, sparse",
    product([128], [1], [50], [True, False]),
)
def test_disk_backed_anntorchdataset(
    save_path: str, batch_size: int, n_batches: int, n_genes: int, sparse: bool
):
    adata_path = os.path.join(
        save_path, f"adata_{batch_size}_{n_batches}_{sparse}.h5ad"
    )
    adata = synthetic_iid(
        batch_size=batch_size, n_batches=n_batches, n_genes=n_genes, sparse=sparse
    )
    adata.write(adata_path)
    del adata

    adata = anndata.read_h5ad(adata_path, backed="r")
    assert adata.isbacked
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    # test init
    dataset = AnnTorchDataset(manager)
    assert hasattr(dataset, "data")  # data should be loaded on init
    assert len(dataset) == n_batches * batch_size

    # getitem_tensors default
    ## single index
    data = dataset[0]
    X, batch, labels = (
        data[REGISTRY_KEYS.X_KEY],
        data[REGISTRY_KEYS.BATCH_KEY],
        data[REGISTRY_KEYS.LABELS_KEY],
    )
    assert X.dtype == np.float32
    assert batch.dtype == np.float32
    assert labels.dtype == np.float32
    assert X.shape == (1, n_genes)  # single observation should be a row vector
    assert batch.shape == (1, 1)
    assert labels.shape == (1, 1)


@pytest.mark.parametrize(
    "batch_size, n_batches, n_genes, sparse",
    product([128], [1], [50], [True, False]),
)
def test_cuda_backed_anntorchdataset(
    cuda: bool, batch_size: int, n_batches: int, n_genes: int, sparse: bool
):
    assert cuda

    adata = synthetic_iid(
        batch_size=batch_size, n_batches=n_batches, n_genes=n_genes, sparse=sparse
    )
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    # test init
    dataset = AnnTorchDataset(manager, accelerator="cuda", device_backed=True)
    assert hasattr(dataset, "data")  # data should be loaded on init
    assert len(dataset) == n_batches * batch_size
    data = dataset[0]
    X, batch, labels = (
        data[REGISTRY_KEYS.X_KEY],
        data[REGISTRY_KEYS.BATCH_KEY],
        data[REGISTRY_KEYS.LABELS_KEY],
    )
    assert isinstance(X, torch.Tensor)
    assert isinstance(batch, torch.Tensor)
    assert isinstance(labels, torch.Tensor)
    assert X.device.type == "cuda"
    assert batch.device.type == "cuda"
    assert labels.device.type == "cuda"
    assert X.dtype == torch.float32
    assert batch.dtype == torch.float32
    assert labels.dtype == torch.float32


@pytest.mark.parametrize(
    "batch_size, n_batches, n_genes, n_proteins, sparse",
    product([128], [1], [50], [20], [True, False]),
)
def test_protein_expression_obsm_anntorchdataset(
    batch_size: int,
    n_batches: int,
    n_genes: int,
    n_proteins: int,
    sparse: bool,
):
    adata = synthetic_iid(
        batch_size=batch_size,
        n_batches=n_batches,
        n_genes=n_genes,
        n_proteins=n_proteins,
        sparse=sparse,
    )
    manager = generic_setup_adata_manager(
        adata, batch_key="batch", protein_expression_obsm_key="protein_expression"
    )
    dataset = AnnTorchDataset(manager)
    for value in dataset[0].values():
        assert isinstance(value, np.ndarray)


@pytest.mark.parametrize(
    "batch_size, n_batches, n_genes, n_proteins, sparse",
    product([128], [1], [50], [20], [True, False]),
)
def test_mudata_anntorchdataset(
    batch_size: int,
    n_batches: int,
    n_genes: int,
    n_proteins: int,
    sparse: bool,
):
    adata = synthetic_iid(
        batch_size=batch_size,
        n_batches=n_batches,
        n_genes=n_genes,
        sparse=sparse,
    )
    bdata = synthetic_iid(
        batch_size=batch_size,
        n_batches=n_batches,
        n_genes=n_proteins,
        sparse=sparse,
    )
    mdata = mudata.MuData({"rna": adata, "protein": bdata})
    manager = generic_setup_mudata_manager(
        mdata,
        layer_mod="rna",
        layer=None,
        batch_mod="rna",
        batch_key="batch",
        protein_expression_mod="protein",
        protein_expression_layer=None,
    )

    # test init
    dataset = AnnTorchDataset(manager)
    assert hasattr(dataset, "data")  # data should be loaded on init
    assert len(dataset) == n_batches * batch_size
    attr_and_types = {k: np.float32 for k in manager.data_registry}
    assert dataset.attributes_and_types == attr_and_types
