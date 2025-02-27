import numpy as np
import pandas as pd
import pytest

from scvi.data import synthetic_iid
from scvi.external import CellAssign


def get_test_adata_marker_mat():
    adata = synthetic_iid(
        n_labels=5,
    )
    adata.obs["size_factor"] = adata.X.sum(1)

    marker_df = pd.DataFrame(data=np.random.randint(2, size=(100, 5)))
    marker_df.index = pd.Index([f"gene_{i}" for i in range(100)])

    return adata, marker_df


def test_cellassign():
    adata, marker_df = get_test_adata_marker_mat()
    CellAssign.setup_anndata(
        adata,
        "size_factor",
        batch_key="batch",
    )
    model = CellAssign(adata, marker_df)
    model.train(max_epochs=1)
    model.predict()
    model.get_normalized_expression()
    model.get_normalized_expression(transform_batch="batch_1")
    # model.get_normalized_expression(n_samples=2)


def test_cellassign_error_duplicates():
    adata, marker_df = get_test_adata_marker_mat()
    # Add a duplicate for internal test
    marker_df = marker_df._append(marker_df.iloc[1], ignore_index=False)
    CellAssign.setup_anndata(
        adata,
        "size_factor",
        batch_key="batch",
    )
    with pytest.raises(AssertionError):
        CellAssign(adata, marker_df)


def test_cellassign_covariates():
    adata, marker_df = get_test_adata_marker_mat()
    adata.obs["test"] = np.ones((adata.n_obs,))
    CellAssign.setup_anndata(
        adata,
        "size_factor",
        batch_key="batch",
        categorical_covariate_keys=["batch"],
        continuous_covariate_keys=["test"],
    )
    model = CellAssign(adata, marker_df)
    model.train(max_epochs=1)
    model.predict()
