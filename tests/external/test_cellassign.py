import numpy as np
import pandas as pd

from scvi.data import synthetic_iid
from scvi.external import CellAssign


def get_test_adata_marker_mat():
    adata = synthetic_iid(
        n_labels=5,
    )
    adata.obs["size_factor"] = adata.X.sum(1)

    marker_df = pd.DataFrame(data=np.random.randint(2, size=(100, 5)))
    marker_df.index = marker_df.index.map(str)

    return adata, marker_df


def test_cellassign(save_path):
    adata, marker_df = get_test_adata_marker_mat()
    CellAssign.setup_anndata(
        adata,
        "size_factor",
        batch_key="batch",
    )
    model = CellAssign(adata, marker_df)
    model.train(max_epochs=1)
    model.predict()


def test_cellassign_covariates(save_path):
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
