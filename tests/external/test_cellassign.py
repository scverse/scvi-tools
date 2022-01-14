import numpy as np
import pandas as pd

from scvi.data import synthetic_iid
from scvi.external import CellAssign


def test_cellassign(save_path):
    adata = synthetic_iid(
        n_labels=5,
    )
    adata.obs["size_factor"] = adata.X.sum(1)
    CellAssign.setup_anndata(
        adata,
        "size_factor",
        batch_key="batch",
    )
    marker_df = pd.DataFrame(data=np.random.randint(2, size=(100, 5)))
    marker_df.index = marker_df.index.map(str)

    model = CellAssign(adata, marker_df)
    model.train(max_epochs=1)
    model.predict()
