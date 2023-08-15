import numpy as np

import scvi


def test_scvi_batch_embedding(n_batches: int = 3, batch_embedding_dim: int = 5):
    adata = scvi.data.synthetic_iid(n_batches=n_batches)
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")

    model = scvi.model.SCVI(
        adata, use_batch_embedding=True, batch_embedding_dim=batch_embedding_dim
    )
    model.train(max_epochs=1)

    batch_representation = model.get_batch_representation()
    assert isinstance(batch_representation, np.ndarray)
    assert batch_representation.shape == (n_batches, batch_embedding_dim)
