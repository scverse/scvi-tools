import numpy as np
import pytest

import scvi


@pytest.mark.parametrize("n_batches", [1, 2, 3])
def test_scvi_batch_embedding(
    n_batches: int, batch_size: int = 10, batch_embedding_dim: int = 5
):
    adata = scvi.data.synthetic_iid(n_batches=n_batches, batch_size=batch_size)
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")

    if n_batches == 1:
        # cannot use batch embeddings with only one batch
        with pytest.raises(ValueError):
            scvi.model.SCVI(adata, use_batch_embedding=True)
        return

    model = scvi.model.SCVI(
        adata, use_batch_embedding=True, batch_embedding_dim=batch_embedding_dim
    )
    model.train(max_epochs=1)

    assert hasattr(model.module, "batch_embedding")
    assert model.module.batch_embedding is not None
    assert model.module.batch_embedding.weight.shape == (n_batches, batch_embedding_dim)

    batch_representation = model.get_batch_representation()
    assert isinstance(batch_representation, np.ndarray)
    assert batch_representation.shape == (n_batches, batch_embedding_dim)

    batch_representation = model.get_batch_representation(batch_keys=["batch_0"])
    assert isinstance(batch_representation, np.ndarray)
    assert batch_representation.shape == (1, batch_embedding_dim)

    batch_representation = model.get_batch_representation(
        batch_keys=["batch_0", "batch_1"]
    )
    assert isinstance(batch_representation, np.ndarray)
    assert batch_representation.shape == (2, batch_embedding_dim)

    with pytest.raises(ValueError):
        model.get_batch_representation(["not_a_batch_key"])

    model_no_batch = scvi.model.SCVI(adata)
    model_no_batch.train(max_epochs=1)

    with pytest.raises(NotImplementedError):
        _ = model_no_batch.get_batch_representation()
