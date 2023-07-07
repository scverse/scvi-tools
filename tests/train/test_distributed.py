import pytest

import scvi


@pytest.mark.optional
def test_scvi_train_ddp():
    adata = scvi.data.synthetic_iid()
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata)

    model.train(
        max_epochs=1,
        accelerator="cpu",
        devices=2,
        strategy="ddp_find_unused_parameters_true",
    )
