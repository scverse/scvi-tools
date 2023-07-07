import pytest

import scvi


@pytest.mark.cuda
def test_scvi_train_ddp(devices: int = -1):
    adata = scvi.data.synthetic_iid()
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata)

    model.train(
        max_epochs=1,
        check_val_every_n_epoch=1,
        accelerator="gpu",
        devices=devices,
        strategy="ddp_find_unused_parameters_true",
    )

    assert model.is_trained
    assert len(model.history) > 0
