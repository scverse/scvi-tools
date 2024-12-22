import pytest

from scvi.data import synthetic_iid
from scvi.external import Decipher


@pytest.fixture(scope="session")
def adata():
    adata = synthetic_iid()
    return adata


def test_decipher_train(adata):
    Decipher.setup_anndata(adata)
    model = Decipher(adata)
    model.train(
        max_epochs=2,
        check_val_every_n_epoch=1,
        train_size=0.5,
        early_stopping=True,
    )
    model.get_latent_representation(give_z=False)
    model.get_latent_representation(give_z=True)
