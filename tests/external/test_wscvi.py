from scvi.external.wscvi import WVAE, WSCVI
from scvi.data import synthetic_iid


def test_wscvi():
    adata = synthetic_iid()
    model = WSCVI(adata=adata, n_latent=10, loss_type="IWELBO", n_particles=25)
    model.train(
        max_epochs=5,
    )
    pass