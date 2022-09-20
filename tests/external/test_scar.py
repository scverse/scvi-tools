from scvi.data import synthetic_iid
from scvi.external import SCAR


def test_scar():
    n_latent = 5
    adata = synthetic_iid()
    SCAR.setup_anndata(adata, batch_key="batch", labels_key="labels")

    model = SCAR(adata, ambient_profile=None, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)
    model.get_elbo()
    model.get_latent_representation()
    model.get_marginal_ll(n_mc_samples=5)
    model.get_reconstruction_error()
    # testing both flavors of count removal
    model.get_denoised_counts(adata, n_samples=1)
    model.get_denoised_counts(adata, n_samples=1, flavor="remove_ambient_counts")
    model.history

    # tests __repr__
    print(model)
