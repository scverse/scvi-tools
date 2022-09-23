import scipy

from scvi.data import synthetic_iid
from scvi.external import SCAR


def test_scar():
    n_latent = 5
    adata = synthetic_iid()
    adata.X = scipy.sparse.csr_matrix(adata.X)
    SCAR.setup_anndata(adata)

    _ = SCAR.get_ambient_profile(adata, adata, prob=0.0, iterations=1, sample=100)
    model = SCAR(adata, ambient_profile=None, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)
    model.get_elbo()
    model.get_latent_representation()
    model.get_marginal_ll(n_mc_samples=5)
    model.get_reconstruction_error()
    model.get_denoised_counts(adata, n_samples=1)
    model.history

    # tests __repr__
    print(model)
