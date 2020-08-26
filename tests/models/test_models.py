import pandas as pd

from unittest import TestCase
from scvi.dataset import synthetic_iid, transfer_anndata_setup
from scvi.models import SCVI, SCANVI, GIMVI, TOTALVI, LinearSCVI


class TestModels(TestCase):
    def test_SCVI(self):
        adata = synthetic_iid()
        model = SCVI(adata, n_latent=10)
        model.train(1)
        z = model.get_latent_representation()
        assert z.shape == (adata.shape[0], 10)
        model.get_elbo()
        model.get_marginal_ll()
        model.get_reconstruction_error()
        model.get_normalized_expression()

        adata2 = synthetic_iid()
        model.get_elbo(adata2)
        model.get_marginal_ll(adata2)
        model.get_reconstruction_error(adata2)
        model.get_normalized_expression(adata2)

        # test transfer_anndata_setup
        adata2 = synthetic_iid(run_setup_anndata=False)
        transfer_anndata_setup(adata, adata2)
        model.get_elbo(adata2)

        # test automatic transfer_anndata_setup
        adata = synthetic_iid()
        model = SCVI(adata)
        adata2 = synthetic_iid(run_setup_anndata=False)
        model.get_elbo(adata2)

        # test that we catch incorrect mappings
        adata = synthetic_iid()
        adata2 = synthetic_iid(run_setup_anndata=False)
        transfer_anndata_setup(adata, adata2)
        adata2.uns["_scvi"]["categorical_mappings"]["_scvi_labels"][
            "mapping"
        ] = pd.Index(data=["undefined_1", "undefined_0", "undefined_2"])
        with self.assertRaises(AssertionError):
            model.get_elbo(adata2)

    def test_SCANVI(self):
        adata = synthetic_iid()
        model = SCANVI(adata, n_latent=10)
        model.train(1)

    def test_LinearSCVI(self):
        adata = synthetic_iid()
        model = LinearSCVI(adata, n_latent=10)
        model.train(1)

    def test_GIMVI(self):
        adata = synthetic_iid()
        model = GIMVI(adata, n_latent=10)
        model.train(1)

    def test_TOTALVI(self):
        adata = synthetic_iid()
        n_obs = adata.n_obs
        n_vars = adata.n_vars
        n_proteins = adata.obsm["protein_expression"].shape[1]
        n_latent = 10

        model = TOTALVI(adata, n_latent=n_latent)
        model.train(1)
        z = model.get_latent_representation()
        assert z.shape == (n_obs, n_latent)
        model.get_elbo()
        model.get_marginal_ll()
        model.get_reconstruction_error()
        model.get_normalized_expression()
        model.get_latent_library_size()
        model.get_protein_foreground_probability()
        post_pred = model.posterior_predictive_sample(n_samples=2)
        assert post_pred.shape == (n_obs, n_vars + n_proteins, 2)
        post_pred = model.posterior_predictive_sample(n_samples=1)
        assert post_pred.shape == (n_obs, n_vars + n_proteins)
        feature_correlation_matrix = model.get_feature_correlation_matrix(
            correlation_type="pearson"
        )
        assert feature_correlation_matrix.shape == (
            n_vars + n_proteins,
            n_vars + n_proteins,
        )
        # model.get_likelihood_parameters()

        adata2 = synthetic_iid()
        model.get_elbo(adata2)
        model.get_marginal_ll(adata2)
        model.get_reconstruction_error(adata2)

        norm_exp = model.get_normalized_expression(adata2, indices=[1, 2, 3])
        assert norm_exp[0].shape == (3, adata2.n_vars)
        assert norm_exp[1].shape == (3, adata2.obsm["protein_expression"].shape[1])

        latent_lib_size = model.get_latent_library_size(adata2, indices=[1, 2, 3])
        assert latent_lib_size.shape == (3, 1)

        pro_foreground_prob = model.get_protein_foreground_probability(
            adata2, indices=[1, 2, 3], protein_list=["1", "2"]
        )
        assert pro_foreground_prob.shape == (3, 2)
        model.posterior_predictive_sample(adata2)
        # model.get_feature_correlation_matrix(adata2)
        # model.get_likelihood_parameters(adata2)

        # test transfer_anndata_setup
        adata2 = synthetic_iid(run_setup_anndata=False)
        transfer_anndata_setup(adata, adata2)
        model.get_elbo(adata2)

        # test automatic transfer_anndata_setup
        adata = synthetic_iid()
        model = SCVI(adata)
        adata2 = synthetic_iid(run_setup_anndata=False)
        model.get_elbo(adata2)

        # test that we catch incorrect mappings
        adata = synthetic_iid()
        adata2 = synthetic_iid(run_setup_anndata=False)
        transfer_anndata_setup(adata, adata2)
        adata2.uns["_scvi"]["categorical_mappings"]["_scvi_labels"][
            "mapping"
        ] = pd.Index(data=["undefined_1", "undefined_0", "undefined_2"])
        with self.assertRaises(AssertionError):
            model.get_elbo(adata2)
