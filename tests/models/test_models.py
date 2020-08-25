import numpy as np
import pandas as pd

from unittest import TestCase
from scvi.dataset import synthetic_iid, transfer_anndata_setup
from scvi.models import SCVI, SCANVI, GIMVI, TOTALVI, LinearSCVI


class TestSCVI(TestCase):
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

    def test_DE(self):
        pass
