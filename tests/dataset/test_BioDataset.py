import numpy as np
import random
import pandas as pd

from unittest import TestCase
from scvi.dataset._biodataset import BioDataset
from scvi.dataset._datasets import synthetic_iid
from scvi.dataset import setup_anndata
from scvi.dataset._anndata import get_from_registry


class TestBioDataset(TestCase):
    def test_setup_anndata(self,):
        # test regular setup
        adata = synthetic_iid()
        setup_anndata(
            adata,
            batch_key="batch",
            labels_key="labels",
            protein_expression_obsm_key="protein_expression",
            protein_names_uns_key="protein_names",
        )
        np.testing.assert_array_equal(
            get_from_registry(adata, "batch_indices"),
            np.array(adata.obs["batch"]).reshape((-1, 1)),
        )
        np.testing.assert_array_equal(
            get_from_registry(adata, "labels"),
            np.array(adata.obs["labels"].cat.codes).reshape((-1, 1)),
        )
        np.testing.assert_array_equal(get_from_registry(adata, "X"), adata.X)
        np.testing.assert_array_equal(
            get_from_registry(adata, "protein_expression"),
            adata.obsm["protein_expression"],
        )
        np.testing.assert_array_equal(
            adata.uns["scvi_summary_stats"]["protein_names"],
            adata.uns["protein_names"],
        )

        # test that error is thrown if its a view:
        adata = synthetic_iid()
        with self.assertRaises(ValueError):
            setup_anndata(adata[1])

        # If obsm is a df and protein_names_uns_key is None, protein names should be grabbed from column of df
        adata = synthetic_iid()
        new_protein_names = np.array(random.sample(range(100), 100)).astype("str")
        df = pd.DataFrame(
            adata.obsm["protein_expression"],
            index=adata.obs_names,
            columns=new_protein_names,
        )
        adata.obsm["protein_expression"] = df
        setup_anndata(adata, protein_expression_obsm_key="protein_expression")
        np.testing.assert_array_equal(
            adata.uns["scvi_summary_stats"]["protein_names"], new_protein_names
        )

        # test that X_layers_key is working properly
        adata = synthetic_iid()
        true_X = adata.X
        adata.layers["X"] = true_X
        adata.X = np.ones_like(adata.X)
        setup_anndata(adata, X_layers_key="X")
        np.testing.assert_array_equal(get_from_registry(adata, "X"), true_X)

        # test that it creates layers and batch if no layers_key is passed
        adata = synthetic_iid()
        setup_anndata(
            adata,
            protein_expression_obsm_key="protein_expression",
            protein_names_uns_key="protein_names",
        )
        np.testing.assert_array_equal(
            get_from_registry(adata, "batch_indices"), np.zeros((adata.shape[0], 1))
        )
        np.testing.assert_array_equal(
            get_from_registry(adata, "labels"), np.zeros((adata.shape[0], 1))
        )

    def test_BioDataset_getitem(self,):
        adata = synthetic_iid()
        setup_anndata(
            adata,
            batch_key="batch",
            labels_key="labels",
            protein_expression_obsm_key="protein_expression",
            protein_names_uns_key="protein_names",
        )
        # check that we can successfully pass in a list of tensors to get
        tensors_to_get = ["batch_indices", "local_l_var"]
        bd = BioDataset(adata, getitem_tensors=tensors_to_get)
        np.testing.assert_array_equal(tensors_to_get, list(bd[1].keys()))

        # check that we can successfully pass in a dict of tensors and their associated types
        bd = BioDataset(adata, getitem_tensors={"X": np.int, "local_l_var": np.float64})
        assert bd[1]["X"].dtype == np.int64
        assert bd[1]["local_l_var"].dtype == np.float64

        # check that by default we get all the registered tensors
        bd = BioDataset(adata)
        all_registered_tensors = list(adata.uns["scvi_data_registry"].keys())
        np.testing.assert_array_equal(all_registered_tensors, list(bd[1].keys()))
        assert bd[1]["X"].shape[0] == bd.n_genes
