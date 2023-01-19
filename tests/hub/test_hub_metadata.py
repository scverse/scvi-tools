import os

import pytest

import scvi
from scvi.hub import HubMetadata, HubModelCardHelper


def prep_model():
    adata = scvi.data.synthetic_iid()
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata)
    model.train(1)
    return model


def test_hub_metadata(request, save_path):
    hm = HubMetadata("0.17.4", "0.8.0")
    assert hm.scvi_version == "0.17.4"
    assert hm.anndata_version == "0.8.0"
    assert hm.training_data_url is None
    assert hm.model_parent_module == "scvi.model"

    d = {
        "scvi_version": "0.15.4",
        "anndata_version": "0.8.1",
        "training_data_url": None,
        "model_parent_module": "bar",
    }
    hm = HubMetadata(**d)
    assert hm.scvi_version == "0.15.4"
    assert hm.anndata_version == "0.8.1"
    assert hm.training_data_url is None
    assert hm.model_parent_module == "bar"

    d = {
        "scvi_version": "0.15.4",
        "anndata_version": "0.8.1",
        "foo": "bar",
    }
    with pytest.raises(TypeError):
        hm = HubMetadata(**d)

    d = {
        "scvi_version": "0.15.4",
    }
    with pytest.raises(TypeError):
        hm = HubMetadata(**d)

    model = prep_model()
    test_save_path = os.path.join(save_path, request.node.name)
    model.save(test_save_path, overwrite=True)
    hm = HubMetadata.from_dir(
        test_save_path, anndata_version="0.9.0", model_parent_module="foo"
    )
    assert hm.scvi_version == scvi.__version__
    assert hm.anndata_version == "0.9.0"
    assert hm.training_data_url is None
    assert hm.model_parent_module == "foo"


def test_hub_modelcardhelper(request, save_path):
    model = prep_model()

    hmch = HubModelCardHelper(
        license_info="cc-by-4.0",
        model_cls_name="SCVI",
        model_init_params=model.init_params_,
        model_setup_anndata_args=model.adata_manager._get_setup_method_args()[
            "setup_args"
        ],
        model_summary_stats=model.summary_stats,
        model_data_registry=model.adata_manager.data_registry,
        scvi_version="0.17.8",
        anndata_version="0.8.0",
        tissues=["eye"],
    )

    assert hmch.license_info == "cc-by-4.0"
    assert hmch.model_cls_name == "SCVI"
    assert hmch.model_init_params == {
        "kwargs": {"model_kwargs": {}},
        "non_kwargs": {
            "n_hidden": 128,
            "n_latent": 10,
            "n_layers": 1,
            "dropout_rate": 0.1,
            "dispersion": "gene",
            "gene_likelihood": "zinb",
            "latent_distribution": "normal",
        },
    }
    assert hmch.model_setup_anndata_args == {
        "layer": None,
        "batch_key": None,
        "labels_key": None,
        "size_factor_key": None,
        "categorical_covariate_keys": None,
        "continuous_covariate_keys": None,
    }
    assert dict(hmch.model_summary_stats) == {
        "n_batch": 1,
        "n_cells": 400,
        "n_extra_categorical_covs": 0,
        "n_extra_continuous_covs": 0,
        "n_labels": 1,
        "n_vars": 100,
    }
    assert hmch.model_data_registry.keys() == ["X", "batch", "labels"]
    assert dict(hmch.model_data_registry["X"]) == {"attr_key": None, "attr_name": "X"}
    assert dict(hmch.model_data_registry["batch"]) == {
        "attr_key": "_scvi_batch",
        "attr_name": "obs",
    }
    assert dict(hmch.model_data_registry["labels"]) == {
        "attr_key": "_scvi_labels",
        "attr_name": "obs",
    }
    assert hmch.scvi_version == "0.17.8"
    assert hmch.anndata_version == "0.8.0"
    assert hmch.data_modalities == []
    assert hmch.tissues == ["eye"]
    assert hmch.data_is_annotated is None
    assert hmch.data_is_latent is None
    assert hmch.training_data_url is None
    assert hmch.training_code_url is None
    assert hmch.model_parent_module == "scvi.model"
    assert hmch.description == "To be added..."
    assert hmch.references == "To be added..."
    assert hmch.model_card.data.to_dict() == {
        "license": "cc-by-4.0",
        "library_name": "scvi-tools",
        "tags": [
            "biology",
            "genomics",
            "single-cell",
            "model_cls_name:SCVI",
            "scvi_version:0.17.8",
            "anndata_version:0.8.0",
            "tissue:eye",
        ],
    }

    test_save_path = os.path.join(save_path, request.node.name)
    model.save(test_save_path, overwrite=True, save_anndata=True)
    hmch = HubModelCardHelper.from_dir(
        test_save_path,
        license_info="cc-by-4.0",
        anndata_version="0.8.0",
        model_parent_module="other_module",
    )

    assert hmch.license_info == "cc-by-4.0"
    assert hmch.model_cls_name == "SCVI"
    assert hmch.model_init_params == model.init_params_
    assert (
        hmch.model_setup_anndata_args
        == model.adata_manager._get_setup_method_args()["setup_args"]
    )
    assert hmch.model_summary_stats == dict(model.summary_stats)
    assert hmch.model_data_registry == dict(model.adata_manager.data_registry)
    assert hmch.scvi_version == scvi.__version__
    assert hmch.anndata_version == "0.8.0"
    assert hmch.data_modalities == []
    assert hmch.tissues == []
    assert hmch.data_is_annotated is None
    assert hmch.data_is_latent is False
    assert hmch.training_data_url is None
    assert hmch.training_code_url is None
    assert hmch.model_parent_module == "other_module"
    assert hmch.description == "To be added..."
    assert hmch.references == "To be added..."
    assert hmch.model_card.data.to_dict() == {
        "license": "cc-by-4.0",
        "library_name": "scvi-tools",
        "tags": [
            "biology",
            "genomics",
            "single-cell",
            "model_cls_name:SCVI",
            f"scvi_version:{scvi.__version__}",
            "anndata_version:0.8.0",
        ],
    }
