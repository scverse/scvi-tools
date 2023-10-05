import os
import pickle

import numpy as np
import pytest
import torch

import scvi
from scvi.data import synthetic_iid
from scvi.model import (
    AUTOZI,
    PEAKVI,
    SCANVI,
    SCVI,
    TOTALVI,
    LinearSCVI,
)
from scvi.utils import attrdict


def test_saving_and_loading(save_path):
    def legacy_save(
        model,
        dir_path,
        prefix=None,
        overwrite=False,
        save_anndata=False,
        **anndata_write_kwargs,
    ):
        if not os.path.exists(dir_path) or overwrite:
            os.makedirs(dir_path, exist_ok=overwrite)
        else:
            raise ValueError(
                f"{dir_path} already exists. Please provide an unexisting directory for saving."
            )

        file_name_prefix = prefix or ""

        if save_anndata:
            model.adata.write(
                os.path.join(dir_path, f"{file_name_prefix}adata.h5ad"),
                **anndata_write_kwargs,
            )

        model_save_path = os.path.join(dir_path, f"{file_name_prefix}model_params.pt")
        attr_save_path = os.path.join(dir_path, f"{file_name_prefix}attr.pkl")
        varnames_save_path = os.path.join(dir_path, f"{file_name_prefix}var_names.csv")

        torch.save(model.module.state_dict(), model_save_path)

        var_names = model.adata.var_names.astype(str)
        var_names = var_names.to_numpy()
        np.savetxt(varnames_save_path, var_names, fmt="%s")

        # get all the user attributes
        user_attributes = model._get_user_attributes()
        # only save the public attributes with _ at the very end
        user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}

        with open(attr_save_path, "wb") as f:
            pickle.dump(user_attributes, f)

    def test_save_load_model(cls, adata, save_path, prefix=None):
        if cls is TOTALVI:
            cls.setup_anndata(
                adata,
                batch_key="batch",
                protein_expression_obsm_key="protein_expression",
                protein_names_uns_key="protein_names",
            )
        else:
            cls.setup_anndata(adata, batch_key="batch", labels_key="labels")
        model = cls(adata, latent_distribution="normal")
        model.train(1, train_size=0.2)
        z1 = model.get_latent_representation(adata)
        test_idx1 = model.validation_indices
        model.save(save_path, overwrite=True, save_anndata=True, prefix=prefix)
        model.view_setup_args(save_path, prefix=prefix)
        model = cls.load(save_path, prefix=prefix)
        model.get_latent_representation()

        # Load with mismatched genes.
        tmp_adata = synthetic_iid(
            n_genes=200,
        )
        with pytest.raises(ValueError):
            cls.load(save_path, adata=tmp_adata, prefix=prefix)

        # Load with different batches.
        tmp_adata = synthetic_iid()
        tmp_adata.obs["batch"] = tmp_adata.obs["batch"].cat.rename_categories(
            ["batch_2", "batch_3"]
        )
        with pytest.raises(ValueError):
            cls.load(save_path, adata=tmp_adata, prefix=prefix)

        model = cls.load(save_path, adata=adata, prefix=prefix)
        assert "batch" in model.adata_manager.data_registry
        assert model.adata_manager.data_registry.batch == attrdict(
            {"attr_name": "obs", "attr_key": "_scvi_batch"}
        )

        z2 = model.get_latent_representation()
        test_idx2 = model.validation_indices
        np.testing.assert_array_equal(z1, z2)
        np.testing.assert_array_equal(test_idx1, test_idx2)
        assert model.is_trained is True

        # Test legacy loading
        legacy_save_path = os.path.join(save_path, "legacy/")
        legacy_save(
            model, legacy_save_path, overwrite=True, save_anndata=True, prefix=prefix
        )
        with pytest.raises(ValueError):
            cls.load(legacy_save_path, adata=adata, prefix=prefix)
        cls.convert_legacy_save(
            legacy_save_path,
            legacy_save_path,
            overwrite=True,
            prefix=prefix,
        )
        m = cls.load(legacy_save_path, adata=adata, prefix=prefix)
        m.train(1)

    save_path = os.path.join(save_path, "tmp")
    adata = synthetic_iid()

    for cls in [SCVI, LinearSCVI, TOTALVI, PEAKVI]:
        test_save_load_model(cls, adata, save_path, prefix=f"{cls.__name__}_")

    # AUTOZI
    prefix = "AUTOZI_"
    AUTOZI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = AUTOZI(adata, latent_distribution="normal")
    model.train(1, train_size=0.5)
    ab1 = model.get_alphas_betas()
    model.save(save_path, overwrite=True, save_anndata=True, prefix=prefix)
    model.view_setup_args(save_path, prefix=prefix)
    model = AUTOZI.load(save_path, prefix=prefix)
    model.get_latent_representation()
    tmp_adata = scvi.data.synthetic_iid(n_genes=200)
    with pytest.raises(ValueError):
        AUTOZI.load(save_path, adata=tmp_adata, prefix=prefix)
    model = AUTOZI.load(save_path, adata=adata, prefix=prefix)
    assert "batch" in model.adata_manager.data_registry
    assert model.adata_manager.data_registry.batch == attrdict(
        {"attr_name": "obs", "attr_key": "_scvi_batch"}
    )

    ab2 = model.get_alphas_betas()
    np.testing.assert_array_equal(ab1["alpha_posterior"], ab2["alpha_posterior"])
    np.testing.assert_array_equal(ab1["beta_posterior"], ab2["beta_posterior"])
    assert model.is_trained is True

    # Test legacy loading
    legacy_save_path = os.path.join(save_path, "legacy/")
    legacy_save(
        model, legacy_save_path, overwrite=True, save_anndata=True, prefix=prefix
    )
    with pytest.raises(ValueError):
        AUTOZI.load(legacy_save_path, adata=adata, prefix=prefix)
    AUTOZI.convert_legacy_save(
        legacy_save_path, legacy_save_path, overwrite=True, prefix=prefix
    )
    m = AUTOZI.load(legacy_save_path, adata=adata, prefix=prefix)
    m.train(1)

    # SCANVI
    prefix = "SCANVI_"
    SCANVI.setup_anndata(adata, "labels", "label_0", batch_key="batch")
    model = SCANVI(adata)
    model.train(max_epochs=1, train_size=0.5)
    p1 = model.predict()
    model.save(save_path, overwrite=True, save_anndata=True, prefix=prefix)
    model.view_setup_args(save_path, prefix=prefix)
    model = SCANVI.load(save_path, prefix=prefix)
    model.get_latent_representation()
    tmp_adata = scvi.data.synthetic_iid(n_genes=200)
    with pytest.raises(ValueError):
        SCANVI.load(save_path, adata=tmp_adata, prefix=prefix)
    model = SCANVI.load(save_path, adata=adata, prefix=prefix)
    assert "batch" in model.adata_manager.data_registry
    assert model.adata_manager.data_registry.batch == attrdict(
        {"attr_name": "obs", "attr_key": "_scvi_batch"}
    )

    p2 = model.predict()
    np.testing.assert_array_equal(p1, p2)
    assert model.is_trained is True

    # Test legacy loading
    legacy_save_path = os.path.join(save_path, "legacy/")
    legacy_save(
        model, legacy_save_path, overwrite=True, save_anndata=True, prefix=prefix
    )
    with pytest.raises(ValueError):
        SCANVI.load(legacy_save_path, adata=adata, prefix=prefix)
    SCANVI.convert_legacy_save(
        legacy_save_path, legacy_save_path, overwrite=True, prefix=prefix
    )
    m = SCANVI.load(legacy_save_path, adata=adata, prefix=prefix)
    m.train(1)


def test_peakvi():
    data = synthetic_iid()
    PEAKVI.setup_anndata(
        data,
        batch_key="batch",
    )
    vae = PEAKVI(
        data,
        model_depth=False,
    )
    vae.train(1, save_best=False)
    vae = PEAKVI(
        data,
        region_factors=False,
    )
    vae.train(1, save_best=False)
    vae = PEAKVI(
        data,
    )
    vae.train(3)
    vae.get_elbo(indices=vae.validation_indices)
    vae.get_accessibility_estimates()
    vae.get_accessibility_estimates(normalize_cells=True)
    vae.get_accessibility_estimates(normalize_regions=True)
    vae.get_library_size_factors()
    vae.get_region_factors()
    vae.get_reconstruction_error(indices=vae.validation_indices)
    vae.get_latent_representation()
    vae.differential_accessibility(groupby="labels", group1="label_1")
