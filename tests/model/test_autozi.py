import os
import pickle

import numpy as np
import pytest
import torch

from scvi.data import synthetic_iid
from scvi.model import AUTOZI
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

    save_path = os.path.join(save_path, "tmp")
    adata = synthetic_iid()

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
    tmp_adata = synthetic_iid(n_genes=200)
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
    legacy_save(model, legacy_save_path, overwrite=True, save_anndata=True, prefix=prefix)
    with pytest.raises(ValueError):
        AUTOZI.load(legacy_save_path, adata=adata, prefix=prefix)
    AUTOZI.convert_legacy_save(legacy_save_path, legacy_save_path, overwrite=True, prefix=prefix)
    m = AUTOZI.load(legacy_save_path, adata=adata, prefix=prefix)
    m.train(1)


def test_autozi():
    data = synthetic_iid(
        n_batches=2,
    )
    AUTOZI.setup_anndata(
        data,
        batch_key="batch",
        labels_key="labels",
    )

    for disp_zi in ["gene", "gene-label"]:
        autozivae = AUTOZI(
            data,
            dispersion=disp_zi,
            zero_inflation=disp_zi,
        )
        autozivae.train(1, plan_kwargs={"lr": 1e-2}, check_val_every_n_epoch=1)
        assert len(autozivae.history["elbo_train"]) == 1
        assert len(autozivae.history["elbo_validation"]) == 1
        autozivae.get_elbo(indices=autozivae.validation_indices)
        autozivae.get_reconstruction_error(indices=autozivae.validation_indices)
        autozivae.get_marginal_ll(indices=autozivae.validation_indices, n_mc_samples=3)
        autozivae.get_alphas_betas()
        autozivae.get_normalized_expression()
        autozivae.get_normalized_expression(transform_batch="batch_1")
        autozivae.get_normalized_expression(n_samples=2)

    # Model library size.
    for disp_zi in ["gene", "gene-label"]:
        autozivae = AUTOZI(
            data,
            dispersion=disp_zi,
            zero_inflation=disp_zi,
            use_observed_lib_size=False,
        )
        autozivae.train(1, plan_kwargs={"lr": 1e-2}, check_val_every_n_epoch=1)
        assert hasattr(autozivae.module, "library_log_means")
        assert hasattr(autozivae.module, "library_log_vars")
        assert len(autozivae.history["elbo_train"]) == 1
        assert len(autozivae.history["elbo_validation"]) == 1
        autozivae.get_elbo(indices=autozivae.validation_indices)
        autozivae.get_reconstruction_error(indices=autozivae.validation_indices)
        autozivae.get_marginal_ll(indices=autozivae.validation_indices, n_mc_samples=3)
        autozivae.get_alphas_betas()
        autozivae.get_normalized_expression()
        autozivae.get_normalized_expression(transform_batch="batch_1")
        autozivae.get_normalized_expression(n_samples=2)
