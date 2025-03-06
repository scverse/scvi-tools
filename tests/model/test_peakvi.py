import os
import pickle

import numpy as np
import pytest
import torch

from scvi.data import synthetic_iid
from scvi.model import PEAKVI
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
        legacy_save(model, legacy_save_path, overwrite=True, save_anndata=True, prefix=prefix)
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

    test_save_load_model(PEAKVI, adata, save_path, prefix=f"{PEAKVI.__name__}_")


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
    vae.get_normalized_accessibility()
    vae.get_normalized_accessibility(normalize_cells=True)
    vae.get_normalized_accessibility(normalize_regions=True)
    vae.get_library_size_factors()
    vae.get_region_factors()
    vae.get_reconstruction_error(indices=vae.validation_indices)
    vae.get_latent_representation()
    vae.differential_accessibility(groupby="labels", group1="label_1")
    vae.get_normalized_expression()
    vae.get_normalized_expression(transform_batch="batch_1")
    vae.get_normalized_expression(n_samples=2)


def single_pass_for_online_update(model):
    dl = model._make_data_loader(model.adata, indices=range(0, 10))
    for tensors in dl:
        _, _, scvi_loss = model.module(tensors)
    scvi_loss.loss.backward()


def test_peakvi_online_update(save_path):
    n_latent = 5
    adata1 = synthetic_iid()
    PEAKVI.setup_anndata(adata1, batch_key="batch", labels_key="labels")
    model = PEAKVI(adata1, n_latent=n_latent)
    model.train(1, save_best=False)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # also test subset var option
    adata2 = synthetic_iid()
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])

    model2 = PEAKVI.load_query_data(adata2, dir_path)
    model2.train(max_epochs=1, weight_decay=0.0, save_best=False)
    model2.get_latent_representation()
    single_pass_for_online_update(model2)

    # encoder linear layer equal for peak features
    one = (
        model.module.z_encoder.encoder.fc_layers[0][0]
        .weight.detach()
        .cpu()
        .numpy()[:, : adata1.shape[1]]
    )
    two = (
        model2.module.z_encoder.encoder.fc_layers[0][0]
        .weight.detach()
        .cpu()
        .numpy()[:, : adata1.shape[1]]
    )
    np.testing.assert_equal(one, two)
    assert (
        np.sum(
            model2.module.z_encoder.encoder.fc_layers[0][0]
            .weight.grad.cpu()
            .numpy()[:, : adata1.shape[1]]
        )
        == 0
    )

    # test options
    adata1 = synthetic_iid()
    PEAKVI.setup_anndata(adata1, batch_key="batch", labels_key="labels")
    model = PEAKVI(
        adata1,
        n_latent=n_latent,
        encode_covariates=True,
    )
    model.train(1, check_val_every_n_epoch=1, save_best=False)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    adata2 = synthetic_iid()
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])

    model2 = PEAKVI.load_query_data(adata2, dir_path, freeze_expression=True)
    model2.train(max_epochs=1, weight_decay=0.0, save_best=False)
    # deactivate no grad decorator
    model2.get_latent_representation()
    # pytorch lightning zeros the grad, so this will get a grad to inspect
    single_pass_for_online_update(model2)
    grad = model2.module.z_encoder.encoder.fc_layers[0][0].weight.grad.cpu().numpy()
    # expression part has zero grad
    assert np.sum(grad[:, :-4]) == 0
    # categorical part has non-zero grad
    assert np.count_nonzero(grad[:, -4:]) > 0

    # do not freeze expression
    model3 = PEAKVI.load_query_data(
        adata2,
        dir_path,
        freeze_expression=False,
        freeze_decoder_first_layer=False,
    )
    model3.train(max_epochs=1, save_best=False, weight_decay=0.0)
    model3.get_latent_representation()
    single_pass_for_online_update(model3)
    grad = model3.module.z_encoder.encoder.fc_layers[0][0].weight.grad.cpu().numpy()
    # linear layer weight in encoder layer has non-zero grad
    assert np.count_nonzero(grad[:, :-4]) != 0
    grad = model3.module.z_decoder.px_decoder.fc_layers[0][0].weight.grad.cpu().numpy()
    # linear layer weight in decoder layer has non-zero grad
    assert np.count_nonzero(grad[:, :-4]) != 0
