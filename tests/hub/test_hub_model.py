import json
import os
from dataclasses import asdict

import numpy as np
import pytest

import scvi
from scvi.hub import HubMetadata, HubModel, HubModelCardHelper, get_models_df
from scvi.hub._constants import METADATA_FILE_NAME, MODEL_CARD_FILE_NAME


def prep_model():
    adata = scvi.data.synthetic_iid()
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata)
    model.train(1)
    return model


def test_hub_model_init(request, save_path):
    model = prep_model()
    test_save_path = os.path.join(save_path, request.node.name)
    model.save(test_save_path, overwrite=True, save_anndata=True)

    # no metadata, no model card
    with pytest.raises(ValueError) as e:
        HubModel(test_save_path)
    assert str(e.value) == "No metadata found"

    # metadata, no model card
    hm = HubMetadata("0.17.4", "0.8.0")
    with pytest.raises(ValueError) as e:
        HubModel(test_save_path, metadata=hm)
    assert str(e.value) == "No model card found"

    # model card, no metadata
    hmch = HubModelCardHelper.from_dir(
        test_save_path,
        license_info="cc-by-4.0",
        anndata_version="0.8.0",
        model_parent_module="other_module",
    )
    with pytest.raises(ValueError) as e:
        HubModel(test_save_path, model_card=hmch.model_card)
    assert str(e.value) == "No metadata found"

    # successful inits
    hmo = HubModel(test_save_path, metadata=hm, model_card=hmch.model_card)
    assert hmo.metadata is hm
    assert hmo.model_card is hmch.model_card

    hmo = HubModel(test_save_path, metadata=hm, model_card=hmch)
    assert hmo.metadata is hm
    assert hmo.model_card is hmch.model_card

    hmch.model_card.save(os.path.join(test_save_path, MODEL_CARD_FILE_NAME))
    with open(os.path.join(test_save_path, METADATA_FILE_NAME), "w") as fp:
        json.dump(asdict(hm), fp, indent=4)
    hmo = HubModel(test_save_path)
    assert hmo.metadata == hm
    assert hmo.model_card.content == hmch.model_card.content

    assert hmo._local_dir == test_save_path
    assert hmo._model is None
    assert hmo._adata is None
    assert hmo._adata_large is None


def test_hub_model_load(request, save_path):
    model = prep_model()
    test_save_path = os.path.join(save_path, request.node.name)
    model.save(test_save_path, overwrite=True, save_anndata=True)

    hm = HubMetadata("0.17.4", "0.8.0")
    hmch = HubModelCardHelper.from_dir(
        test_save_path,
        license_info="cc-by-4.0",
        anndata_version="0.8.0",
        model_parent_module="other_module",
    )

    hmo = HubModel(test_save_path, metadata=hm, model_card=hmch.model_card)
    assert isinstance(hmo.model, scvi.model.SCVI)
    assert isinstance(hmo.model.module, scvi.module.VAE)
    assert np.array_equal(hmo.adata.X, model.adata.X)
    assert hmo.adata.obs.equals(model.adata.obs)
    assert hmo.adata.var.equals(model.adata.var)
    assert hmo.adata_large is None

    # with a custom adata
    test_save_path = os.path.join(save_path, request.node.name + "_2")
    model.save(test_save_path, overwrite=True)
    hmo = HubModel(test_save_path, metadata=hm, model_card=hmch.model_card)
    adata2 = scvi.data.synthetic_iid()
    hmo.load_model(adata=adata2)
    assert isinstance(hmo.model, scvi.model.SCVI)
    assert isinstance(hmo.model.module, scvi.module.VAE)
    assert np.array_equal(hmo.model.adata.X, adata2.X)
    assert hmo.adata is None
    assert hmo.adata_large is None

    # no adata
    hmo = HubModel(test_save_path, metadata=hm, model_card=hmch.model_card)
    with pytest.raises(ValueError):
        hmo.model


@pytest.mark.internet
def test_hub_model_adata_large(request, save_path):
    large_data_url = "https://huggingface.co/datasets/scvi-tools/DATASET-FOR-UNIT-TESTING-1/resolve/main/adata.h5ad"
    model = prep_model()
    test_save_path = os.path.join(save_path, request.node.name)
    model.save(test_save_path, overwrite=True)

    hm = HubMetadata("0.17.4", "0.8.0", large_data_url=large_data_url)
    hmch = HubModelCardHelper.from_dir(
        test_save_path,
        license_info="cc-by-4.0",
        anndata_version="0.8.0",
        model_parent_module="other_module",
    )

    hmo = HubModel(test_save_path, metadata=hm, model_card=hmch.model_card)
    assert isinstance(hmo.model, scvi.model.SCVI)
    assert isinstance(hmo.model.module, scvi.module.VAE)
    assert np.array_equal(hmo.adata_large.X, model.adata.X)
    assert hmo.adata_large.obs.equals(model.adata.obs)
    assert hmo.adata_large.var.equals(model.adata.var)
    assert hmo.adata is None


@pytest.mark.internet
def test_hub_model_pull_from_hf(request, save_path):
    # the repo we are pulling from was populated with the contents of
    # `test_save_path` as below
    # model = prep_model()
    # test_save_path = os.path.join(save_path, request.node.name)
    # model.save(test_save_path, overwrite=True, save_anndata=True)
    # hm = HubMetadata("0.17.0", "0.8.0")
    # with open(os.path.join(test_save_path, METADATA_FILE_NAME), "w") as fp:
    #     json.dump(asdict(hm), fp, indent=4)

    hmo = HubModel.pull_from_huggingface_hub(
        repo_name="scvi-tools/MODEL-FOR-UNIT-TESTING-1"
    )
    assert hmo.metadata == HubMetadata("0.17.0", "0.8.0")
    assert hmo.model_card.content == "---\nlicense: cc-by-4.0\n---\n"
    assert isinstance(hmo.model, scvi.model.SCVI)
    assert isinstance(hmo.model.module, scvi.module.VAE)
    assert hmo.adata.shape == (400, 100)
    assert hmo.adata_large is None


# @pytest.mark.internet
# def test_hub_model_push_to_hf(request, save_path):
#     model = prep_model()
#     test_save_path = os.path.join(save_path, request.node.name)
#     model.save(test_save_path, overwrite=True, save_anndata=True)

#     hm = HubMetadata("0.17.0", "0.8.0")
#     with open(os.path.join(test_save_path, METADATA_FILE_NAME), "w") as fp:
#         json.dump(asdict(hm), fp, indent=4)

#     hmch = HubModelCardHelper.from_dir(
#         test_save_path,
#         license_info="cc-by-4.0",
#         anndata_version="0.8.0",
#     )

#     hmo = HubModel(test_save_path, metadata=hm, model_card=hmch.model_card)
#     # TODO replace this
#     token_path = "TBD"
#     repo_name = "scvi-tools/MODEL-FOR-UNIT-TESTING-2"
#     hmo.push_to_huggingface_hub(repo_name=repo_name, repo_token_path=token_path, repo_create=True)

#     # pull back down and validate
#     hmo = HubModel.pull_from_huggingface_hub(repo_name=repo_name)
#     assert hmo.metadata == hm
#     assert hmo.model_card.content == hmch.model_card.content
#     assert isinstance(hmo.model, scvi.model.SCVI)
#     assert isinstance(hmo.model.module, scvi.module.VAE)
#     assert np.array_equal(hmo.adata.X, model.adata.X)
#     assert hmo.adata.obs.equals(model.adata.obs)
#     assert hmo.adata.var.equals(model.adata.var)
#     assert hmo.adata_large is None

#     # delete the HF repo
#     repo_token = Path(token_path).read_text()
#     delete_repo(repo_name, token=repo_token)


@pytest.mark.internet
def test_hub_model_get_models_df():
    df = get_models_df()
    repo = "scvi-tools/test_model"
    assert repo in df.index
    assert df.loc[repo].author == "scvi-tools"
