import json
import os
from dataclasses import asdict

import anndata
import numpy as np
import pytest

import scvi
from scvi.criticism import create_criticism_report
from scvi.data import synthetic_iid
from scvi.hub import HubMetadata, HubModel, HubModelCardHelper
from scvi.hub._constants import _SCVI_HUB


def prep_model() -> scvi.model.SCVI:
    adata = scvi.data.synthetic_iid()
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata)
    model.train(1)
    return model


def prep_scvi_hub_model(save_path: str) -> HubModel:
    model = prep_model()
    model_path = os.path.join(save_path, "test_scvi")
    model.save(model_path, save_anndata=True, overwrite=True)
    create_criticism_report(model, save_folder=model_path, n_samples=2)

    metadata = HubMetadata.from_dir(model_path, anndata_version=anndata.__version__)
    desc = "scVI model trained on synthetic IID data and uploaded with the full training data."
    card = HubModelCardHelper.from_dir(
        model_path,
        license_info="cc-by-4.0",
        anndata_version=anndata.__version__,
        data_modalities=["rna"],
        data_is_annotated=False,
        description=desc,
    )
    return HubModel(model_path, metadata=metadata, model_card=card)


def prep_scvi_no_anndata_hub_model(save_path: str) -> HubModel:
    model = prep_model()
    model_path = os.path.join(save_path, "test_scvi_no_anndata")
    model.save(model_path, save_anndata=False, overwrite=True)

    metadata = HubMetadata.from_dir(model_path, anndata_version=anndata.__version__)
    card = HubModelCardHelper.from_dir(
        model_path,
        license_info="cc-by-4.0",
        anndata_version=anndata.__version__,
        data_modalities=["rna"],
        data_is_annotated=False,
        description="scVI model trained on synthetic IID data and uploaded with no data.",
    )
    return HubModel(model_path, metadata=metadata, model_card=card)


def prep_scvi_minified_hub_model(save_path: str) -> HubModel:
    model = prep_model()
    model_path = os.path.join(save_path, "test_scvi_minified")
    model.save(model_path, save_anndata=False, overwrite=True)

    qzm, qzv = model.get_latent_representation(return_dist=True, give_mean=False)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv
    model.minify_adata()
    model_path = os.path.join(save_path, "test_scvi_minified")
    model.save(model_path, save_anndata=True, overwrite=True)

    metadata = HubMetadata.from_dir(model_path, anndata_version=anndata.__version__)
    desc = "scVI model trained on synthetic IID data and uploaded with the minified data."
    card = HubModelCardHelper.from_dir(
        model_path,
        license_info="cc-by-4.0",
        anndata_version=anndata.__version__,
        data_modalities=["rna"],
        data_is_annotated=False,
        description=desc,
    )
    return HubModel(model_path, metadata=metadata, model_card=card)


def test_hub_model_init(request, save_path):
    model = prep_model()
    test_save_path = os.path.join(save_path, request.node.name)
    model.save(test_save_path, overwrite=True, save_anndata=True)

    # no metadata, no model card
    with pytest.raises(ValueError) as e:
        HubModel(test_save_path)
    assert str(e.value) == "No metadata found"

    # metadata, no model card
    hm = HubMetadata(scvi.__version__, anndata.__version__, "SCVI")
    with pytest.raises(ValueError) as e:
        HubModel(test_save_path, metadata=hm)
    assert str(e.value) == "No model card found"

    # model card, no metadata
    hmch = HubModelCardHelper.from_dir(
        test_save_path,
        license_info="cc-by-4.0",
        anndata_version=anndata.__version__,
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

    hmch.model_card.save(os.path.join(test_save_path, _SCVI_HUB.MODEL_CARD_FILE_NAME))
    with open(os.path.join(test_save_path, _SCVI_HUB.METADATA_FILE_NAME), "w") as fp:
        json.dump(asdict(hm), fp, indent=4)
    hmo = HubModel(test_save_path)
    assert hmo.metadata == hm
    assert hmo.model_card.content == hmch.model_card.content

    assert hmo._local_dir == test_save_path
    assert hmo._model is None
    assert hmo._adata is None
    assert hmo._large_training_adata is None


def test_hub_model_load(request, save_path):
    model = prep_model()
    test_save_path = os.path.join(save_path, request.node.name)
    model.save(test_save_path, overwrite=True, save_anndata=True)

    hm = HubMetadata(scvi.__version__, anndata.__version__, "SCVI")
    hmch = HubModelCardHelper.from_dir(
        test_save_path,
        license_info="cc-by-4.0",
        anndata_version=anndata.__version__,
        model_parent_module="other_module",
    )

    hmo = HubModel(test_save_path, metadata=hm, model_card=hmch.model_card)
    assert isinstance(hmo.model, scvi.model.SCVI)
    assert isinstance(hmo.model.module, scvi.module.VAE)
    assert np.array_equal(hmo.adata.X, model.adata.X)
    assert hmo.adata.obs.equals(model.adata.obs)
    assert hmo.adata.var.equals(model.adata.var)
    assert hmo.large_training_adata is None

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
    assert hmo.large_training_adata is None

    # no adata
    hmo = HubModel(test_save_path, metadata=hm, model_card=hmch.model_card)
    with pytest.raises(ValueError):
        print(hmo.model)


@pytest.mark.parametrize("save_anndata", [True, False])
def test_hub_model_save(save_anndata: bool, save_path: str):
    model = prep_model()
    model_path = os.path.join(save_path, "model_save")
    model.save(model_path, save_anndata=save_anndata, overwrite=True)

    metadata = HubMetadata.from_dir(model_path, anndata_version=anndata.__version__)
    card = HubModelCardHelper.from_dir(
        model_path,
        license_info="cc-by-4.0",
        anndata_version=anndata.__version__,
        data_modalities=["rna"],
        data_is_annotated=False,
    )
    hub_model = HubModel(model_path, metadata=metadata, model_card=card)
    hub_model.save(overwrite=True)

    card_path = os.path.join(model_path, _SCVI_HUB.MODEL_CARD_FILE_NAME)
    assert os.path.exists(card_path)
    assert os.path.isfile(card_path)
    metadata_path = os.path.join(model_path, _SCVI_HUB.METADATA_FILE_NAME)
    assert os.path.exists(metadata_path)
    assert os.path.isfile(metadata_path)

    with pytest.raises(FileExistsError):
        hub_model.save(overwrite=False)

    hub_model.save(overwrite=True)


def test_prep_scvi_hub_model(save_path: str) -> HubModel:
    prep_scvi_hub_model(save_path)


@pytest.mark.private
def test_hub_model_large_training_adata(request, save_path):
    training_data_url = "https://huggingface.co/datasets/scvi-tools/DATASET-FOR-UNIT-TESTING-1/resolve/main/adata.h5ad"
    model = prep_model()
    test_save_path = os.path.join(save_path, request.node.name)
    model.save(test_save_path, overwrite=True)

    hm = HubMetadata(
        scvi.__version__, anndata.__version__, "SCVI", training_data_url=training_data_url
    )
    hmch = HubModelCardHelper.from_dir(
        test_save_path,
        license_info="cc-by-4.0",
        anndata_version=anndata.__version__,
        model_parent_module="other_module",
    )

    hmo = HubModel(test_save_path, metadata=hm, model_card=hmch.model_card)
    assert isinstance(hmo.model, scvi.model.SCVI)
    assert isinstance(hmo.model.module, scvi.module.VAE)
    assert hmo.large_training_adata is hmo.model.adata
    assert hmo.adata is None


@pytest.mark.private
def test_hub_model_create_repo_hf(save_path: str):
    from huggingface_hub import delete_repo, repo_exists

    if repo_exists("scvi-tools/test-scvi-create"):
        delete_repo("scvi-tools/test-scvi-create", token=os.environ.get("HF_API_TOKEN", None))

    hub_model = prep_scvi_hub_model(save_path)
    hub_model.push_to_huggingface_hub(
        "scvi-tools/test-scvi-create",
        os.environ.get("HF_API_TOKEN", None),
        collection_name="test",
        repo_create=True,
    )
    delete_repo("scvi-tools/test-scvi-create", token=os.environ.get("HF_API_TOKEN", None))


@pytest.mark.private
def test_hub_model_push_to_hf(save_path: str):
    hub_model = prep_scvi_hub_model(save_path)
    hub_model.push_to_huggingface_hub(
        "scvi-tools/test-scvi",
        os.environ.get("HF_API_TOKEN", None),
        repo_create=False,
        collection_name="test",
    )

    hub_model = prep_scvi_no_anndata_hub_model(save_path)
    hub_model.push_to_huggingface_hub(
        "scvi-tools/test-scvi-no-anndata",
        os.environ.get("HF_API_TOKEN", None),
        repo_create=False,
        push_anndata=False,
        collection_name="test",
    )

    hub_model = prep_scvi_minified_hub_model(save_path)
    hub_model.push_to_huggingface_hub(
        "scvi-tools/test-scvi-minified",
        os.environ.get("HF_API_TOKEN", None),
        repo_create=False,
        collection_name="test",
    )


@pytest.mark.private
def test_hub_model_pull_from_hf():
    hub_model = HubModel.pull_from_huggingface_hub(repo_name="scvi-tools/test-scvi")
    assert hub_model.model is not None
    assert hub_model.adata is not None

    hub_model = HubModel.pull_from_huggingface_hub(repo_name="scvi-tools/test-scvi-minified")
    assert hub_model.model is not None
    assert hub_model.adata is not None

    hub_model = HubModel.pull_from_huggingface_hub(repo_name="scvi-tools/test-scvi-no-anndata")
    with pytest.raises(FileNotFoundError):
        _ = hub_model.model

    adata = synthetic_iid()
    hub_model.load_model(adata=adata)
    assert hub_model.model is not None
    assert hub_model.adata is None


@pytest.mark.private
def test_hub_model_push_to_s3(save_path: str):
    hub_model = prep_scvi_hub_model(save_path)
    hub_model.push_to_s3("scvi-tools-wis", "tests/hub/test-scvi")

    hub_model = prep_scvi_no_anndata_hub_model(save_path)
    with pytest.raises(ValueError):
        hub_model.push_to_s3("scvi-tools-wis", "tests/hub/test-scvi-no-anndata", push_anndata=True)
    hub_model.push_to_s3("scvi-tools-wis", "tests/hub/test-scvi-no-anndata", push_anndata=False)

    hub_model = prep_scvi_minified_hub_model(save_path)
    hub_model.push_to_s3("scvi-tools-wis", "tests/hub/test-scvi-minified")


@pytest.mark.private
def test_hub_model_pull_from_s3():
    from botocore.exceptions import ClientError

    hub_model = HubModel.pull_from_s3(
        "scvi-tools-wis",
        "tests/hub/test-scvi",
    )
    assert hub_model.model is not None
    assert hub_model.adata is not None

    hub_model = HubModel.pull_from_s3("scvi-tools-wis", "tests/hub/test-scvi-minified")
    assert hub_model.model is not None
    assert hub_model.adata is not None

    with pytest.raises(ClientError):
        hub_model = HubModel.pull_from_s3("scvi-tools-wis", "tests/hub/test-scvi-no-anndata")

    hub_model = HubModel.pull_from_s3(
        "scvi-tools-wis",
        "tests/hub/test-scvi-no-anndata",
        pull_anndata=False,
    )
    with pytest.raises(ValueError):
        _ = hub_model.model

    adata = synthetic_iid()
    hub_model.load_model(adata=adata)
    assert hub_model.model is not None
    assert hub_model.adata is None
