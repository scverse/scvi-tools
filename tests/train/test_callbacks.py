import os
from shutil import rmtree

import pytest
import torch
from lightning.pytorch.callbacks import ModelCheckpoint

import scvi
from scvi.data import synthetic_iid
from scvi.model import MULTIVI, SCANVI, SCVI, TOTALVI
from scvi.train._callbacks import SaveCheckpoint, ScibCallback


@pytest.mark.parametrize("load_best_on_end", [False, True])
def test_savecheckpoint(load_best_on_end: bool, save_path: str):
    from anndata import AnnData

    from scvi.model.base import BaseModelClass

    def check_checkpoint_logging(model: BaseModelClass, adata: AnnData):
        assert any(isinstance(c, SaveCheckpoint) for c in model.trainer.callbacks)

        callback = [c for c in model.trainer.callbacks if isinstance(c, SaveCheckpoint)]
        assert len(callback) == 1

        callback = callback[0]
        assert callback.best_model_path is not None
        assert callback.best_model_score is not None
        assert os.path.exists(callback.best_model_path)

        log_dirs = os.listdir(scvi.settings.logging_dir)
        assert len(log_dirs) >= 1

        checkpoints_dir = os.path.join(scvi.settings.logging_dir, log_dirs[0])
        checkpoint_dirs = os.listdir(checkpoints_dir)
        assert len(checkpoint_dirs) >= 1

        checkpoint_dir = os.path.join(checkpoints_dir, checkpoint_dirs[0])
        checkpoint = model.__class__.load(checkpoint_dir, adata=adata)
        assert checkpoint.is_trained_

        if load_best_on_end:
            best_model = model.__class__.load(callback.best_model_path, adata=adata)
            assert best_model.is_trained_

            # we assume that the current model state is equal to the best one (not true always...)
            current_state_dict = model.module.state_dict()
            best_state_dict = best_model.module.state_dict()
            assert len(current_state_dict) == len(best_state_dict)
            for k, v in current_state_dict.items():
                assert torch.equal(v, best_state_dict[k])
                assert v.device == best_state_dict[k].device

    def test_model_cls(model_cls, adata: AnnData):
        scvi.settings.logging_dir = os.path.join(save_path, model_cls.__name__)

        # Clean the folder so this test run is isolated/deterministic.
        if os.path.exists(scvi.settings.logging_dir):
            rmtree(scvi.settings.logging_dir)
        os.makedirs(scvi.settings.logging_dir, exist_ok=True)

        # enable_checkpointing=True, default monitor ("validation_loss" for this test)
        model = model_cls(adata)
        model.train(max_epochs=5, enable_checkpointing=True)
        check_checkpoint_logging(model, adata)

        # enable_checkpointing=True, custom monitor
        model = model_cls(adata)
        model.train(
            max_epochs=5,
            enable_checkpointing=True,
            checkpointing_monitor="elbo_validation",
        )
        check_checkpoint_logging(model, adata)

        # manual callback
        model = model_cls(adata)
        model.train(
            max_epochs=5,
            callbacks=[
                SaveCheckpoint(monitor="elbo_validation", load_best_on_end=load_best_on_end)
            ],
        )
        check_checkpoint_logging(model, adata)

    old_logging_dir = scvi.settings.logging_dir
    adata = scvi.data.synthetic_iid()

    SCVI.setup_anndata(
        adata,
        batch_key="batch",
    )
    test_model_cls(SCVI, adata)

    SCANVI.setup_anndata(
        adata, batch_key="batch", labels_key="labels", unlabeled_category="label_0"
    )
    test_model_cls(SCANVI, adata)

    load_best_on_end = False  # hard coded for TOTALVI
    TOTALVI.setup_anndata(
        adata,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    test_model_cls(TOTALVI, adata)

    mdata = synthetic_iid(return_mudata=True)
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities={"rna_layer": "rna", "protein_layer": "protein_expression"},
    )
    test_model_cls(TOTALVI, adata)

    scvi.settings.logging_dir = old_logging_dir


@pytest.mark.parametrize("model_cls", [SCVI, TOTALVI])
def test_exception_callback(model_cls):
    torch.set_float32_matmul_precision("high")
    scvi.settings.seed = 0

    # we still need to find a proper way to simulate an adata that fail quickly during training
    adata = synthetic_iid()

    if model_cls == SCVI:
        model_cls.setup_anndata(adata, batch_key="batch")
    if model_cls == TOTALVI:
        model_cls.setup_anndata(
            adata,
            batch_key="batch",
            protein_expression_obsm_key="protein_expression",
            protein_names_uns_key="protein_names",
        )
    model = model_cls(adata)
    model.train(max_epochs=5)

    ckpt_cb = SaveCheckpoint(
        dirpath="checkpoints/",
        monitor="elbo_validation",
        mode="min",
        save_top_k=1,
        save_last=True,
        load_best_on_end=True,
        check_nan_gradients=True,
    )

    model.train(
        max_epochs=5,
        check_val_every_n_epoch=1,
        callbacks=[ckpt_cb],
        enable_checkpointing=True,
    )


@pytest.mark.parametrize("metric", ["Total", "Bio conservation", "iLISI"])
@pytest.mark.parametrize("model_cls", [SCVI, SCANVI, TOTALVI])
def test_scib_callback_adata(model_cls, metric: str):
    adata = synthetic_iid()
    if model_cls == TOTALVI:
        model_cls.setup_anndata(
            adata,
            batch_key="batch",
            protein_expression_obsm_key="protein_expression",
            protein_names_uns_key="protein_names",
        )
    else:
        if model_cls == SCANVI:
            model_cls.setup_anndata(
                adata,
                labels_key="labels",
                unlabeled_category="unknown",
                batch_key="batch",
            )
        else:
            model_cls.setup_anndata(
                adata,
                labels_key="labels",
                batch_key="batch",
            )
    model = model_cls(adata)
    model.train(
        1,
        train_size=0.5,
        callbacks=[ScibCallback()],
    )


@pytest.mark.parametrize("metric", ["Total", "Bio conservation", "iLISI"])
@pytest.mark.parametrize("model_cls", [MULTIVI, TOTALVI])
def test_scib_callback_mudata(model_cls, metric: str):
    mdata = synthetic_iid(return_mudata=True)
    if model_cls == MULTIVI:
        model_cls.setup_mudata(
            mdata,
            batch_key="batch",
            modalities={
                "rna_layer": "rna",
                "atac_layer": "accessibility",
                "protein_layer": "protein_expression",
            },
        )
    else:
        model_cls.setup_mudata(
            mdata,
            batch_key="batch",
            modalities={"rna_layer": "rna", "protein_layer": "protein_expression"},
        )
    model = model_cls(mdata)
    model.train(
        1,
        train_size=0.5,
        callbacks=[ScibCallback()],
    )


@pytest.mark.parametrize("model_cls", [SCVI, TOTALVI])
def test_lightning_checkpoint(model_cls):
    adata = synthetic_iid()

    if model_cls == SCVI:
        model_cls.setup_anndata(adata, batch_key="batch")
    if model_cls == TOTALVI:
        model_cls.setup_anndata(
            adata,
            batch_key="batch",
            protein_expression_obsm_key="protein_expression",
            protein_names_uns_key="protein_names",
        )
    model = model_cls(adata)
    model.train(max_epochs=5)

    ckpt_cb = ModelCheckpoint(
        dirpath="checkpoints/",
        monitor="elbo_validation",
        mode="min",
        save_top_k=1,  # keep the best model
        save_last=True,  # also keep a "last.ckpt"
    )
    model.train(
        max_epochs=1, check_val_every_n_epoch=1, callbacks=[ckpt_cb], enable_checkpointing=True
    )

    print("Best ckpt:", ckpt_cb.best_model_path)

    ckpt = torch.load(ckpt_cb.best_model_path)

    assert "optimizer_states" in ckpt
    assert "lr_schedulers" in ckpt
    assert "state_dict" in ckpt

    model.train(
        max_epochs=1,
        check_val_every_n_epoch=1,
        callbacks=[ckpt_cb],
        enable_checkpointing=True,
        ckpt_path=ckpt_cb.best_model_path,
    )

    ckpt = torch.load(ckpt_cb.best_model_path)

    assert "optimizer_states" in ckpt
    assert "lr_schedulers" in ckpt
    assert "state_dict" in ckpt
