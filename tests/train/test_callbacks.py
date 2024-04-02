import os

import scvi
from scvi.train._callbacks import SaveCheckpoint


def test_savecheckpoint(save_path: str):
    def check_checkpoint_logging(model, adata):
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

    def test_model_cls(model_cls, adata):
        scvi.settings.logging_dir = os.path.join(save_path, model_cls.__name__)

        # enable_checkpointing=True, default monitor
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
            callbacks=[SaveCheckpoint(monitor="elbo_validation")],
        )
        check_checkpoint_logging(model, adata)

    old_logging_dir = scvi.settings.logging_dir
    adata = scvi.data.synthetic_iid()

    scvi.model.SCVI.setup_anndata(adata)
    test_model_cls(scvi.model.SCVI, adata)

    scvi.model.SCANVI.setup_anndata(adata, "labels", "label_0")
    test_model_cls(scvi.model.SCANVI, adata)

    scvi.settings.logging_dir = old_logging_dir
