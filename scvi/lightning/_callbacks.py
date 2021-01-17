from pytorch_lightning.callbacks import Callback


class SubSampleLabels(Callback):
    def __init__(self):
        super().__init__()

    def on_epoch_start(self, trainer, pl_module):
        trainer.train_dataloader.resample_labels()
        super().on_epoch_start(trainer, pl_module)
