import logging

from pytorch_lightning.callbacks import ProgressBarBase

from scvi._utils import track

logger = logging.getLogger(__name__)


class ProgressBar(ProgressBarBase):
    """Custom progress bar for scvi-tools models."""

    def __init__(self):
        super().__init__()
        self._enabled = True

    def __getstate__(self):
        # can't pickle the tqdm objects
        state = self.__dict__.copy()
        state["main_progress_bar"] = None
        return state

    @property
    def is_enabled(self) -> bool:
        return self._enabled

    @property
    def is_disabled(self) -> bool:
        return not self.is_enabled

    def disable(self) -> None:
        self._enabled = False

    def enable(self) -> None:
        self._enabled = True

    def init_train_tqdm(self, trainer):
        """Override this to customize the tqdm bar for training."""
        bar = track(
            None,
            total=trainer.max_epochs,
            description="Training",
            style="tqdm",
            initial=self.train_batch_idx,
            disable=self.is_disabled,
        )
        return bar

    def on_sanity_check_start(self, trainer, pl_module):
        super().on_sanity_check_start(trainer, pl_module)
        logger.info("Running sanity check on val set...")

    def on_train_start(self, trainer, pl_module):
        super().on_train_start(trainer, pl_module)
        self.main_progress_bar = self.init_train_tqdm(trainer)

    def on_epoch_start(self, trainer, pl_module):
        super().on_epoch_start(trainer, pl_module)
        epoch = trainer.current_epoch + 1
        self.main_progress_bar.set_description(f"Epoch {epoch}/{trainer.max_epochs}")

    def on_train_epoch_end(self, trainer, pl_module, outputs):
        super().on_train_epoch_end(trainer, pl_module, outputs)
        if self.is_enabled:
            self.main_progress_bar.update()
            self.main_progress_bar.set_postfix(trainer.progress_bar_dict)

    def on_train_end(self, trainer, pl_module):
        super().on_train_end(trainer, pl_module)
        self.main_progress_bar.close()


def convert_inf(x):
    """The tqdm doesn't support inf values. We have to convert it to None."""
    if x == float("inf"):
        return None
    return x
