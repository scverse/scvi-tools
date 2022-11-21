import logging

from pytorch_lightning.callbacks import ProgressBarBase

from scvi.utils import track

logger = logging.getLogger(__name__)


class ProgressBar(ProgressBarBase):
    """
    Custom progress bar for scvi-tools models.

    Parameters
    ----------
    refresh_rate
        Determines at which rate (in number of epochs) the progress bars get updated.
        Set it to ``0`` to disable the display.
    """

    def __init__(self, refresh_rate: int = 1):
        super().__init__()
        if refresh_rate > 1:
            raise ValueError(
                "scvi-tools progress bar only supports a value of 0 of 1 for `progress_bar_refresh_rate`"
            )
        self._refresh_rate = refresh_rate
        self._enabled = True

    def __getstate__(self):
        # can't pickle the tqdm objects
        state = self.__dict__.copy()
        state["main_progress_bar"] = None
        return state

    @property
    def is_enabled(self) -> bool:  # noqa: D102
        return self._enabled and self.refresh_rate > 0

    @property
    def refresh_rate(self) -> int:  # noqa: D102
        return self._refresh_rate

    @property
    def is_disabled(self) -> bool:  # noqa: D102
        return not self.is_enabled

    def disable(self) -> None:  # noqa: D102
        self._enabled = False

    def enable(self) -> None:  # noqa: D102
        self._enabled = True

    def init_train_tqdm(self, trainer):
        """Override this to customize the tqdm bar for training."""
        bar = track(
            None,
            total=trainer.max_epochs,
            description="Training",
            style="tqdm",
            disable=self.is_disabled,
        )
        return bar

    def on_sanity_check_start(self, trainer, pl_module):  # noqa: D102
        super().on_sanity_check_start(trainer, pl_module)
        logger.info("Running sanity check on val set...")

    def on_train_start(self, trainer, pl_module):  # noqa: D102
        super().on_train_start(trainer, pl_module)
        self.main_progress_bar = self.init_train_tqdm(trainer)

    def on_train_epoch_start(self, trainer, pl_module):  # noqa: D102
        super().on_train_epoch_start(trainer, pl_module)
        if self._should_update(self.trainer.current_epoch, self.trainer.max_epochs):
            epoch = trainer.current_epoch + 1
            self.main_progress_bar.set_description(
                f"Epoch {epoch}/{trainer.max_epochs}"
            )

    def _should_update(self, current, total):
        return self.is_enabled and (
            current % self.refresh_rate == 0 or current == total
        )

    def on_train_epoch_end(self, trainer, pl_module):  # noqa: D102
        super().on_train_epoch_end(trainer, pl_module)
        if self._should_update(self.trainer.current_epoch, self.trainer.max_epochs):
            self.main_progress_bar.update()
            self.main_progress_bar.set_postfix(self.get_metrics(trainer, pl_module))

    def on_train_end(self, trainer, pl_module):  # noqa: D102
        super().on_train_end(trainer, pl_module)
        if self.is_enabled:
            self.main_progress_bar.close()


def convert_inf(x):
    """The tqdm doesn't support inf values. We have to convert it to None."""
    if x == float("inf"):
        return None
    return x
