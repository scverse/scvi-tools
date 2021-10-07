from typing import Optional, Tuple

import pytorch_lightning as pl
import torch
from pytorch_lightning.callbacks.early_stopping import EarlyStopping


class LoudEarlyStopping(EarlyStopping):
    """
    Wrapper of Pytorch Lightning EarlyStopping callback that prints the reason for stopping on teardown.

    When the early stopping condition is met, the reason is saved to the callback instance,
    then printed on teardown. By printing on teardown, we do not interfere with the progress bar callback.
    """

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.early_stopping_reason = None

    def _evaluate_stopping_criteria(self, current: torch.Tensor) -> Tuple[bool, str]:
        should_stop, reason = super()._evaluate_stopping_criteria(current)
        if should_stop:
            self.early_stopping_reason = reason
        return should_stop, reason

    def teardown(
        self,
        _trainer: pl.Trainer,
        _pl_module: pl.LightningModule,
        stage: Optional[str] = None,
    ) -> None:
        if self.early_stopping_reason is not None:
            print(self.early_stopping_reason)
