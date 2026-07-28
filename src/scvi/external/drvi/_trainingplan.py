from __future__ import annotations

import torch

from scvi.train import TrainingPlan
from scvi.train._trainingplans import _compute_kl_weight


class DRVITrainingPlan(TrainingPlan):
    """:class:`~scvi.train.TrainingPlan` that anneals the KL weight over the whole training run.

    DRVI's disentanglement objective benefits from spreading the KL warmup across the entire run,
    especially when the number of epochs is small (compare to the default of ``400``).
    With the default ``n_epochs_kl_warmup="auto"``, the warmup length is derived from the
    trainer's ``max_epochs`` at train time, so it matches the requested number of epochs.
    Passing an explicit ``n_epochs_kl_warmup`` (an ``int`` or ``None``) restores the standard
    behavior.
    """

    def __init__(self, *args, n_epochs_kl_warmup: int | str | None = "auto", **kwargs):
        self._auto_kl_warmup = n_epochs_kl_warmup == "auto"
        super().__init__(
            *args,
            n_epochs_kl_warmup=None if self._auto_kl_warmup else n_epochs_kl_warmup,
            **kwargs,
        )

    @property
    def kl_weight(self):
        """Scaling factor on KL divergence, warming up over ``max_epochs`` when ``"auto"``."""
        n_epochs_kl_warmup = self.n_epochs_kl_warmup
        if self._auto_kl_warmup:
            try:
                n_epochs_kl_warmup = self.trainer.max_epochs
            except RuntimeError:
                # Not yet attached to a trainer (``kl_weight`` is read once during __init__).
                # The exact warmup length is irrelevant here since ``current_epoch`` is 0; any
                # positive value yields ``min_kl_weight``.
                n_epochs_kl_warmup = 1
        klw = _compute_kl_weight(
            self.current_epoch,
            self.global_step,
            n_epochs_kl_warmup,
            self.n_steps_kl_warmup,
            self.max_kl_weight,
            self.min_kl_weight,
        )
        return torch.tensor(klw).to(self.device)
