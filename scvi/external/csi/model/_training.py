from scvi.external.csi.train import TrainingPlanCustom
from scvi.model.base import UnsupervisedTrainingMixin


class TrainingCustom(UnsupervisedTrainingMixin):
    """Train method with custom TrainingPlan."""

    # TODO could make custom Trainer (in a custom TrainRunner) to have in init params for early stopping
    #  "loss" rather than "elbo" components in available param specifications - for now just use
    #  a loss that is against the param specification

    # TODO run and log val before training - already tried some solutions by calling trainer.validate before
    #  fit and num_sanity_val_steps (designed not to log)
    _training_plan_cls = TrainingPlanCustom
