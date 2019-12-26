from .posterior import Posterior
from .trainer import Trainer
from .inference import UnsupervisedTrainer, AdapterTrainer
from .annotation import (
    JointSemiSupervisedTrainer,
    SemiSupervisedTrainer,
    AlternateSemiSupervisedTrainer,
    ClassifierTrainer,
)
from .jvae_trainer import JVAETrainer
from .total_inference import TotalPosterior, TotalTrainer
import torch
import numpy as np

__all__ = [
    "Trainer",
    "Posterior",
    "UnsupervisedTrainer",
    "AdapterTrainer",
    "JointSemiSupervisedTrainer",
    "SemiSupervisedTrainer",
    "AlternateSemiSupervisedTrainer",
    "ClassifierTrainer",
    "JVAETrainer",
    "TotalPosterior",
    "TotalTrainer",
]


def set_seed(seed: int = 0):
    torch.manual_seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    np.random.seed(seed)
