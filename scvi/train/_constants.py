from typing import NamedTuple


class _METRIC_KEYS_NT(NamedTuple):
    TRAINING_KEY: str = "training"
    VALIDATION_KEY: str = "validation"
    # classification
    ACCURACY_KEY: str = "accuracy"
    F1_SCORE_KEY: str = "f1_score"
    CALIBRATION_ERROR_KEY: str = "calibration_error"
    AUROC_KEY: str = "auroc"
    CLASSIFICATION_LOSS_KEY: str = "classification_loss"
    TRUE_LABELS_KEY: str = "true_labels"
    LOGITS_KEY: str = "logits"


METRIC_KEYS = _METRIC_KEYS_NT()
