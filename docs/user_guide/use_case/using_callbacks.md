# Train SCVI model with callbacks

In PyTorch Lightning, callbacks are special functions that let you execute custom actions during the training process, like saving checkpoints, adjusting learning rates, or early stopping based on performance. To use callbacks, you first create a callback class that defines what should happen at specific points in training (e.g., at the end of an epoch or batch). Then, when setting up your Trainer, you simply pass a list of these callbacks. For example, if you want to save the best model during training, you can use ModelCheckpoint to automatically store the model when it achieves the best validation score.

In scvi-tools we provide two custom callbacks based on that:
1. Early Stopping callback: {class}`scvi.train._callbacks.LoudEarlyStopping`

With this callback, training will stop when its monitored metric is reaching an optimal point for and did not improve for several epochs.
There are several common parameters that can be controlled when running this callback:
- early_stopping: A boolean on whether to activate the early stopping callback.
- monitor: There are several metric keys we automatically record while training:
 "elbo_validation",
 "reconstruction_loss_validation",
 "kl_local_validation",
 "elbo_train",
 "reconstruction_loss_train",
 "kl_local_train",
 "validation_classification_loss",
 "validation_accuracy",
 "validation_f1_score",
 "validation_calibration_error". Those, per model, can be selected as a monitored metric.
- min_delta: minimum change in the monitored quantity to qualify as an improvement, i.e. an absolute
    change of less than or equal to `min_delta`, will count as no improvement.
- patience: number of checks with no improvement
    after which training will be stopped. Under the default configuration, one check happens after
    every training epoch. However, the frequency of validation can be modified by setting various parameters on
    the ``Trainer``, for example ``check_val_every_n_epoch`` and ``val_check_interval``. by default, if not set, ``check_val_every_n_epoch`` will be 1, thus adding computation overhead to the training step.
- mode: one of ``'min'``, ``'max'``. In ``'min'`` mode, training will stop when the quantity
    monitored has stopped decreasing and in ``'max'`` mode it will stop when the quantity
    monitored has stopped increasing.

Example of usage:
```python
early_stopping_kwargs = {
    "early_stopping": True,
    "early_stopping_monitor": "elbo_validation",
    "early_stopping_patience": 10,
    "early_stopping_min_delta": 0.0,
    "check_val_every_n_epoch": 1,
}
model.train(..., **early_stopping_kwargs)
```

Several models will be trained with early stopping by default, e.g {class}`~scvi.model.TOTALVI`, {class}`~scvi.model.MULTIVI` and others. Consider disabling it when there is no need.

2. Model Checkpoint Callback: {class}`scvi.train._callbacks.SaveCheckpoint`

Saves model checkpoints based on a monitored metric. The best model saved and best model score based on ``monitor`` can be accessed post-training
    with the ``best_model_path`` and ``best_model_score`` attributes, respectively. Starting in scvi-tools 1.3.0, we added the on_exception option to this callback, which in case of model error exceptions during training, resulting from Nan's in loss or gradients, will save the best model ("best" in terms of what is the monitored metric). It does not gracefully shutdown, but it is the user responsibility to load this model and continue the analysis, e.g., user can take the optimal model that was saved, or continue training it with perhaps different training parameters, to prevent the model from failing to train. Note this option will add some overhead to the training time.

It can be used by adding the following parameter to the train function:
```python
model.train(
    ...,
    callbacks=[
        SaveCheckpoint(
            monitor="elbo_validation", load_best_on_end=True, check_nan_gradients=True
        )
    ],
)
```
