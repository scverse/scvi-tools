# Training configuration

scvi-tools lets you configure training with structured configuration objects in addition to the
existing `plan_kwargs` / `**trainer_kwargs` pattern.

## When to use which config

- `plan_config`: configures the training plan (optimizer, KL warmup, compile, etc.).
- `trainer_config`: configures the Lightning trainer (early stopping, logging, callbacks, etc.).

Both are optional and can be mixed with kwargs:

- `plan_kwargs` overrides values from `plan_config`.
- `trainer_kwargs` overrides values from `trainer_config`.

## Example (SCVI)

```python
import scvi
from scvi.train import TrainingPlanConfig, TrainerConfig

model = scvi.model.SCVI(adata)

plan_config = TrainingPlanConfig(lr=1e-3, weight_decay=0.0, compile=False)
trainer_config = TrainerConfig(early_stopping=True, early_stopping_patience=10)

model.train(
    max_epochs=100,
    plan_config=plan_config,
    trainer_config=trainer_config,
)
```

## Example (SCANVI)

SCANVI uses a semisupervised training plan, so use
`SemiSupervisedTrainingPlanConfig`:

```python
from scvi.train import SemiSupervisedTrainingPlanConfig

plan_config = SemiSupervisedTrainingPlanConfig(lr=1e-3, compile=False)
model_scanvi.train(max_epochs=20, plan_config=plan_config)
```

## Choosing the right plan config

Use the plan config that matches the training plan behind your model:

- `TrainingPlanConfig` → most unsupervised models (e.g. SCVI).
- `SemiSupervisedTrainingPlanConfig` → SCANVI and other semi‑supervised models.
- `AdversarialTrainingPlanConfig` / `SemiSupervisedAdversarialTrainingPlanConfig` →
  models with adversarial mixing.
- `PyroTrainingPlanConfig` / `LowLevelPyroTrainingPlanConfig` → Pyro‑based models.
- `ClassifierTrainingPlanConfig` → classifier training plans.
- `JaxTrainingPlanConfig` → Jax training plans.

If you’re unsure, you can still use `plan_kwargs` exactly as before.
