"""Callbacks for DiagVI training."""

from collections import defaultdict

from lightning.pytorch.callbacks import Callback
from scvi import REGISTRY_KEYS


class CellTypeProportionsPrinter(Callback):
    """Print cell type proportions per batch and modality during validation."""

    def __init__(self):
        super().__init__()
        self.label_names = {}  # {modality_name: {label_id: label_name}}

    def on_train_start(self, trainer, pl_module):
        """Extract label names from the training module at training start."""
        # Try to access adata_managers from the module
        if hasattr(pl_module, "adata_managers"):
            adata_managers = pl_module.adata_managers
        elif hasattr(pl_module, "module") and hasattr(pl_module.module, "adata_managers"):
            adata_managers = pl_module.module.adata_managers
        else:
            # Cannot access adata_managers, will use label_0, label_1, etc.
            return

        # Extract label names from each modality's registry
        for modality_name, adata_manager in adata_managers.items():
            try:
                labels_registry = adata_manager.get_state_registry(REGISTRY_KEYS.LABELS_KEY)
                # The categorical_mapping should be a dict or list of label names
                mapping = labels_registry.get("categorical_mapping", {})
                self.label_names[modality_name] = mapping
            except Exception:
                self.label_names[modality_name] = {}

    def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
        """Print cell type proportions at end of first 4 training steps per epoch."""
        # Only print for the first 4 batches of each epoch
        if batch_idx >= 4:
            return

        epoch = trainer.current_epoch
        step = trainer.global_step

        proportions_by_modality = defaultdict(
            lambda: defaultdict(lambda: defaultdict(int))
        )
        batch_counts = defaultdict(lambda: defaultdict(int))

        # Process each modality in the batch
        for modality_name, tensors in batch.items():
            batch_indices = tensors[REGISTRY_KEYS.BATCH_KEY].cpu().numpy()
            label_indices = tensors[REGISTRY_KEYS.LABELS_KEY].cpu().numpy()

            # Count cell types per batch
            for batch_id, label_id in zip(batch_indices, label_indices):
                batch_id = int(batch_id)
                label_id = int(label_id)
                proportions_by_modality[modality_name][batch_id][label_id] += 1
                batch_counts[modality_name][batch_id] += 1

        # Print proportions
        print(f"\nTraining Step {step} (Epoch {epoch}, Batch {batch_idx}):")
        for modality_name in sorted(proportions_by_modality.keys()):
            print(f"  {modality_name}:", end="")
            for batch_id in sorted(proportions_by_modality[modality_name].keys()):
                total = batch_counts[modality_name][batch_id]
                props = []
                for label_id in sorted(
                    proportions_by_modality[modality_name][batch_id].keys()
                ):
                    count = proportions_by_modality[modality_name][batch_id][label_id]
                    prop = count / total

                    # Try to get label name
                    if (
                        modality_name in self.label_names
                        and label_id in self.label_names[modality_name]
                    ):
                        label_name = self.label_names[modality_name][label_id]
                    else:
                        label_name = f"label_{label_id}"

                    props.append(f"{label_name}={prop:.3f}")

                print(f" batch_{batch_id}: {', '.join(props)}", end="")
            print()  # Newline

    def on_validation_epoch_end(self, trainer, pl_module):
        """Print cell type proportions at end of validation epoch."""
        epoch = trainer.current_epoch

        if not trainer.val_dataloaders:
            return

        # Collect proportions across all validation batches
        proportions_by_modality = defaultdict(
            lambda: defaultdict(lambda: defaultdict(int))
        )
        batch_counts = defaultdict(lambda: defaultdict(int))

        # Iterate through validation dataloader (TrainDL yields dict of modality batches)
        for batch in trainer.val_dataloaders:
            for modality_name, tensors in batch.items():
                batch_indices = tensors[REGISTRY_KEYS.BATCH_KEY].cpu().numpy()
                label_indices = tensors[REGISTRY_KEYS.LABELS_KEY].cpu().numpy()

                # Count cell types per batch
                for batch_id, label_id in zip(batch_indices, label_indices):
                    batch_id = int(batch_id)
                    label_id = int(label_id)
                    proportions_by_modality[modality_name][batch_id][label_id] += 1
                    batch_counts[modality_name][batch_id] += 1

        # Convert counts to proportions and print
        print(f"\n=== Epoch {epoch} Validation Cell Type Proportions ===")
        for modality_name in sorted(proportions_by_modality.keys()):
            print(f"\n  {modality_name}:")
            for batch_id in sorted(proportions_by_modality[modality_name].keys()):
                total = batch_counts[modality_name][batch_id]
                props = []
                for label_id in sorted(
                    proportions_by_modality[modality_name][batch_id].keys()
                ):
                    count = proportions_by_modality[modality_name][batch_id][label_id]
                    prop = count / total

                    # Try to get label name
                    if (
                        modality_name in self.label_names
                        and label_id in self.label_names[modality_name]
                    ):
                        label_name = self.label_names[modality_name][label_id]
                    else:
                        label_name = f"label_{label_id}"

                    props.append(f"{label_name}={prop:.3f}")

                print(f"    batch_{batch_id}: {', '.join(props)}")
