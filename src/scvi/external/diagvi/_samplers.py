"""Samplers for DIAGVI model with cell type stratified sampling."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
from torch.utils.data import Sampler

if TYPE_CHECKING:
    from scvi.data import AnnDataManager

logger = logging.getLogger(__name__)


class StratifiedLabelSampler(Sampler):
    """Sampler that yields batches with stratified cell type composition.

    This sampler ensures that each batch contains cells with similar cell type
    proportions across modalities. It samples labels first, then independently
    samples cells from each modality based on those labels.

    Parameters
    ----------
    labels_per_modality
        Dictionary mapping modality names to arrays of label indices for each cell.
        Labels should be integer-encoded (0, 1, 2, ...).
    indices_per_modality
        Dictionary mapping modality names to arrays of cell indices (into the original
        AnnData) that are available for sampling (e.g., training indices).
    batch_size
        Number of cells to sample per batch per modality.
    drop_last
        If True, drop the last incomplete batch.
    shuffle
        If True, shuffle cells within each label group at the start of each epoch.
    """

    def __init__(
        self,
        labels_per_modality: dict[str, np.ndarray],
        indices_per_modality: dict[str, np.ndarray],
        batch_size: int = 128,
        drop_last: bool = False,
        shuffle: bool = True,
    ):
        self.modality_names = list(labels_per_modality.keys())
        if len(self.modality_names) != 2:
            raise ValueError("StratifiedLabelSampler requires exactly 2 modalities.")

        self.batch_size = batch_size
        self.drop_last = drop_last
        self.shuffle = shuffle

        # Build label-to-cell-indices mapping for each modality
        # Only include cells that are in the provided indices (e.g., training set)
        self.label_to_cells: dict[str, dict[int, list[int]]] = {}
        self.all_labels: set[int] = set()

        for mod_name in self.modality_names:
            labels = labels_per_modality[mod_name]
            indices = indices_per_modality[mod_name]

            self.label_to_cells[mod_name] = {}
            for idx in indices:
                label = labels[idx]
                if label not in self.label_to_cells[mod_name]:
                    self.label_to_cells[mod_name][label] = []
                self.label_to_cells[mod_name][label].append(idx)
                self.all_labels.add(label)

        self.all_labels = sorted(self.all_labels)

        # Compute label counts (use minimum across modalities for each label)
        # This ensures we have enough cells from both modalities for each label
        self.label_counts = {}
        for label in self.all_labels:
            counts = []
            for mod_name in self.modality_names:
                if label in self.label_to_cells[mod_name]:
                    counts.append(len(self.label_to_cells[mod_name][label]))
                else:
                    counts.append(0)
            self.label_counts[label] = min(counts)

        # Total number of cells we can sample per epoch (limited by minimum across modalities)
        self.total_cells = sum(self.label_counts.values())

        if self.total_cells == 0:
            raise ValueError(
                "No cells available for stratified sampling. "
                "Check that labels match across modalities."
            )

        # Compute label proportions based on total available cells
        self.label_proportions = {
            label: count / self.total_cells for label, count in self.label_counts.items()
        }

        logger.info(
            f"StratifiedLabelSampler initialized with {len(self.all_labels)} labels, "
            f"{self.total_cells} cells per modality per epoch."
        )

    def __len__(self) -> int:
        """Return number of batches per epoch."""
        if self.drop_last:
            return self.total_cells // self.batch_size
        return (self.total_cells + self.batch_size - 1) // self.batch_size

    def __iter__(self):
        """Yield batches of cell indices for each modality.

        Each iteration yields a dictionary mapping modality names to lists of
        cell indices for that batch.
        """
        # Reset and shuffle cell pools for each label and modality
        cell_pools: dict[str, dict[int, list[int]]] = {}
        for mod_name in self.modality_names:
            cell_pools[mod_name] = {}
            for label in self.all_labels:
                if label in self.label_to_cells[mod_name]:
                    cells = self.label_to_cells[mod_name][label].copy()
                    if self.shuffle:
                        np.random.shuffle(cells)
                    # Only take up to label_counts[label] cells to ensure balance
                    cell_pools[mod_name][label] = cells[: self.label_counts[label]]
                else:
                    cell_pools[mod_name][label] = []

        # Track current position in each label pool
        pool_positions = {label: 0 for label in self.all_labels}

        # Create label sequence proportional to label counts
        label_sequence = []
        for label, count in self.label_counts.items():
            label_sequence.extend([label] * count)
        if self.shuffle:
            np.random.shuffle(label_sequence)

        # Generate batches
        n_batches = len(self)
        label_idx = 0

        for batch_idx in range(n_batches):
            batch_indices = {mod_name: [] for mod_name in self.modality_names}

            # Determine actual batch size (may be smaller for last batch)
            if batch_idx == n_batches - 1 and not self.drop_last:
                current_batch_size = min(
                    self.batch_size, len(label_sequence) - label_idx
                )
            else:
                current_batch_size = self.batch_size

            # Sample cells for this batch
            for _ in range(current_batch_size):
                if label_idx >= len(label_sequence):
                    break

                label = label_sequence[label_idx]
                label_idx += 1

                # Sample one cell from each modality for this label
                for mod_name in self.modality_names:
                    pool = cell_pools[mod_name][label]
                    pos = pool_positions[label]
                    if pos < len(pool):
                        batch_indices[mod_name].append(pool[pos])

                pool_positions[label] += 1

            # Only yield if we have cells in the batch
            if all(len(indices) > 0 for indices in batch_indices.values()):
                yield batch_indices


class StratifiedBatchSampler(Sampler):
    """Batch sampler wrapper for a single modality using stratified sampling.

    This is used internally by StratifiedTrainDL to create individual
    data loaders that are coordinated through shared label sampling.

    Parameters
    ----------
    indices_iterator
        An iterator that yields lists of indices for each batch.
    """

    def __init__(self, indices_iterator):
        self._indices_iterator = indices_iterator
        self._cached_batches = None

    def set_batches(self, batches: list[list[int]]):
        """Set the pre-computed batches for this sampler."""
        self._cached_batches = batches

    def __iter__(self):
        if self._cached_batches is None:
            raise RuntimeError("Batches not set. Call set_batches() first.")
        for batch in self._cached_batches:
            yield batch

    def __len__(self):
        if self._cached_batches is None:
            return 0
        return len(self._cached_batches)
