"""
Specialized data loading and processing components for the SENADVAE

Implementing dynamic control-perturbation matching and efficient batch processing for
single-cell perturbation analysis.
"""

import logging
from collections.abc import Sequence

import numpy as np
import pandas as pd
import torch
from lightning.pytorch.callbacks import Callback

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.dataloaders import AnnDataLoader, DataSplitter

logger = logging.getLogger(__name__)


class SENADataLoader(AnnDataLoader):
    """
    Custom dataloader for SENADVAE with dynamic preprocessing during training.

    This dataloader processes single-cell gene expression data with perturbations, handling
    the complex task of matching control cells to perturbed cells and converting perturbation
    annotations into a format suitable for neural network training. It inherits from AnnDataLoader
    to maintain compatibility with the scvi-tools framework.

    Parameters
    ----------
    adata_manager : AnnDataManager
        AnnDataManager instance containing the original adata.
    perturbation_key : str
        Key in adata.obs containing perturbation information.
    shuffle : bool, default True
        Whether to shuffle the data.
    batch_size : int, default 128
        Batch size for training.
    mode : str, default "training"
        Mode for dataloader: "training" or "prediction".
    indices : Sequence[int], optional
        Specific indices for prediction mode.
    prediction_indices : Sequence[int], optional
        Alternative parameter name for indices (used by SENADataSplitter).
    **kwargs
        Additional arguments passed to DataLoader.

    Attributes
    ----------
    adata_manager : AnnDataManager
        AnnDataManager instance containing the original adata.
    adata : AnnData
        The AnnData object from adata_manager.
    perturbation_key : str
        Field name where the perturbation gene names are stored.
    mode : str
        Mode for dataloader: "training" or "prediction".
    prediction_indices : Sequence[int] or None
        Specific indices for prediction mode.
    matched_ctrl_X : array_like
        Gene expression data from matched control cells.
    perturbed_X : array_like
        Gene expression data from perturbed cells.
    perturbation_strings : array_like
        Perturbation annotation strings.
    intervention_matrix : np.ndarray
        Numerical intervention matrix encoding perturbation states.
    intervention_genes : list
        List of all intervention genes.
    gene_to_intervention_idx : dict
        Mapping from gene names to indices.
    n_intervention_genes : int
        Number of unique intervention genes.

    Notes
    -----
    Key Features:
    1. Control-Perturbation Pair Management:
       - Identifies and separates control cells from perturbed cells
       - Creates matched pairs by randomly sampling controls for each perturbation
       - Implements dynamic control reassignment between epochs to improve generalization

    2. Perturbation Processing:
       - Converts gene perturbation annotations into numerical matrices
       - Supports both single-gene and multi-gene perturbations
       - Creates efficient sparse representation of intervention states

    3. Data Format:
       Processes and returns data in a standardized format for training:
       - X: Gene expression data from matched control cells
       - labels: Gene expression data from perturbed cells
       - cat_covs: Intervention matrix encoding perturbation states:
         * -1: Control state
         * 0: Non-perturbed gene
         * 1: Perturbed gene

    The dataloader implements efficient batch processing and ensures consistent
    data preprocessing across training, validation, and inference phases.
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        perturbation_key: str,
        shuffle: bool = True,
        batch_size: int = 128,
        mode: str = "training",  # "training" or "prediction"
        indices: Sequence[int] | None = None,  # For prediction mode
        prediction_indices: Sequence[int] | None = None,  # Alternative name for indices
        **kwargs,
    ):
        """
        Initialize SENA dataloader with dynamic preprocessing capabilities.

        Sets up the complete data processing pipeline for single-cell perturbation analysis,
        including control-perturbation matching, intervention matrix creation, and batch
        processing compatible with scvi-tools framework.

        Parameters
        ----------
        adata_manager : AnnDataManager
            AnnDataManager instance containing the original adata.
        perturbation_key : str
            Key in adata.obs containing perturbation information.
        shuffle : bool, default True
            Whether to shuffle the data.
        batch_size : int, default 128
            Batch size for training.
        mode : str, default "training"
            Mode for dataloader: "training" or "prediction".
        indices : Sequence[int], optional
            Specific indices for prediction mode.
        prediction_indices : Sequence[int], optional
            Alternative parameter name for indices (used by SENADataSplitter).
        **kwargs
            Additional arguments passed to DataLoader.
        """
        self.adata_manager = adata_manager
        self.adata = adata_manager.adata
        self.perturbation_key = (
            perturbation_key  # field name where the perturbation gene names are stored
        )
        self.mode = mode  # "training" or "prediction"

        # Handle both parameter names for indices
        self.prediction_indices = prediction_indices if prediction_indices is not None else indices

        # Store these parameters explicitly since we need them later
        self._shuffle = shuffle
        self._batch_size = batch_size

        # Preprocess the data to create the matched pairs (same logic as setup_anndata)
        self._create_matched_dataset()

        # Create indices for the processed dataset
        self.dataset_indices = np.arange(self.matched_ctrl_X.shape[0])

        # Filter out custom parameters before passing to parent
        parent_kwargs = {
            k: v for k, v in kwargs.items() if k not in ["prediction_indices", "mode"]
        }

        # Initialize parent with our custom dataset
        super().__init__(
            adata_manager=adata_manager,
            indices=self.dataset_indices,
            batch_size=batch_size,
            shuffle=shuffle,
            **parent_kwargs,
        )

        # Store these parameters for our custom iterator (don't override PyTorch attributes)
        self._custom_shuffle = shuffle
        self._custom_batch_size = batch_size

    def _create_matched_dataset(self):
        """
        Create matched control-perturbation pairs for training or prediction.

        This method processes the AnnData object to:
        1. Identify control and perturbed cells
        2. Filter based on mode:
           - Training mode: Only single perturbations
           - Prediction mode: Use specified indices (can include double perturbations)
        3. Create matched pairs by random sampling from control cells
        4. Prepare the data structure for batch processing

        Raises
        ------
        ValueError
            If no perturbed cells or control cells are found in AnnData object.
        ValueError
            If no single perturbation cells found in training mode.
        ValueError
            If unknown mode is specified.

        Notes
        -----
        The method implements different filtering strategies based on mode:
        - Training mode: Restricts to single gene perturbations only for stable learning
        - Prediction mode: Allows any perturbation complexity including combinations
        """
        logger.info(
            f"Creating matched control-perturbation dataset in dataloader (mode: {self.mode})."
        )

        # Separate control and perturbed cells
        is_control = self.adata.obs[self.perturbation_key] == ""
        ctrl_indices = np.where(is_control)[0]
        ptb_indices = np.where(~is_control)[0]

        if len(ptb_indices) == 0:
            raise ValueError("No perturbed cells found in AnnData object.")
        if len(ctrl_indices) == 0:
            raise ValueError("No control cells found in AnnData object.")

        logger.info(
            f"Found {len(ctrl_indices)} control cells and {len(ptb_indices)} perturbed cells."
        )

        # Filter based on mode
        if self.mode == "training":
            # Training mode: Filter to only include single perturbations
            single_perturbation_mask = []
            for idx in ptb_indices:
                pert_str = self.adata.obs[self.perturbation_key].iloc[idx]
                # Single perturbation if it doesn't contain a comma
                is_single = pert_str and pert_str != "" and "," not in pert_str
                single_perturbation_mask.append(is_single)

            single_perturbation_mask = np.array(single_perturbation_mask)
            final_ptb_indices = ptb_indices[single_perturbation_mask]

            if len(final_ptb_indices) == 0:
                raise ValueError("No single perturbation cells found in AnnData object.")

            logger.info(
                f"Training mode: Filtered to {len(final_ptb_indices)} single perturbation cells."
            )

        elif self.mode == "prediction":
            # Prediction mode: Use specified indices (can include double perturbations)
            if self.prediction_indices is not None:
                # If indices are provided, they should already be filtered appropriately
                # by the caller (e.g., SENADataSplitter), so trust them
                final_ptb_indices = np.array(self.prediction_indices)

                # Verify these are actually perturbation indices
                valid_ptb_indices = np.intersect1d(final_ptb_indices, ptb_indices)
                if len(valid_ptb_indices) != len(final_ptb_indices):
                    if valid_ptb_indices.size > 0:
                        logger.warning(
                            f"Some provided indices are not perturbation cells. "
                            f"Using {len(valid_ptb_indices)}/{len(final_ptb_indices)}"
                            f"valid indices."
                        )
                        final_ptb_indices = valid_ptb_indices

                if len(final_ptb_indices) == 0:
                    logger.warning("No pertubred cells in indices")

                # Count single vs double perturbations for logging
                single_count = 0
                double_count = 0
                for idx in final_ptb_indices:
                    pert_str = self.adata.obs[self.perturbation_key].iloc[idx]
                    if pert_str and pert_str != "":
                        if "," in pert_str:
                            double_count += 1
                        else:
                            single_count += 1

                logger.info(
                    f"Prediction mode: Using {len(final_ptb_indices)} pre-filtered cells"
                    f"({single_count} single, {double_count} double perturbations)."
                )
            else:
                # No specific indices provided, use all perturbed cells
                final_ptb_indices = ptb_indices
                logger.info(
                    f"Prediction mode: Using all {len(final_ptb_indices)} perturbed cells."
                )
        else:
            raise ValueError(f"Unknown mode: {self.mode}. Must be 'training' or 'prediction'.")

        self.actual_training_indices = final_ptb_indices

        # Create matched pairs by sampling controls for the selected perturbations
        n_ptb = len(final_ptb_indices)
        random_ctrl_indices = np.random.choice(ctrl_indices, n_ptb, replace=True)

        # Store the matched data gene expression and the gene names of the perturbations
        self.matched_ctrl_X = self.adata.X[random_ctrl_indices]  # Control expression
        self.perturbed_X = self.adata.X[final_ptb_indices]  # Perturbed expression
        self.perturbation_strings = (
            self.adata.obs[self.perturbation_key].iloc[final_ptb_indices].values
        )

        # Process perturbation strings to numerical indices, considering all
        # genes (single + double)
        self._process_perturbation_strings()

    def _process_perturbation_strings(self):
        """
        Convert perturbation annotations into numerical intervention matrices.

        This method:
        1. Extracts unique perturbed genes from ALL perturbed annotation strings
        (single and double) for vocabulary
        2. Creates a mapping between genes and numerical indices
        3. Constructs intervention matrices representing perturbation states for training data
        (single only)
        4. Handles single gene perturbations only in the final matrix

        Notes
        -----
        The intervention matrix format:
        - Shape: (n_single_perturbations, n_intervention_genes)
        - Values: -1 for controls, 0 for non-perturbed genes, 1 for perturbed genes

        Vocabulary Construction:
        The method builds a complete vocabulary from ALL perturbations (single and double)
        to ensure compatibility during inference, even though training uses only single
        perturbations.
        """
        logger.info("Processing perturbation strings to numerical format.")

        # Extract all unique genes from ALL perturbation strings (including double perturbations)
        # This ensures the vocabulary includes genes from double perturbations even though
        # we only train on single perturbations
        all_perturbed_genes = set()

        # Get all perturbation strings from the entire dataset (not just training subset)
        all_perturbation_strings = self.adata.obs[self.perturbation_key].values

        for pert_str in all_perturbation_strings:
            if pert_str and pert_str != "":
                genes = [g.strip() for g in pert_str.split(",")]
                all_perturbed_genes.update(genes)

        # Create gene-to-index mapping using the complete vocabulary(Ex: {ATF4:[15]})
        sorted_genes = sorted(all_perturbed_genes)
        gene_to_idx = {gene: idx for idx, gene in enumerate(sorted_genes)}
        n_intervention_genes = len(sorted_genes)

        # this is improtant as it shapes the size of one of the encoder´s first layer.
        logger.info(
            f"Found {n_intervention_genes} unique intervention genes from all perturbations:"
            f"{sorted_genes[:10]}{'...' if n_intervention_genes > 10 else ''}"
        )

        # Count single vs double perturbations for logging
        single_count = sum(1 for s in all_perturbation_strings if s and s != "" and "," not in s)
        double_count = sum(1 for s in all_perturbation_strings if s and s != "" and "," in s)
        logger.info(
            f"Total perturbations in dataset: {single_count} single, {double_count} double"
        )
        if self.mode == "training":
            logger.info(
                f"Training will use only {len(self.perturbation_strings)} single perturbations"
            )

        # Create numerical intervention matrix for training data (single perturbations only)
        # Shape: (n_single_perturbations, n_intervention_genes)
        # Values: -1 for controls, 0 for non-perturbed genes, 1 for perturbed genes
        intervention_matrix = np.full(
            (len(self.perturbation_strings), n_intervention_genes), -1, dtype=np.int32
        )

        # Process only the single perturbations that will be used for training
        for cell_idx, pert_str in enumerate(self.perturbation_strings):
            if pert_str and pert_str != "":  # Should all be single perturbations
                # Verify it's actually a single perturbation
                if "," in pert_str and self.mode == "training":
                    logger.warning(
                        f"Double perturbation '{pert_str}' found"
                        f"in training data - this should not happen"
                    )
                    # continue;

                # Set all genes to 0 first (non-perturbed)
                intervention_matrix[cell_idx, :] = 0

                # Set perturbed gene to 1
                gene = pert_str.strip()
                if "," in gene:
                    gene = gene.split(",")
                    for g in gene:
                        if g in gene_to_idx:
                            intervention_matrix[cell_idx, gene_to_idx[g]] = 1
                        else:
                            logger.warning(f"Gene '{g}' in perturbation not found in gene mapping")
                else:
                    if gene in gene_to_idx:
                        intervention_matrix[cell_idx, gene_to_idx[gene]] = 1
                    else:
                        logger.warning(f"Gene '{gene}' in perturbation not found in gene mapping")

        # Store the processed data
        self.intervention_matrix = intervention_matrix
        self.intervention_genes = sorted_genes
        self.gene_to_intervention_idx = gene_to_idx
        self.n_intervention_genes = n_intervention_genes

    def reshuffle_controls(self):
        """
        Reshuffle control assignments for the next epoch, respecting current mode.

        Re-creates the matched pairs with new random sampling to enhance model
        generalization by preventing memorization of specific control-perturbation
        relationships.

        Notes
        -----
        This method applies the same filtering logic as in `_create_matched_dataset`
        to maintain consistency between initial setup and reshuffling. The method
        prevents overfitting by ensuring different control-perturbation pairings
        across epochs while maintaining the same perturbation distribution.
        """
        # Re-create the matched pairs with new random sampling
        is_control = self.adata.obs[self.perturbation_key] == ""
        ctrl_indices = np.where(is_control)[0]
        ptb_indices = np.where(~is_control)[0]

        # Apply same filtering logic as in _create_matched_dataset. Both methods correspond
        # to the same class so the self atributes are the same.We could use the control indices
        # created during the _create_matched_dataset, without needed to reshufle.
        if self.mode == "training":
            # Filter to only include single perturbations
            single_perturbation_mask = []
            for idx in ptb_indices:
                pert_str = self.adata.obs[self.perturbation_key].iloc[idx]
                is_single = pert_str and pert_str != "" and "," not in pert_str
                single_perturbation_mask.append(is_single)

            single_perturbation_mask = np.array(single_perturbation_mask)
            final_ptb_indices = ptb_indices[single_perturbation_mask]

        elif self.mode == "prediction":
            if self.prediction_indices is not None:
                final_ptb_indices = np.intersect1d(ptb_indices, self.prediction_indices)
            else:
                final_ptb_indices = ptb_indices

        n_ptb = len(final_ptb_indices)
        random_ctrl_indices = np.random.choice(ctrl_indices, n_ptb, replace=True)
        self.matched_ctrl_X = self.adata.X[random_ctrl_indices]

    def __iter__(self):
        """
        Make the dataloader iterable for batch processing.

        When you write 'for batch in dataloader:' in your training code,
        Python calls this __iter__ method to provide batches of processed data.

        Yields
        ------
        dict of str to torch.Tensor
            Dictionary containing:
            - X_KEY : torch.Tensor
                Control expression data, shape (batch_size, n_genes).
            - LABELS_KEY : torch.Tensor
                Perturbed expression data, shape (batch_size, n_genes).
            - CAT_COVS_KEY : torch.Tensor
                Intervention matrix, shape (batch_size, n_intervention_genes).

        Notes
        -----
        This method has two phases:

        PHASE 1 - SETUP (runs once per epoch):
        - Gets the total number of samples (e.g., 52,047 training pairs)
        - Creates a list of indices [0, 1, 2, ..., 52046]
        - Shuffles these indices randomly for this epoch
        - Sets the batch size (e.g., 32 samples per batch)

        PHASE 2 - BATCH CREATION (runs multiple times per epoch):
        - The 'for' loop below divides data into batches of 32 samples each
        - Each time your training code asks for the next batch, this loop runs ONE iteration
        - The 'yield' statement pauses here and returns the current batch
        - When the next batch is requested, the loop resumes from where it left off
        - This continues until all ~1627 batches (52047÷32) are processed
        """
        # Use shape[0] for sparse matrix compatibility
        n_samples = self.matched_ctrl_X.shape[0]
        indices = np.arange(n_samples)

        # Use our stored shuffle parameter
        if getattr(self, "_custom_shuffle", self._shuffle):
            np.random.shuffle(indices)

        # Use our stored batch_size
        # Get the batch size
        batch_size = getattr(self, "_custom_batch_size", self._batch_size)
        if batch_size is None:
            batch_size = self._batch_size

        batch_count = 0

        # Once we have the matched control and perturbations we need to divide them in batches
        for i in range(0, n_samples, batch_size):
            batch_indices = indices[i : i + batch_size]

            # Get batch data and obtain the control and pertubed expresion and
            # the intervention matrix
            batch_ctrl_X = self.matched_ctrl_X[batch_indices]
            batch_perturbed_X = self.perturbed_X[batch_indices]
            batch_intervention_matrix = self.intervention_matrix[batch_indices]

            # We create a distioanry of three tensors
            tensors = {
                REGISTRY_KEYS.X_KEY: torch.from_numpy(
                    batch_ctrl_X.toarray() if hasattr(batch_ctrl_X, "toarray") else batch_ctrl_X
                ).float(),
                REGISTRY_KEYS.LABELS_KEY: torch.from_numpy(
                    batch_perturbed_X.toarray()
                    if hasattr(batch_perturbed_X, "toarray")
                    else batch_perturbed_X
                ).float(),
                REGISTRY_KEYS.CAT_COVS_KEY: torch.from_numpy(batch_intervention_matrix).long(),
            }

            # Debug output for first batch
            if batch_count == 0:
                batch_count += 1
            # This freeses the for loop an returns the current batch until the generator
            # funciton its called again.
            yield tensors

    def __len__(self):
        """
        Return number of batches in the dataloader.

        This is a special method that Python calls when we do len(dataloader).
        It calculates the total number of batches based on the number of matched
        control-perturbation pairs and the batch size.

        Returns
        -------
        int
            Number of batches in the dataloader.

        Notes
        -----
        The calculation uses ceiling division to ensure all samples are included:
        (n_samples + batch_size - 1) // batch_size
        This handles cases where n_samples is not evenly divisible by batch_size.
        """
        n_samples = self.matched_ctrl_X.shape[0]
        batch_size = getattr(self, "_custom_batch_size", self._batch_size)
        if batch_size is None:
            batch_size = self._batch_size
        return (n_samples + batch_size - 1) // batch_size

    def get_intervention_vocabulary(self):
        """
        Return the complete intervention vocabulary including genes from double perturbations.

        Provides access to the complete gene vocabulary used for intervention encoding,
        including genes that appear only in combination perturbations. This is essential
        for inference on novel perturbation combinations.

        Returns
        -------
        dict
            Dictionary containing:
            - 'genes' : list
                List of all intervention genes.
            - 'gene_to_idx' : dict
                Mapping from gene names to indices.
            - 'n_genes' : int
                Number of unique intervention genes.

        Notes
        -----
        The vocabulary includes all genes from both single and double perturbations
        to ensure compatibility during inference, even though training may use only
        single perturbations.
        """
        return {
            "genes": self.intervention_genes,
            "gene_to_idx": self.gene_to_intervention_idx,
            "n_genes": self.n_intervention_genes,
        }

    def create_intervention_vector(self, perturbation_string):
        """
        Create an intervention vector for any perturbation string (single or double).

        This method allows creating intervention vectors for double perturbations
        at inference time, even though the model was trained only on single perturbations.
        It enables prediction of combination effects by encoding multiple simultaneous
        gene perturbations.

        Parameters
        ----------
        perturbation_string : str
            Perturbation string (e.g., "GENE1" or "GENE1,GENE2").

        Returns
        -------
        np.ndarray
            Intervention vector of shape (n_intervention_genes,).
            Values: 0 for non-perturbed genes, 1 for perturbed genes.

        Notes
        -----
        The method handles various perturbation formats:
        - Empty string or "": Returns zero vector (control condition)
        - Single gene: "GENE1" -> one-hot vector with 1 at GENE1 position
        - Multiple genes: "GENE1,GENE2" -> multi-hot vector with 1s at both positions
        """
        if not perturbation_string or perturbation_string == "":
            # Control case - return all zeros
            return np.zeros(self.n_intervention_genes, dtype=np.int32)

        # Initialize intervention vector
        intervention_vector = np.zeros(self.n_intervention_genes, dtype=np.int32)

        # Parse perturbation string
        genes = [g.strip() for g in perturbation_string.split(",")]

        # Set perturbed genes to 1
        for gene in genes:
            if gene in self.gene_to_intervention_idx:
                intervention_vector[self.gene_to_intervention_idx[gene]] = 1
            else:
                logger.warning(f"Gene '{gene}' not found in intervention vocabulary")

        return intervention_vector


class ControlReshuffleCallback(Callback):
    """
    A PyTorch Lightning callback that performs dynamic cell reassignment during training.

    This callback enhances model generalization by randomizing the control-perturbation
    pairs at the start of each training epoch. It prevents the model from memorizing
    specific control-perturbation relationships and ensures learning of general
    perturbation effects.

    Notes
    -----
    Key Functions:
    - Prevents overfitting by randomizing control cell assignments
    - Enhances learning of generalizable perturbation effects
    - Reduces bias from specific control-perturbation pairings
    - Seamlessly integrates with PyTorch Lightning's training loop

    Implementation:
    The callback operates by calling the reshuffle_controls() method of the
    dataloader at the start of each epoch, ensuring random reassignment of
    control cells while maintaining the overall data structure.

    Scientific Rationale:
    In single-cell perturbation experiments, the choice of control cells paired
    with each perturbed cell can introduce systematic biases. By reshuffling
    these pairings each epoch, the model learns perturbation effects that
    generalize across different cellular contexts.
    """

    def __init__(self):
        """
        Initialize control reshuffle callback for SENA training workflows.

        Sets up the callback infrastructure for dynamic cell reassignment during
        training epochs, enabling enhanced generalization through randomized
        control-perturbation pairings.
        """
        super().__init__()

    def on_train_epoch_start(self, trainer, pl_module):
        """
        Reshuffle control assignments in the SENA dataloaders.

        Function that is called at the start of each epoch to perform dynamic
        control cell reassignment for enhanced generalization. Accesses the training
        dataloader and calls its reshuffle_controls method if available.

        Parameters
        ----------
        trainer : pytorch_lightning.Trainer
            PyTorch Lightning trainer managing the training process.
        pl_module : pytorch_lightning.LightningModule
            SENA model being trained.

        Notes
        -----
        This method attempts to access the train_dataloader through the trainer
        and calls its reshuffle_controls method if available. It also logs the
        current epoch for monitoring purposes. The reshuffling ensures that
        control-perturbation pairings vary across epochs to prevent overfitting.
        """
        # I am not ssure this would wok because reshuffle controls is a custom function
        # from our dataloader. But maybe Claudio just copied the name form an
        # exsitng function in the trainer
        # Method 2: Try to access through trainer's train_dataloader
        if hasattr(trainer, "train_dataloader") and trainer.train_dataloader is not None:
            train_dl = trainer.train_dataloader
            if hasattr(train_dl, "reshuffle_controls"):
                train_dl.reshuffle_controls()
                if hasattr(pl_module, "log"):
                    pl_module.log(
                        "control_reshuffle_epoch", float(trainer.current_epoch), on_epoch=True
                    )


def create_sena_dataloader(
    adata_manager: AnnDataManager,
    perturbation_key: str,
    indices: Sequence[int] | None = None,
    batch_size: int = 128,
    shuffle: bool = True,
    mode: str = "training",  # "training" or "prediction"
    **kwargs,
) -> SENADataLoader:
    """
    Factory function to create a SENA dataloader with specified configuration.

    Provides a convenient interface for creating SENADataLoader instances with
    standardized parameters and proper error handling. This factory function
    ensures consistent dataloader configuration across different parts of the
    SENA framework.

    Parameters
    ----------
    adata_manager : AnnDataManager
        AnnDataManager instance containing processed single-cell data.
    perturbation_key : str
        Key in adata.obs containing perturbation information.
    indices : Sequence[int], optional
        Indices to use (for prediction mode).
    batch_size : int, default 128
        Batch size for training.
    shuffle : bool, default True
        Whether to shuffle the data.
    mode : str, default "training"
        Mode for dataloader: "training" (single perturbations only) or
        "prediction" (all specified indices).
    **kwargs
        Additional arguments passed to SENADataLoader.

    Returns
    -------
    SENADataLoader
        Configured SENADataLoader instance ready for training or inference.

    Notes
    -----
    This factory function standardizes the creation of SENA dataloaders and
    provides a single point for configuration validation and error handling.
    It ensures consistent behavior across different usage contexts.
    """
    return SENADataLoader(
        adata_manager=adata_manager,
        perturbation_key=perturbation_key,
        batch_size=batch_size,
        shuffle=shuffle,
        mode=mode,
        indices=indices,
        **kwargs,
    )


class SENADataSplitter(DataSplitter):
    """
    Custom data splitter for SENA perturbation analysis with stratified splitting.

    This splitter ensures proper train/validation separation for perturbation experiments
    while maintaining scientific validity and preventing data leakage. It implements
    stratified sampling to preserve perturbation type distributions across splits
    and handles the unique requirements of control-perturbation experimental designs.

    Parameters
    ----------
    adata_manager : AnnDataManager
        AnnDataManager instance containing single-cell perturbation data.
    perturbation_key : str
        Key in adata.obs containing perturbation information.
    train_size : float, default 0.9
        Size of training set (default 0.9 = 90%).
    validation_size : float, optional
        Size of validation set (auto-calculated if None).
    shuffle_set_split : bool, default True
        Whether to shuffle when splitting.
    batch_size : int, default 128
        Batch size for dataloaders.
    **kwargs
        Additional arguments passed to DataSplitter.

    Attributes
    ----------
    perturbation_key : str
        Key in adata.obs containing perturbation information.
    train_ptb_idx : np.ndarray
        Indices of training perturbation cells.
    val_ptb_idx : np.ndarray
        Indices of validation perturbation cells.

    Notes
    -----
    Key Features:
    1. Perturbation-Only Splitting:
       - Only single perturbation cells are split between train/validation
       - Uses stratified sampling to maintain perturbation type distribution
       - Controls are available to both train and validation dataloaders

    2. scvi-tools Integration:
       - Seamlessly integrates with UnsupervisedTrainingMixin
       - Follows standard scvi-tools validation patterns
       - No separate validation function needed - automatic during training

    3. Scientific Validity:
       - Prevents overfitting by using separate perturbation cells for validation
       - Maintains proper evaluation methodology
       - Control sharing is acceptable since controls represent baseline state

    Experimental Design Rationale:
    In perturbation experiments, the key is to evaluate the model's ability to
    predict unseen perturbation effects. Therefore, validation must use different
    perturbed cells than training, while control cells can be shared since they
    represent the common baseline state.
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        perturbation_key: str,
        train_size: float = 0.9,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        batch_size: int = 128,
        **kwargs,
    ):
        """
        Initialize SENA data splitter with stratified perturbation splitting.

        Sets up the data splitting infrastructure for SENA training, ensuring
        proper separation of perturbation cells while allowing control cell
        sharing between training and validation sets.

        Parameters
        ----------
        adata_manager : AnnDataManager
            AnnDataManager instance containing single-cell perturbation data.
        perturbation_key : str
            Key in adata.obs containing perturbation information.
        train_size : float, default 0.9
            Size of training set (default 0.9 = 90%).
        validation_size : float, optional
            Size of validation set (auto-calculated if None).
        shuffle_set_split : bool, default True
            Whether to shuffle when splitting.
        batch_size : int, default 128
            Batch size for dataloaders.
        **kwargs
            Additional arguments passed to DataSplitter.
        """
        self.perturbation_key = perturbation_key
        self._batch_size = batch_size

        # Initialize parent DataSplitter
        super().__init__(
            adata_manager=adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            batch_size=batch_size,
            **kwargs,
        )

    def setup(self, stage: str = None):
        """
        Setup data splits for SENA with proper perturbation-only splitting.

        Performs stratified splitting of single perturbation cells between training
        and validation sets while ensuring control cells remain available to both.
        This approach maintains scientific validity by testing on unseen perturbation
        effects while allowing baseline comparison.

        Parameters
        ----------
        stage : str, optional
            Training stage (unused but required for compatibility).

        Raises
        ------
        ValueError
            If no single perturbation cells found for training.

        Notes
        -----
        Splitting Strategy:
        1. Identifies single perturbation cells (excludes controls and double perturbations)
        2. Performs stratified split to maintain perturbation type proportions
        3. Assigns training and validation indices for perturbation cells only
        4. Control cells remain accessible to both training and validation dataloaders

        The stratification ensures that rare perturbation types are represented
        in both training and validation sets proportionally.
        """
        adata = self.adata_manager.adata

        # Get perturbation information
        perturbations = adata.obs[self.perturbation_key].values

        # Identify single perturbation cells (training data)
        is_control = perturbations == ""
        is_single_perturbation = (~is_control) & (~pd.Series(perturbations).str.contains(","))

        # Get indices for single perturbation cells only
        single_perturbation_indices = np.where(is_single_perturbation)[0]

        if len(single_perturbation_indices) == 0:
            raise ValueError("No single perturbation cells found for training")

        # Perform stratified (same proportion of each perturbation on train and validation)
        # split on single perturbation cells only
        single_perturbations = perturbations[single_perturbation_indices]

        from sklearn.model_selection import train_test_split

        train_ptb_idx, val_ptb_idx = train_test_split(
            single_perturbation_indices,
            test_size=1 - self.train_size,
            stratify=single_perturbations,
            random_state=42,
            shuffle=self.shuffle_set_split,
        )

        # Store the perturbation indices for train/validation
        self.train_ptb_idx = train_ptb_idx
        self.val_ptb_idx = val_ptb_idx

        # For compatibility with parent class, set train_idx and val_idx
        self.train_idx = train_ptb_idx
        self.val_idx = val_ptb_idx
        self.test_idx = np.array([])  # No test set for SENA

        logger.info(
            f"SENA data split: {len(train_ptb_idx)} train perturbations, {len(val_ptb_idx)}"
            f"validation perturbations"
        )
        logger.info("Control cells will be shared between train/validation sets")

        # Log perturbation distribution
        train_ptb_types = np.unique(perturbations[train_ptb_idx], return_counts=True)
        val_ptb_types = np.unique(perturbations[val_ptb_idx], return_counts=True)
        logger.info(f"Train perturbation types: {len(train_ptb_types[0])}")
        logger.info(f"Val perturbation types: {len(val_ptb_types[0])}")

    def train_dataloader(self):
        """
        Create training dataloader with only training perturbation indices.

        Generates a SENADataLoader configured for training mode with the subset
        of perturbation cells assigned to the training split. Control cells are
        accessible through the dataloader's dynamic matching mechanism.

        Returns
        -------
        SENADataLoader
            Training dataloader configured for single perturbations only.

        Notes
        -----
        The dataloader uses training mode which restricts processing to single
        perturbations and enables shuffling for proper stochastic gradient descent.
        Control cells are dynamically matched to perturbation cells during iteration.
        """
        return SENADataLoader(
            adata_manager=self.adata_manager,
            perturbation_key=self.perturbation_key,
            batch_size=self._batch_size,
            shuffle=True,  # Always shuffle for training
            mode="training",  # Training mode - only single perturbations
            prediction_indices=self.train_ptb_idx,  # Use only training perturbation indices
        )

    def val_dataloader(self):
        """
        Create validation dataloader with only validation perturbation indices.

        Generates a SENADataLoader configured for validation with the subset
        of perturbation cells assigned to the validation split. Uses the same
        control cells as training but with different perturbation targets.

        Returns
        -------
        SENADataLoader
            Validation dataloader configured for single perturbations only.

        Notes
        -----
        The validation dataloader uses training mode (single perturbations only)
        but with shuffle=False for consistent validation metrics. Control cells
        are shared with training, which is scientifically appropriate since
        controls represent the baseline cellular state.
        """
        # Check val_ptb_idx instead of val_idx since that's what we actually use for validation
        return SENADataLoader(
            adata_manager=self.adata_manager,
            perturbation_key=self.perturbation_key,
            batch_size=self._batch_size,
            shuffle=False,  # No shuffling for validation
            mode="training",  # Training mode - only single perturbations
            prediction_indices=self.val_ptb_idx,  # Use only validation perturbation indices
        )
