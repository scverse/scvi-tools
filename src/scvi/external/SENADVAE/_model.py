import logging
import os
import random
from collections import defaultdict
from collections.abc import Sequence

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from anndata import AnnData
from scipy import stats
from statsmodels.stats.multitest import multipletests

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import LayerField
from scvi.model.base import BaseModelClass, UnsupervisedTrainingMixin
from scvi.train import TrainingPlan, TrainRunner
from scvi.utils import setup_anndata_dsp

from ._dataloader import ControlReshuffleCallback, SENADataLoader, SENADataSplitter
from ._module import SENAModule
from ._training_plan import (
    LossWeightScheduler,
    SENABatchProgressBar,
    SENATrainingPlan,
    TemperatureScheduler,
)

logger = logging.getLogger(__name__)


class SENADVAE(UnsupervisedTrainingMixin, BaseModelClass):
    """
    Structural Equation Network Analysis Variational AutoEncoder for causal perturbation modeling.

    SENADVAE implements a novel approach to understanding causal relationships in single-cell
    perturbation experiments by combining Variational AutoEncoders with biological pathway
    constraints and causal graph learning. This model is specifically designed for analyzing
    CRISPR screens, drug perturbations, and other intervention studies in single-cell genomics.

    Parameters
    ----------
    adata : AnnData
        Annotated data object containing single-cell expression data with perturbation annotations.
        Must be pre-registered via setup_anndata with control/perturbation labels.
    go_file_path : str
        Path to Gene Ontology annotation file containing pathway definitions for biological
        constraints.
        Typically in format: pathway_id<tab>gene_name for GO term memberships.
    go_gene_map_path : str
        Path to gene-to-GO mapping file enabling pathway constraint enforcement
        in NetworkActivity layers.
        Format: gene_name<tab>go_term_id for comprehensive gene-pathway associations.
    n_latent : int, default 64
        Dimensionality of VAE latent space representing pathway activity levels.
        Should balance expressiveness with interpretability (32-128 typical range).
    sena_lambda : float, default 0.1
        L1 regularization strength for causal graph sparsity, promoting interpretable
        pathway interactions.
        Higher values encourage sparser causal structures (0.01-1.0 typical range).
    n_hidden_encoder : int, default 256
        Hidden layer width in pathway-constrained encoder networks.
        Should accommodate pathway complexity while maintaining computational efficiency.
    n_hidden_decoder : int, default 256
        Hidden layer width in intervention-aware decoder networks.
        Larger values may improve reconstruction quality for complex perturbation effects.
    n_hidden_interv : int, default 256
        Hidden layer width in intervention network responsible for modeling perturbation-specific
        effects.
        Controls capacity for learning complex intervention-response relationships.
    mode : str, default "sena"
        Training mode selection: "sena" enables pathway constraints and causal learning,
        "normal" uses standard VAE training without biological constraints.

    Attributes
    ----------
    perturbation_key : str
        Key in adata.obs containing perturbation information.
    rel : dict
        Gene-to-pathway relationship mapping dictionary.
    go_map : pd.DataFrame
        Gene-to-GO mapping dataframe.
    go_dict : dict
        GO term to index mapping.
    gen_dict : dict
        Gene to index mapping.
    mapping_dict : dict
        Gene symbol to ensemble ID mapping.
    module : SENAModule
        The underlying neural network module.

    Notes
    -----
    Scientific Innovation:
    - **Biological Pathway Integration**: Enforces Gene Ontology constraints via
    NetworkActivity layers
    - **Causal Graph Learning**: Discovers intervention-specific causal relationships between
    pathways
    - **Intervention Modeling**: Explicitly models control vs. perturbation differences
    - **Distribution Matching**: Uses MMD or MSE to ensure realistic perturbation predictions

    Model Architecture:
    1. **Pathway-Constrained Encoder**: Maps gene expression to biologically meaningful
    latent space
    2. **Causal Graph Network**: Learns intervention-specific pathway interactions
    3. **Intervention-Aware Decoder**: Reconstructs expression under control and
    perturbation conditions
    4. **Multi-Component Loss**: Balances reconstruction, prediction, and biological constraints

    Key Applications:
    - CRISPR knockout/knockdown effect prediction
    - Drug mechanism of action discovery
    - Pathway interaction mapping under perturbations
    - Single-cell intervention response modeling
    - Causal biomarker identification

    Examples
    --------
    >>> import scanpy as sc
    >>> from scvi.external import SENADVAE
    >>> # Load perturbation screen data
    >>> adata = sc.read_h5ad("crispr_screen.h5ad")
    >>> # adata.obs["perturbation"] contains: "", "KRAS", "TP53", "KRAS,TP53", etc.
    >>> # Register data with SENA-specific preprocessing
    >>> SENADVAE.setup_anndata(adata, perturbation_key="perturbation")
    >>> # Initialize model with pathway constraints
    >>> model = SENADVAE(
    ...     adata,
    ...     go_file_path="GO_pathways.tsv",
    ...     go_gene_map_path="gene_GO_mapping.tsv",
    ...     n_latent=64,
    ...     sena_lambda=0.1,
    ... )
    >>> # Train with curriculum learning and progress monitoring
    >>> model.train(max_epochs=100, batch_size=256)
    >>> # Extract pathway activity representations
    >>> latent_pathways = model.get_latent_representation()
    >>> # Predict perturbation effects
    >>> perturb_pred = model.predict_perturbation(["KRAS", "TP53"])
    """

    # Core scvi-tools integration components for SENA training infrastructure
    _module_cls = SENAModule  # Neural architecture with pathway constraints
    _training_plan_cls = SENATrainingPlan  # Multi-component loss and curriculum learning
    _train_runner_cls = TrainRunner  # Standard scvi trainer managing training loops
    _data_splitter_cls = SENADataSplitter  # Control-perturbation aware data splitting

    @staticmethod
    def set_seeds(seed: int) -> None:
        """
        Configure deterministic random number generation for reproducible SENA training.

        Ensures reproducible results across different runs of perturbation modeling experiments,
        critical for validating biological discoveries and comparing model configurations.
        Sets seeds for all relevant random number generators used in training pipeline.

        Parameters
        ----------
        seed : int
            Random seed for deterministic model initialization and training.
        """
        torch.manual_seed(seed)
        np.random.seed(seed)
        random.seed(seed)

        # Configure GPU determinism if available
        if torch.cuda.is_available():
            torch.cuda.manual_seed(seed)
            torch.cuda.manual_seed_all(seed)

        # Set PyTorch Lightning seeds for comprehensive determinism
        try:
            import lightning.pytorch as pl

            pl.seed_everything(seed, workers=True)
        except ImportError:
            pass

    def __init__(
        self,
        adata: AnnData,
        go_file_path: str,
        go_gene_map_path: str,
        gene_symb_ensemble_path: str,
        n_latent: int = 105,
        n_go_thresh: int = 5,
        sena_lambda: float = 0.0,
        n_hidden_encoder: int = 128,
        n_hidden_decoder: int = 128,
        n_hidden_interv: int = 128,
        seed: int | None = None,
        mode: str = "sena",
        **model_kwargs,
    ):
        """
        Initialize SENA model for causal perturbation analysis with biological constraints.

        Sets up the complete SENA architecture including pathway-constrained encoders,
        causal graph learning networks, and intervention-aware decoders. Configures
        biological constraints from Gene Ontology annotations and prepares model
        for training on single-cell perturbation data.

        Parameters
        ----------
        adata : AnnData
            Pre-registered single-cell perturbation dataset with control/treatment annotations.
            Must contain perturbation metadata in .obs and normalized expression in .X.
        go_file_path : str
            Path to Gene Ontology pathway definition file for biological constraint enforcement.
            Expected format: pathway_id<tab>gene_symbol for comprehensive pathway coverage.
        go_gene_map_path : str
            Path to gene-to-GO term mapping file enabling pathway constraint application.
            Expected format: gene_symbol<tab>go_term_id for complete gene annotation.
        gene_symb_ensemble_path : str
            Path to gene symbol to ensemble ID mapping file.
        n_latent : int, default 105
            Latent space dimensionality representing pathway activity levels.
            Should approximate number of relevant biological pathways (50-200 typical).
        n_go_thresh : int, default 5
            Minimum gene count threshold for including GO pathways in constraints.
            Filters out very small pathways that may introduce noise (3-10 typical range).
        sena_lambda : float, default 0.0
            L1 regularization strength for causal graph sparsity in pathway interactions.
            Higher values promote sparser, more interpretable causal structures.
        n_hidden_encoder : int, default 512
            Hidden layer width in pathway-constrained encoder networks.
            Larger values increase model capacity for complex expression patterns.
        n_hidden_decoder : int, default 128
            Hidden layer width in intervention-aware decoder networks.
            Should balance reconstruction quality with computational efficiency.
        n_hidden_interv : int, default 256
            Hidden layer width in intervention prediction networks.
            Controls capacity for modeling perturbation-specific effects.
        seed : int, optional
            Random seed for reproducible model initialization and training.
            Essential for validating biological discoveries across experiments.
        mode : str, default "sena"
            Training mode: "sena" enables pathway constraints and causal learning,
            "normal" uses standard VAE without biological constraints.
        **model_kwargs
            Additional arguments passed to parent BaseModelClass initialization.

        Notes
        -----
        Architecture Initialization:
        - Loads Gene Ontology pathway constraints for NetworkActivity layers
        - Configures VAE latent dimensionality for pathway representation
        - Sets up multi-component loss function with biological regularization
        - Initializes causal graph learning with sparsity constraints
        """
        # Set deterministic training if seed provided
        if seed is not None:
            self.set_seeds(seed)

        # Initialize scvi-tools base model infrastructure
        super().__init__(adata)

        logger.info("Initializing SENA model with Gene Ontology pathway constraints.")

        # Load and process Gene Ontology pathway annotations for biological constraints

        # Load gene-to-GO mapping file containing pathway memberships
        go_map = pd.read_csv(go_gene_map_path, sep="\t")
        go_map.columns = ["GO_id", "ensembl_id"]

        # Filter to only genes present in the expression dataset
        go_map = go_map[go_map["ensembl_id"].isin(adata.var_names)]

        # Load selected GO pathways for biological constraint enforcement
        selected_gos = pd.read_csv(go_file_path, sep="\t")["PathwayID"].values.tolist()

        # Restrict to pathways from the curated pathway selection
        go_map = go_map[go_map["GO_id"].isin(selected_gos)]

        # Apply pathway size filtering to ensure robust constraint learning
        go_counts = go_map["GO_id"].value_counts()
        genesets_in = go_counts[go_counts >= n_go_thresh].index
        logger.info(
            f"Applying pathway size filter: {len(genesets_in)} pathways with ≥{n_go_thresh} genes."
        )

        go_map = go_map[go_map["GO_id"].isin(genesets_in)]
        # We obtain a sorted list of the GOs that pases the filtering
        gos = sorted(go_map["GO_id"].unique())

        # Create gene-pathway relationship dictionary
        genes = adata.var.index.values
        # dictinaries with gos and genes as their keys, and indices as their valies
        go_dict = dict(zip(gos, range(len(gos)), strict=False))
        gen_dict = dict(zip(genes, range(len(genes)), strict=False))
        # We create a dictionary that relates the gos and ensemble ids
        # Create filtered pathway-gene relationship mapping for NetworkActivity constraints
        rel_dict = defaultdict(list)
        gene_set, go_set = set(genes), set(gos)
        self.go_dict = go_dict
        self.gen_dict = gen_dict

        mapping = pd.read_csv(gene_symb_ensemble_path, sep="\t")
        self.mapping_dict = dict(
            zip(mapping["external_gene_name"], mapping["ensembl_gene_id"], strict=False)
        )

        # Build gene-to-pathway mapping dictionary for biological constraint enforcement
        for go, gen in zip(go_map["GO_id"], go_map["ensembl_id"], strict=False):
            if (gen in gene_set) and (go in go_set):
                # Map each gene index to its associated pathway indices for NetworkActivity layers
                rel_dict[gen_dict[gen]].append(go_dict[go])

        # Extract intervention metadata from dataset preprocessing
        n_intervention_genes = adata.uns.get("n_intervention_genes", 1)

        # Retrieve perturbation annotation key from setup_anndata configuration
        self.perturbation_key = adata.uns.get("_sena_perturbation_key")
        if self.perturbation_key is None:
            raise ValueError(
                "Perturbation key not found in AnnData.uns. "
                "Ensure setup_anndata() was called with perturbation_key parameter."
            )

        # Set latent dimensionality to match number of biological pathways
        n_latent = n_intervention_genes
        logger.info(f"Configuring SENA with {n_latent} pathway-constrained latent dimensions.")

        self.rel = rel_dict

        self.go_map = go_map
        self.gos = gos

        # Initialize SENA neural architecture with biological constraints
        self.module = self._module_cls(
            n_input=adata.n_vars,  # Number of genes in expression matrix
            n_latent=n_latent,  # Pathway activity latent dimensions
            n_cat_covs=n_intervention_genes,  # Number of intervention categories
            n_categories_interv=n_intervention_genes,  # Intervention encoding dimensionality
            gos=gos,  # Pathway identifiers for constraints
            rel_dict=rel_dict,  # Gene-pathway relationship mapping
            n_hidden_encoder=n_hidden_encoder,  # Encoder network capacity
            n_hidden_decoder=n_hidden_decoder,  # Decoder network capacity
            n_hidden_interv=n_hidden_interv,  # Intervention network capacity
            mode=mode,  # Training mode (sena vs normal)
            sena_lambda=sena_lambda,  # L1 sparsity regularization strength
            **model_kwargs,
        )

        # Configure model summary for inspection and debugging
        self._model_summary_string = (
            f"SENA Model Summary:\n"
            f"├── Latent Pathways: {n_latent}\n"
            f"├── Encoder Hidden: {n_hidden_encoder}\n"
            f"├── Decoder Hidden: {n_hidden_decoder}\n"
            f"├── Intervention Hidden: {n_hidden_interv}\n"
            f"├── Training Mode: {mode}\n"
            f"├── L1 Regularization: {sena_lambda}\n"
            f"├── Input Genes: {adata.n_vars}\n"
            f"└── Intervention Categories: {n_intervention_genes}"
        )

        # Store initialization parameters for model persistence and reproducibility
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        perturbation_key: str,
        layer: str | None = None,
        **kwargs,
    ) -> AnnData | None:
        """
        Configure AnnData object for SENA perturbation modeling with automatic preprocessing.

        This class method prepares single-cell perturbation data for SENA training by:
        1. Parsing perturbation annotations into numerical intervention matrices
        2. Validating data quality and perturbation coverage
        3. Extracting intervention gene metadata for model configuration
        4. Registering data fields with scvi-tools infrastructure

        The actual control-perturbation matching and batch preprocessing is handled by
        the custom SENADataLoader during training, preserving the original AnnData structure
        while enabling sophisticated intervention modeling workflows.

        Parameters
        ----------
        %(param_adata)s
        perturbation_key : str
            Column name in adata.obs containing perturbation identifiers.
            Must contain control ("") and perturbation gene name entries.
            Example: "perturbation", "treatment", "intervention"
        %(param_layer)s
        **kwargs
            Additional arguments passed to parent setup_anndata method.

        Returns
        -------
        AnnData or None
            Modified AnnData object with registered fields, or None if in-place modification.

        Notes
        -----
        Data Requirements:
        - Expression data in adata.X (normalized log counts recommended)
        - Perturbation annotations in adata.obs[perturbation_key]
        - Both control ("") and perturbed cells must be present
        - Gene names in adata.var_names matching GO annotation files

        Perturbation Annotation Format:
        - Control cells: empty string ("")
        - Single perturbations: gene name ("KRAS")
        - Multiple perturbations: comma-separated ("KRAS,TP53")
        - Supports any combination of single/multiple perturbations

        Examples
        --------
        >>> # Setup CRISPR screen data
        >>> adata.obs["perturbation"] = ["", "KRAS", "TP53", "KRAS,TP53", ...]
        >>> SENADVAE.setup_anndata(adata, perturbation_key="perturbation")
        >>> # Use specific expression layer
        >>> SENADVAE.setup_anndata(adata, perturbation_key="treatment", layer="log1p")
        """
        logger.info("Configuring AnnData for SENA perturbation modeling.")

        # Validate perturbation annotation presence
        if perturbation_key not in adata.obs.columns:
            raise ValueError(
                f"Perturbation key '{perturbation_key}' not found in adata.obs columns."
            )

        # Validate control and perturbation cell presence
        is_control = adata.obs[perturbation_key] == ""
        n_controls = is_control.sum()
        n_perturbed = (~is_control).sum()

        if n_perturbed == 0:
            raise ValueError(
                "No perturbed cells detected. Ensure perturbation annotations "
                "contain non-empty strings for perturbed cells."
            )
        if n_controls == 0:
            raise ValueError(
                "No control cells detected. Ensure control cells are annotated "
                "with empty strings ('') in perturbation column."
            )

        logger.info(
            f"Data validation complete: {n_controls} controls, {n_perturbed} perturbed cells."
        )

        # Store perturbation key for dataloader access
        adata.uns["_sena_perturbation_key"] = perturbation_key

        # Extract and process intervention gene metadata
        all_perturbed_genes = set()
        perturbation_strings = adata.obs[perturbation_key].values

        # Parse perturbation strings to extract unique intervention genes
        for pert_str in perturbation_strings:
            if pert_str and pert_str != "":
                # Split comma-separated gene names and clean whitespace
                genes = [g.strip() for g in pert_str.split(",")]
                all_perturbed_genes.update(genes)

        # Create sorted intervention gene list for consistent model configuration
        sorted_genes = sorted(all_perturbed_genes)
        n_intervention_genes = len(sorted_genes)

        logger.info(
            f"Extracted {n_intervention_genes} unique intervention targets: {sorted_genes[:5]}..."
        )

        # Store intervention metadata for model initialization
        adata.uns["intervention_genes"] = sorted_genes
        adata.uns["n_intervention_genes"] = n_intervention_genes

        # Configure scvi-tools data registration with SENA-specific fields
        setup_method_args = cls._get_setup_method_args(**locals())

        # Register expression data with scvi-tools infrastructure
        # Note: Perturbation annotations are handled by custom dataloader
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
        ]

        # Complete data registration with AnnDataManager
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

        return adata

    @staticmethod
    def analyze_perturbations(adata: AnnData, perturbation_key: str) -> dict:
        """
        Analyze perturbation structure and provide comprehensive intervention statistics.

        This utility function examines the distribution and complexity of perturbations
        in single-cell intervention datasets, providing essential quality control metrics
        for SENA model configuration and experimental design validation.

        Parameters
        ----------
        adata : AnnData
            Single-cell perturbation dataset with intervention annotations.
            Must contain perturbation metadata in .obs[perturbation_key].
        perturbation_key : str
            Column name in adata.obs containing perturbation identifiers.
            Expected format: "", "gene1", "gene1,gene2" for controls and perturbations.

        Returns
        -------
        dict
            Comprehensive perturbation statistics including:
            - Cell counts by perturbation type (control, single, multiple)
            - Most frequent intervention targets
            - Unique gene counts and intervention complexity metrics
            - Data quality recommendations for SENA training

        Notes
        -----
        Statistical Analysis:
        - Control vs. perturbed cell counts and ratios
        - Single vs. multiple perturbation distribution
        - Most frequent intervention targets and combinations
        - Unique gene coverage and intervention complexity

        Quality Control Insights:
        - Sufficient control cell representation for robust modeling
        - Perturbation target diversity for comprehensive causal learning
        - Multi-gene intervention complexity for pathway interaction discovery
        - Data balance assessment for training stability

        Examples
        --------
        >>> # Analyze CRISPR screen structure
        >>> stats = SENADVAE.analyze_perturbations(adata, "perturbation")
        >>> print(f"Controls: {stats['n_controls']}, Perturbed: {stats['n_perturbed']}")
        >>> print(f"Unique targets: {stats['n_unique_genes']}")
        """
        if perturbation_key not in adata.obs.columns:
            raise ValueError(
                f"Perturbation key '{perturbation_key}' not found in adata.obs columns."
            )

        # Process perturbation annotations
        perturbations = adata.obs[perturbation_key].fillna("")

        # Categorize cells by perturbation type
        controls = perturbations == ""
        single_perturb = (perturbations != "") & (~perturbations.str.contains(","))
        multi_perturb = perturbations.str.contains(",")

        # Extract all unique intervention genes
        all_genes = set()
        for pert_str in perturbations[perturbations != ""]:
            all_genes.update([g.strip() for g in pert_str.split(",")])

        # Compile comprehensive perturbation statistics
        summary = {
            "n_total_cells": len(adata),
            "n_controls": controls.sum(),
            "n_single_perturbations": single_perturb.sum(),
            "n_multi_perturbations": multi_perturb.sum(),
            "n_unique_combinations": len(perturbations[perturbations != ""].unique()),
            "n_unique_genes": len(all_genes),
            "most_common_perturbations": perturbations.value_counts().head(10).to_dict(),
            "intervention_genes": sorted(all_genes),
            "control_ratio": controls.sum() / len(adata),
            "perturbation_complexity": multi_perturb.sum()
            / (single_perturb.sum() + multi_perturb.sum())
            if (single_perturb.sum() + multi_perturb.sum()) > 0
            else 0,
        }

        return summary

    def train(
        self,
        max_epochs: int = 100,
        alpha_max: float = 1.0,
        beta_max: float = 1.0,
        temp_max: float = 100.0,
        alpha_start_epoch: int = 5,
        beta_start_epoch: int = 10,
        temp_start_epoch: int = 5,
        accelerator: str = "auto",
        devices: int | str = "auto",
        train_size: float = 0.9,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        batch_size: int = 32,
        early_stopping: bool = False,
        check_val_every_n_epoch: int = 1,
        training_plan: TrainingPlan
        | None = None,  # Remember that this one is not used, as we have a custom training plan
        plan_kwargs: dict | None = None,
        **trainer_kwargs,
    ):
        """
        Train the SENA model with scheduled loss weights and temperature annealing.

        Parameters
        ----------
        max_epochs : int, default 100
            Maximum number of training epochs.
        alpha_max : float, default 1.0
            Maximum weight for intervention loss.
        beta_max : float, default 1.0
            Maximum weight for KL divergence loss.
        temp_max : float, default 100.0
            Maximum temperature for intervention softmax.
        alpha_start_epoch : int, default 5
            Epoch to start ramping up alpha.
        beta_start_epoch : int, default 10
            Epoch to start ramping up beta.
        temp_start_epoch : int, default 5
            Epoch to start ramping up temperature.
        accelerator : str, default "auto"
            Accelerator type for training.
        devices : int or str, default "auto"
            Number of devices to use.
        train_size : float, default 0.9
            Proportion of data to use for training.
        validation_size : float, optional
            Proportion of data to use for validation.
        shuffle_set_split : bool, default True
            Whether to shuffle data when splitting.
        batch_size : int, default 32
            Batch size for training.
        early_stopping : bool, default False
            Whether to use early stopping.
        check_val_every_n_epoch : int, default 1
            How often to run validation. Set to 1 to validate every epoch.
        training_plan : TrainingPlan, optional
            Training plan to use (not used as we have custom training plan).
        plan_kwargs : dict, optional
            Additional keyword arguments for the training plan.
        **trainer_kwargs
            Additional keyword arguments for the trainer.
        """
        # Create custom SENA dataloader for reference (callbacks may need it)
        train_dataloader = SENADataLoader(
            adata_manager=self.adata_manager,
            perturbation_key=self.perturbation_key,
            batch_size=batch_size,
            shuffle=True,  # Always shuffle during training
            mode="training",  # Training mode - only single perturbations
        )

        # Store reference for callbacks
        self._current_dataloader = train_dataloader

        # Create schedulers
        # Those schedulers are custom scheduleres that inherit from Callback class

        # This one is form the MSE/MMD and KL divergence regularization parameters (alpha and beta)
        loss_scheduler = LossWeightScheduler(
            alpha_max=alpha_max,
            beta_max=beta_max,
            alpha_start_epoch=alpha_start_epoch,
            beta_start_epoch=beta_start_epoch,
        )

        # This one is for the temp parameter that is used to determin the sharpness
        # the softmax function
        temp_scheduler = TemperatureScheduler(temp_max=temp_max, temp_start_epoch=temp_start_epoch)

        # Create control reshuffle callback for per-epoch control resampling
        control_reshuffle_callback = ControlReshuffleCallback()

        # Create custom batch progress bar for batch-level progress tracking
        batch_progress_bar = SENABatchProgressBar()

        # Import ModelCheckpoint for best model saving
        from lightning.pytorch.callbacks import ModelCheckpoint

        # Create ModelCheckpoint callback to save best model based on validation loss
        best_model_checkpoint = ModelCheckpoint(
            monitor="validation_loss",  # Monitor total validation loss (all components combined)
            mode="min",  # Save model with minimum validation loss
            save_top_k=1,  # Keep only the best model
            save_last=False,  # Don't save the last model automatically
            verbose=False,  # Reduce checkpoint messages
            filename="best_model",  # Name for the best model checkpoint
        )

        # Add schedulers, control reshuffle callback, custom progressbar and best model checkpoint
        callbacks = [
            loss_scheduler,
            temp_scheduler,
            control_reshuffle_callback,
            batch_progress_bar,
            best_model_checkpoint,  # Add best model checkpoint callback
        ]
        if "callbacks" in trainer_kwargs:
            trainer_kwargs["callbacks"].extend(callbacks)
        else:
            trainer_kwargs["callbacks"] = callbacks

        # Setup training plan - don't create it ourselves, let the mixin handle it
        # but pass our custom plan_kwargs if needed
        if plan_kwargs is None:
            plan_kwargs = {}

        # Remove training_plan from trainer_kwargs if it exists to avoid duplicate argument
        if "training_plan" in trainer_kwargs:
            trainer_kwargs.pop("training_plan")

        # Ensure perturbation_key is passed to the data splitter
        datasplitter_kwargs = trainer_kwargs.get("datasplitter_kwargs", {})
        datasplitter_kwargs["perturbation_key"] = self.perturbation_key
        trainer_kwargs["datasplitter_kwargs"] = datasplitter_kwargs

        # Configure trainer for batch-level progress display
        trainer_kwargs.update(
            {
                "log_every_n_steps": 1,  # Log metrics every batch for real-time updates
                "enable_progress_bar": False,  # Disable default progress bars (custom one)
                "simple_progress_bar": False,  # Disable scvi-tools progress bar (custom one)
                "enable_model_summary": False,  # Reduce startup messages
                "enable_checkpointing": True,  # Enable checkpointing for best model saving
            }
        )

        # Reduce logging verbosity for cleaner output
        import logging

        logging.getLogger("scvi").setLevel(logging.WARNING)
        logging.getLogger("lightning").setLevel(logging.WARNING)
        logging.getLogger("pytorch_lightning").setLevel(logging.WARNING)

        # Pass the arguments to the parent train method
        # The UnsupervisedTrainingMixin will use our SENADataSplitter automatically
        super().train(
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=devices,
            train_size=train_size,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            batch_size=batch_size,
            early_stopping=early_stopping,
            plan_kwargs=plan_kwargs,  # Let mixin create training_plan with our kwargs
            check_val_every_n_epoch=check_val_every_n_epoch,  # Enable validation
            **trainer_kwargs,  # Includes datasplitter_kwargs with perturbation_key
        )

        # Load the best model checkpoint after training completes
        self._load_best_checkpoint()

    def _load_best_checkpoint(self):
        """Load the best model checkpoint after training completes."""
        try:
            # Find the best model checkpoint in the trainer's checkpoint directory
            if (
                hasattr(self.trainer, "checkpoint_callback")
                and self.trainer.checkpoint_callback is not None
            ):
                best_model_path = self.trainer.checkpoint_callback.best_model_path
                if best_model_path and os.path.exists(best_model_path):
                    # Load the best model state
                    checkpoint = torch.load(best_model_path, map_location=self.device)
                    self.module.load_state_dict(checkpoint["state_dict"])
                    logger.info(f"Loaded best model from checkpoint: {best_model_path}")
                else:
                    logger.warning("Best model checkpoint not found, using final epoch model")
            else:
                logger.warning("No checkpoint callback found, using final epoch model")
        except (ValueError, RuntimeError, KeyError, TypeError) as e:
            logger.warning(f"Failed to load best checkpoint: {e}, using final epoch model")

    def get_latent_representation(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        give_mean: bool = True,
        batch_size: int | None = None,
    ) -> np.ndarray:
        """
        Return the latent representation for each cell.

        Parameters
        ----------
        adata : AnnData, optional
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices : Sequence[int], optional
            Indices of cells in adata to use. If `None`, all cells are used.
        give_mean : bool, default True
            Give mean of distribution or sample from it.
        batch_size : int, optional
            Batch size for data loading. If `None`, full data is used.

        Returns
        -------
        np.ndarray
            Latent representation of cells.
        """
        self._check_if_trained(warn=False)
        # checks if the andata has been resisteger with th elatest andata manager.
        # It is a basemoduleclasss methos that inturn has functions from AnnDataManager()
        adata = self._validate_anndata(adata)
        # From the indices of cells given it extracts the latent representation.
        # For sme reason it does so in batches, not all at once
        # This just create an instance of the custom dataloader
        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        latent = []

        # Then it iterates over the tensors in the trained andata object,
        # extracts the input of the gene expresion and performs the mu,sd and z calculations
        for tensors in scdl:
            inference_inputs = self.module._get_inference_input(tensors)
            outputs = self.module.inference(**inference_inputs)
            if give_mean:
                latent_sample = outputs["z_mu"]
            else:
                latent_sample = outputs["z"]
            # detach() isoaltes the tensor from the gradients computation graph and
            # cpu moves the tensor from the gpu to the cpu
            latent += [latent_sample.detach().cpu()]
        # It returns the corresponding tensor (mu or z)
        return torch.cat(latent).numpy()

    def get_metric(self, adata, batch_size, N=100, temp_values=None) -> dict:
        """
        Calculate comprehensive pathway activity analysis and perturbation metrics.

        This method provides a complete analysis similar to the generating_data function,
        extracting all neural network components, pathway activities, intervention effects,
        and traditional metrics for comprehensive perturbation analysis.

        Parameters
        ----------
        adata : AnnData
            AnnData object containing perturbation data.
        batch_size : int
            Batch size for processing.
        N : int, default 100
            Top N pathways to consider for hit ratio calculation.
        temp_values : list, default [1, 100, 1000]
            Temperature values for intervention effect analysis.

        Returns
        -------
        dict
            Comprehensive results dictionary containing:
            - Traditional metrics: da_p, dar_p, hitn
            - Neural network components: fc1, fc_mean, fc_var, z
            - Intervention effects: bc_temp1, bc_temp100, bc_temp1000
            - Causal analysis: u, causal_graph
            - Network weights: mean_delta_matrix, std_delta_matrix
            - Mappings: pert_map
            - GO pathways: gos (list of pathway identifiers)

        Notes
        -----
        This method provides complete access to SENA's internal representations
        and learned parameters, enabling detailed analysis of pathway activities,
        causal relationships, and intervention mechanisms.
        """
        logger.info("Starting comprehensive pathway activity analysis.")

        # Initialize default temperature values if not provided
        if temp_values is None:
            temp_values = [1, 100, 1000]

        # Initialize results dictionary with all components
        results_dict = {
            # Traditional perturbation metrics
            "da_p": {},
            "dar_p": {},
            "hitn": {},
            # Neural network pathway activities
            "fc1": {},  # Raw pathway scores from encoder
            "fc_mean": {},  # VAE mean parameters
            "fc_var": {},  # VAE variance parameters
            "z": {},
            "z_interv": {},
            "z_interv_temp_1": {},
            "z_interv_temp_100": {},
            "z_interv_temp_1000": {},  # Sampled latents (reparameterized)
            "activity_score": {},  # Pathway activity scores for traditional metrics
            # Intervention effects at different temperatures
            "bc_temp1": {},  # Temperature=1 intervention gates
            "bc_temp100": {},  # Temperature=100 intervention gates
            "bc_temp1000": {},  # Temperature=1000 intervention gates
            # Final pathway activities after causal propagation
            "u": {},
            # Causal graph and network weights
            "causal_graph": None,
            "mean_delta_matrix": None,
            "std_delta_matrix": None,
            "pert_map": {},
            # GO pathway information
            "gos": self.gos,  # List of GO pathway identifiers
        }

        # Get indices for control and perturbed cells
        indices = np.where(adata.obs[self.perturbation_key] != "")[0]
        indices_control = np.where(adata.obs[self.perturbation_key] == "")[0]

        # Create dataloaders for perturbed and control cells
        perturbed_dataloader = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )
        controls_dataloader = self._make_data_loader(
            adata=adata, indices=indices_control, batch_size=batch_size
        )

        # Extract neural network weights
        self._extract_network_weights(results_dict)

        # Process perturbed cells
        self._process_perturbed_cells(perturbed_dataloader, results_dict, temp_values)

        # Process control cells
        self._process_control_cells(controls_dataloader, results_dict)

        # Calculate traditional metrics (da_p, dar_p, hitn)

        # Convert pert_map to DataFrame format (equivalent to original code)

        # Convert network weights to DataFrames with GO pathway indices
        if "mean_delta_matrix" in results_dict and results_dict["mean_delta_matrix"] is not None:
            results_dict["mean_delta_matrix"] = pd.DataFrame(
                results_dict["mean_delta_matrix"].T, index=self.gos
            )

        if "std_delta_matrix" in results_dict and results_dict["std_delta_matrix"] is not None:
            results_dict["std_delta_matrix"] = pd.DataFrame(
                results_dict["std_delta_matrix"].T, index=self.gos
            )

        # Convert stored activation lists to pandas DataFrames (equivalent to original processing)
        logger.info("Converting activation data to pandas DataFrames.")

        # Process each layer's data into DataFrames (including activity_score_fdr)
        layers_to_process = [
            "fc1",
            "fc_mean",
            "fc_var",
            "z",
            "z_interv",
            "u",
            "bc_temp1",
            "bc_temp100",
            "bc_temp1000",
            "activity_score",
            "z_interv_temp1",
            "z_interv_temp100",
            "z_interv_temp1000",
        ]
        for layer in layers_to_process:
            if layer in results_dict and results_dict[layer] is not None:
                # Check if it's already a DataFrame (for activity_score_fdr)
                if isinstance(results_dict[layer], pd.DataFrame):
                    continue

                temp_df = []
                for gene in results_dict[layer]:
                    gene_data_list = results_dict[layer][gene]
                    if len(gene_data_list) > 0:  # Check if list is not empty
                        # Convert list of arrays to DataFrame
                        gene_data = pd.DataFrame(np.vstack(gene_data_list))
                        gene_data.index = [gene] * gene_data.shape[0]
                        temp_df.append(gene_data)

                if temp_df:  # Only concatenate if we have data
                    # Substitute the dictionary with concatenated DataFrame
                    results_dict[layer] = pd.concat(temp_df)

                    # Add GO pathway column names for layers with pathway dimensions
                    if layer in ["fc1", "activity_score", "activity_score_fdr"]:
                        results_dict[layer].columns = self.gos
        self._calculate_traditional_metrics(results_dict, N)

        # Perform statistical testing on activity scores after DataFrame conversion
        if "activity_score" in results_dict and isinstance(
            results_dict["activity_score"], pd.DataFrame
        ):
            try:
                logger.info("Performing statistical testing on pathway activity scores...")
                fdr_results = self._test_activity_score_significance(
                    results_dict["activity_score"]
                )
                results_dict["activity_score_fdr"] = fdr_results
                logger.info("Statistical testing completed successfully.")
            except (ValueError, RuntimeError, KeyError, TypeError) as e:
                logger.warning(f"Statistical testing failed: {e}")
                results_dict["activity_score_fdr"] = None

        logger.info("Comprehensive analysis completed successfully.")
        return results_dict

    def _extract_network_weights(self, results_dict):
        """Extract and store neural network weights."""
        try:
            # Extract encoder weights
            if hasattr(self.module.fc_mean, "weight"):
                weight_tensor = self.module.fc_mean.weight.detach().cpu().numpy()
                results_dict["mean_delta_matrix"] = weight_tensor
            if hasattr(self.module.fc_var, "weight"):
                weight_tensor = self.module.fc_var.weight.detach().cpu().numpy()
                results_dict["std_delta_matrix"] = weight_tensor

            # Extract causal graph
            if hasattr(self.module, "G"):
                results_dict["causal_graph"] = self.module.G.detach().cpu().numpy()

        except (AttributeError, RuntimeError) as e:
            logger.warning(f"Could not extract network weights: {e}")

    def _process_perturbed_cells(self, dataloader, results_dict, temp_values):
        """Process perturbed cells to extract all neural network activations."""
        logger.info("Processing perturbed cells for comprehensive analysis.")

        with torch.no_grad():
            for tensors in dataloader:
                # Get intervention information (one-hot encoded)
                perturbations = tensors["extra_categorical_covs"]

                # Standard inference pass
                inference_inputs = self.module._get_inference_input(tensors, metrics=True)
                # inference_outputs has: "z": z, "z_mu": z_mu, "z_var": z_var, "activity_score": h
                inference_outputs = self.module.inference(**inference_inputs)

                # Extract pathway activities and latent representations
                # gene_expression
                x = inference_inputs["x"]

                # Get fc1 outputs (raw pathway activities) we need the softmax also
                fc1_output = self.module.fc1(x)

                # Get VAE parameters
                fc_mean_output = self.module.fc_mean(fc1_output)
                fc_var_output = F.softplus(self.module.fc_var(fc1_output))

                # Get sampled latents
                z = inference_outputs["z"]

                # Process interventions at different temperatures
                generative_inputs = self.module._get_generative_input(tensors, inference_outputs)

                # Get internal variables at different temperatures including bc values
                z_interv_by_temp = {}
                bc1_by_temp = {}
                bc2_by_temp = {}

                for t in temp_values:
                    generative_outputs = self.module.generative(**generative_inputs, temp=t)
                    z_interv_by_temp[f"z_interv_temp{t}"] = generative_outputs["z_interv"]
                    bc1_by_temp[f"bc_temp{t}"] = generative_outputs["bc1"]
                    bc2_by_temp[f"bc_temp{t}"] = generative_outputs["bc2"]

                # Use default temp for main variables (maintaining backward compatibility)
                generative_outputs = self.module.generative(**generative_inputs, temp=1.0)

                # Extract all internal variables from generative outputs
                u = generative_outputs["u"]
                z_interv = generative_outputs["z_interv"]

                # Store all activations grouped by perturbation including bc values
                activations_dict = {
                    "fc1": fc1_output,
                    "fc_mean": fc_mean_output,
                    "fc_var": fc_var_output,
                    "z": z,
                    "z_interv": z_interv,
                    "u": u,
                    "activity_score": inference_outputs["activity_score"],
                }

                # Add temperature-specific z_interv
                for temp_key, z_interv_temp in z_interv_by_temp.items():
                    activations_dict[temp_key] = z_interv_temp

                # Add temperature-specific bc values
                for bc_key, bc_temp in bc1_by_temp.items():
                    if bc_temp is not None:
                        activations_dict[bc_key] = bc_temp

                # Store all activations grouped by perturbation
                self._store_activations_by_perturbation(
                    perturbations,
                    dataloader,
                    activations_dict,
                    results_dict,
                )

    def _store_activations_by_perturbation(
        self, perturbations, dataloader, activations, results_dict
    ):
        """Store neural network activations grouped by perturbation type."""
        # Get perturbation strings from dataloader
        onehot_gene_dict = dataloader.gene_to_intervention_idx

        # Convert one-hot to perturbation strings
        perturbation_strings = []
        for _, pert_vector in enumerate(perturbations):
            genes = []
            for j, val in enumerate(pert_vector):
                if val == 1:
                    # Find gene name from index
                    for gene_name, gene_idx in onehot_gene_dict.items():
                        if gene_idx == j:
                            genes.append(gene_name)
                            break

            if genes:
                pert_string = ",".join(genes)
            else:
                pert_string = "ctrl"
            perturbation_strings.append(pert_string)

        # Group activations by perturbation
        for i, pert_string in enumerate(perturbation_strings):
            for key, tensor in activations.items():
                if key not in results_dict:
                    continue

                if pert_string not in results_dict[key]:
                    results_dict[key][pert_string] = []

                # Store single cell activation
                cell_activation = tensor[i : i + 1].detach().cpu().numpy()
                results_dict[key][pert_string].append(cell_activation)

    def _process_control_cells(self, dataloader, results_dict):
        """Process control cells for baseline measurements."""
        logger.info("Processing control cells for baseline analysis.")

        with torch.no_grad():
            for tensors in dataloader:
                inference_inputs = self.module._get_inference_input(tensors, metrics=True)
                inference_outputs = self.module.inference(**inference_inputs)

                # Extract control cell activations
                x = inference_inputs["x"]
                fc1_output = self.module.fc1(x)
                fc_mean_output = self.module.fc_mean(fc1_output)
                fc_var_output = F.softplus(self.module.fc_var(fc1_output))
                z = inference_outputs["z"]

                # Calculate u for controls (no intervention)
                I = torch.eye(self.module.n_latent, device=z.device)
                dag_matrix = torch.inverse(I - torch.triu(self.module.G, diagonal=1))
                u = z @ dag_matrix

                # Store control activations
                control_activations = {
                    "fc1": fc1_output,
                    "fc_mean": fc_mean_output,
                    "fc_var": fc_var_output,
                    "z": z,
                    "z_interv": z,  # For controls, z_interv = z (no intervention)
                    "u": u,
                    "activity_score": inference_outputs["activity_score"],
                }

                for key, tensor in control_activations.items():
                    if key not in results_dict:
                        continue

                    if "ctrl" not in results_dict[key]:
                        results_dict[key]["ctrl"] = []

                    results_dict[key]["ctrl"].append(tensor.detach().cpu().numpy())

    def _calculate_traditional_metrics(self, results_dict, N):
        """Calculate traditional perturbation metrics (da_p, dar_p, hitn)."""
        logger.info("Calculating traditional perturbation metrics.")

        # Use activity_score for traditional metrics calculation
        if "activity_score" not in results_dict:
            logger.warning("No activity scores found for traditional metrics calculation.")
            return

        # Get control baseline
        if "ctrl" not in results_dict["activity_score"].index:
            logger.warning("No control cells found for baseline calculation.")
            return

        ctrl_activities = results_dict["activity_score"].loc["ctrl"].values
        ctrl_mean = np.mean(ctrl_activities, axis=0)

        # Initialize empty DataFrame for DA values with NaN values
        df_da = pd.DataFrame(
            data=np.nan,
            index=list(np.unique(results_dict["activity_score"].index)),
            columns=self.gos,
        )

        # Calculate metrics for each perturbation
        for pert_strings in np.unique(results_dict["activity_score"].index):
            if pert_strings == "ctrl":
                continue
            activities = results_dict["activity_score"].loc[pert_strings].values
            pert_mean = np.mean(activities, axis=0)

            # Calculate DA (Differential Activity) as absolute difference per paper definition
            # DA^p_k = |α̅^p_k - α̅^c_k| where α̅^p_k and α̅^c_k are mean activities
            raw_diff = pert_mean - ctrl_mean
            da_values = np.abs(raw_diff)
            df_da.loc[pert_strings] = da_values

            # Calculate DAR and HitN using gene-pathway mappings
        self._calculate_dar_hitn(df_da, results_dict, N)
        results_dict["da_p"] = df_da

    def _calculate_dar_hitn(self, df_da, results_dict, N):
        """Calculate DAR and HitN metrics for a specific perturbation."""
        df_dar = pd.DataFrame(
            data=np.nan, index=list(results_dict["activity_score"].index), columns=["dar"]
        )
        df_hitn = pd.DataFrame(
            data=np.nan, index=list(results_dict["activity_score"].index), columns=["hitn"]
        )

        for pert_string in df_da.index:
            pathway_indices = []
            # get the pathways where this pertubration is included
            gene_symbols = pert_string.split(",")
            for gene_symbol in gene_symbols:
                try:
                    # Convert gene symbol to ensembl ID
                    if hasattr(self, "mapping_dict") and gene_symbol in self.mapping_dict:
                        ensembl_id = self.mapping_dict[gene_symbol]

                        # Get gene index
                        if hasattr(self, "gen_dict") and ensembl_id in self.gen_dict:
                            gene_idx = self.gen_dict[ensembl_id]

                            # Get pathway indices for this gene
                            if hasattr(self, "rel") and gene_idx in self.rel:
                                pathway_indices.extend(self.rel[gene_idx])

                except KeyError:
                    logger.warning(f"Could not map gene {gene_symbol} to pathways")
                    continue
            not_pathway_indices = [
                i for i in range(len(df_da.columns)) if i not in pathway_indices
            ]
            # add a print of the head of df_da

            wp = df_da.loc[pert_string].values[pathway_indices]
            wp_n = df_da.loc[pert_string].values[not_pathway_indices]
            dar = np.mean(wp_n) / np.mean(wp)

            # Now wwe will need to sort the DA values and see how many of the pathways
            # that include the perturbed gene are in the top N
            # quiza habria que añadir una columna al laod con los valores del index y
            # luego ordenarlo por la columna de da
            # Crea un dataframe con los valores de da y su indice
            da = df_da.loc[pert_string].values
            inde = np.arange(0, len(df_da.columns))
            print(df_da.index)
            print(inde)
            print(da)
            print(f"da shape: {da.shape}")
            df = pd.DataFrame({"da": da, "index": inde})
            # now sort the df by the column da
            df = df.sort_values(by="da", ascending=False)
            # now select the number of rows with index in the wp list that are within
            # the first 100 rows
            hits = df.loc[:, "index"].isin(wp)
            # now lets get the numbers of trues in the first N rows
            hitN = sum(hits[:N]) / N

            df_dar.loc["dar"] = dar
            df_hitn.loc["hitn"] = hitN
        results_dict["dar_p"] = df_dar
        results_dict["hitn"] = df_hitn

    def predict_perturbation_response(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
        temp: float = 1.0,
    ) -> np.ndarray:
        """
        Predict gene expression response to perturbations.

        Note: num_interv is now auto-detected based on the intervention matrix.

        Parameters
        ----------
        adata : AnnData, optional
            AnnData object with equivalent structure to initial AnnData.
        indices : Sequence[int], optional
            Indices of cells in adata to use.
        batch_size : int, optional
            Batch size for data loading.
        temp : float, default 1.0
            Temperature for intervention softmax.

        Returns
        -------
        np.ndarray
            Predicted perturbed gene expression.
        """
        # Once the model has been trained with single perturbation we want to know
        # the unknown perturbations
        self._check_if_trained(warn=False)
        adata = self._validate_anndata(adata)
        # Load the
        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        predictions = []

        # In this case it computes the inference and generative parts to obtain all the
        # elements of the nerual network. Then it selects the output predicted expresion

        for tensors in scdl:
            inference_inputs = self.module._get_inference_input(tensors)
            inference_outputs = self.module.inference(**inference_inputs)

            # This returns the one hot vector c1 and c2 along with th z latent vector and
            # auto-detected num_interv
            generative_inputs = self.module._get_generative_input(tensors, inference_outputs)

            # This calculates the shift, computes the activty after the acicclic graph operations
            # and reconstruct both the control and the predicted gene expresion
            # Use auto-detected num_interv from generative_inputs
            generative_outputs = self.module.generative(
                z=generative_inputs["z"],
                c1=generative_inputs["c1"],
                c2=generative_inputs["c2"],
                num_interv=generative_inputs["num_interv"],  # Auto-detected
                temp=temp,
            )
            # It then returns te predicted pertubration gene expresion
            predictions += [generative_outputs["y_hat"].detach().cpu()]
        # This returns the predicted gene exression
        return torch.cat(predictions).numpy()

    def get_causal_graph(self) -> np.ndarray:
        """
        Get the learned causal graph matrix.

        Returns
        -------
        np.ndarray
            Causal graph adjacency matrix.
        """
        self._check_if_trained(warn=False)
        return torch.triu(self.module.G, diagonal=1).detach().cpu().numpy()

    def _test_activity_score_significance(
        self, activity_score_df: pd.DataFrame, alpha: float = 0.05, method: str = "fdr_bh"
    ) -> pd.DataFrame:
        """
        Perform two-tailed t-tests comparing pathway activity scores between perturbations & ctrls.

        For each GO term, compares the activity scores of each perturbation condition against
        control cells, then applies FDR correction across all tests.

        Parameters
        ----------
        activity_score_df : pd.DataFrame
            DataFrame with perturbation conditions as rows and GO terms as columns.
            Must include 'ctrl' rows for control comparisons.
        alpha : float, default 0.05
            Significance level for FDR correction.
        method : str, default 'fdr_bh'
            Method for multiple testing correction.

        Returns
        -------
        pd.DataFrame
            DataFrame with perturbations as rows and GO terms as columns containing FDR values.
        """
        # Validate input
        if not isinstance(activity_score_df, pd.DataFrame):
            raise ValueError("activity_score_df must be a pandas DataFrame")

        # Extract control data
        ctrl_mask = activity_score_df.index == "ctrl"
        if not ctrl_mask.any():
            raise ValueError("No control cells found. Expected 'ctrl' in DataFrame index.")

        ctrl_data = activity_score_df[ctrl_mask]

        # Get unique perturbation conditions (excluding control)
        perturbation_conditions = activity_score_df.index[~ctrl_mask].unique()

        # Initialize results storage
        fdr_results = {}

        for condition in perturbation_conditions:
            # Get perturbation data for this condition
            pert_mask = activity_score_df.index == condition
            pert_data = activity_score_df[pert_mask]

            # Skip if insufficient data
            if len(pert_data) < 2:
                logger.warning(f"Skipping {condition}: insufficient cells (n={len(pert_data)})")
                continue

            # Perform t-tests for each GO term
            p_values = []
            tested_gos = []

            for go_term in activity_score_df.columns:
                try:
                    # Get values for this GO term
                    ctrl_values = ctrl_data[go_term].values
                    pert_values = pert_data[go_term].values

                    # Remove any NaN values
                    ctrl_values = ctrl_values[~np.isnan(ctrl_values)]
                    pert_values = pert_values[~np.isnan(pert_values)]

                    # Check if we have sufficient data
                    if len(ctrl_values) < 2 or len(pert_values) < 2:
                        continue

                    # Perform two-tailed t-test
                    t_stat, p_val = stats.ttest_ind(
                        pert_values,
                        ctrl_values,
                        equal_var=False,  # Welch's t-test (unequal variances)
                    )

                    p_values.append(p_val)
                    tested_gos.append(go_term)

                except (ValueError, RuntimeError, KeyError, TypeError) as e:
                    logger.warning(f"Error testing {condition} vs ctrl for {go_term}: {e}")
                    continue

            # Apply FDR correction if we have tests
            if len(p_values) > 0:
                # Perform multiple testing correction
                rejected, p_adjusted, alpha_sidak, alpha_bonf = multipletests(
                    p_values, alpha=alpha, method=method
                )

                # Store results for this condition
                condition_results = {}
                for go_term, fdr_val in zip(tested_gos, p_adjusted, strict=True):
                    condition_results[go_term] = fdr_val

                # Fill in NaN for GO terms that weren't tested
                for go_term in activity_score_df.columns:
                    if go_term not in condition_results:
                        condition_results[go_term] = np.nan

                fdr_results[condition] = condition_results
            else:
                logger.warning(f"No valid tests for condition {condition}")
                # Fill with NaN if no tests were possible
                fdr_results[condition] = dict.fromkeys(activity_score_df.columns, np.nan)

        # Convert to DataFrame
        fdr_df = pd.DataFrame(fdr_results).T

        # Ensure column order matches input
        fdr_df = fdr_df.reindex(columns=activity_score_df.columns)

        return fdr_df

    def _make_data_loader(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
        shuffle: bool = False,
        data_loader_class=None,
        **data_loader_kwargs,
    ):
        """
        Create a custom SENA data loader that handles control-perturbation matching.

        This method creates a SENADataLoader that performs the exact same preprocessing
        as the original setup_anndata method, but applies it on-the-fly during training.

        Parameters
        ----------
        adata : AnnData, optional
            AnnData object to create dataloader for.
        indices : Sequence[int], optional
            Indices of cells to include.
        batch_size : int, optional
            Batch size for data loading.
        shuffle : bool, default False
            Whether to shuffle the data.
        data_loader_class : optional
            Data loader class to use (ignored, uses SENADataLoader).
        **data_loader_kwargs
            Additional arguments for data loader.

        Returns
        -------
        SENADataLoader
            Custom SENA data loader instance.
        """
        if batch_size is None:
            batch_size = 128

        if adata is None:
            adata_manager = self.adata_manager
        else:
            # Validate that the new adata has the same structure
            adata_manager = self._validate_anndata(adata)
            adata_manager = self.adata_manager

        # Create SENA-specific dataloader that replicates the original preprocessing
        # Use prediction mode if indices are specified, training mode otherwise
        mode = "prediction" if indices is not None else "training"

        sena_dataloader = SENADataLoader(
            adata_manager=adata_manager,
            perturbation_key=self.perturbation_key,
            shuffle=shuffle,
            batch_size=batch_size,
            mode=mode,
            indices=indices,  # Pass indices for prediction mode
            **data_loader_kwargs,
        )

        # Store reference to current dataloader for control reshuffling
        self._current_dataloader = sena_dataloader

        return sena_dataloader
