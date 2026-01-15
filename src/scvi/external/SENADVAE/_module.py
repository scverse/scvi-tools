import logging

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from scvi import REGISTRY_KEYS
from scvi.module.base import BaseModuleClass, auto_move_data

logger = logging.getLogger(__name__)


class NetworkActivity_layer(torch.nn.Module):
    """
    Biologically-constrained neural network layer implementing Gene Ontology pathway constraints.

    This layer represents the core biological innovation of SENA, transforming raw gene expression
    data into biologically meaningful pathway activity scores. Unlike standard neural network
    layers that learn arbitrary gene-gene relationships, this layer enforces
    biological priors by only allowing connections between genes and their experimentally-validated
    Gene Ontology annotations.

    Parameters
    ----------
    input_genes : int
        Number of input genes in expression profiles (typically 5,000-20,000).
    output_gs : int
        Number of output Gene Ontology pathways (typically 100-500 after filtering).
    relation_dict : dict
        Gene-to-pathway mapping from GO annotations.
        Format: {gene_index: [pathway_index1, pathway_index2, ...]}.
    bias : bool, default True
        Whether to include learnable bias terms for each pathway.
    lambda_parameter : float, default 0
        Regularization allowing minimal contribution from non-annotated gene-pathway pairs.
        Prevents complete information loss (typically 0.01-0.1).

    Attributes
    ----------
    n_input_genes : int
        Number of genes in expression profile.
    n_output_gs : int
        Number of GO pathways after filtering.
    relation_dict : dict
        Gene index to list of pathway indices mapping.
    mask : torch.Tensor
        Biological constraint mask enforcing GO annotations.
    weight : nn.Parameter
        Learnable weights for gene contributions to pathway activities.
    bias : nn.Parameter or None
        Learnable bias term for each pathway baseline activity.

    Notes
    -----
    Biological Rationale:
    - Genes function within coordinated biological pathways/processes
    - Gene Ontology provides curated, experimental evidence-based gene-pathway annotations
    - Pathway-level analysis reduces noise and improves interpretability
    - Enforces sparsity consistent with biological network structure

    Mathematical Implementation:
    Creates a masked weight matrix where W[pathway_i, gene_j] = 0 unless gene_j is annotated
    to pathway_i in Gene Ontology. The lambda_parameter allows minimal "leakage" to prevent
    complete information loss for genes with sparse annotations.
    """

    def __init__(self, input_genes, output_gs, relation_dict, bias=True, lambda_parameter=0):
        super().__init__()
        self.n_input_genes = input_genes  # Number of genes in expression profile
        self.n_output_gs = output_gs  # Number of GO pathways after filtering
        self.relation_dict = relation_dict  # Gene index -> list of pathway indices mapping

        # Create biological constraint mask enforcing GO annotations
        # Only genes annotated to a pathway can contribute to its activity
        mask = torch.zeros((self.n_input_genes, self.n_output_gs))
        for gene_idx in range(self.n_input_genes):
            if gene_idx in self.relation_dict:
                for pathway_idx in self.relation_dict[gene_idx]:
                    mask[gene_idx, pathway_idx] = 1.0

        self.mask = mask
        # Apply lambda regularization to prevent complete information loss
        # Allows minimal contribution from non-annotated gene-pathway pairs
        self.mask[self.mask == 0] = lambda_parameter

        # Learnable weights for gene contributions to pathway activities
        # Shape: (n_pathways, n_genes)
        self.weight = nn.Parameter(torch.empty((self.n_output_gs, self.n_input_genes)))

        # Learnable bias term for each pathway baseline activity
        if bias:
            self.bias = nn.Parameter(torch.empty(self.n_output_gs))
        else:
            self.register_parameter("bias", None)
        self.reset_parameters()

    def forward(self, x):
        """
        Transform gene expression to pathway activities using biological constraints.

        Computes pathway activities as a biologically-constrained linear transformation:
        pathway_activity = (gene_expression @ masked_weights) + bias

        The masked weights ensure only biologically relevant gene-pathway connections
        contribute to the pathway activity computation.

        Parameters
        ----------
        x : torch.Tensor
            Gene expression matrix, shape (batch_size, n_genes).
            Values typically log-normalized single-cell counts.

        Returns
        -------
        torch.Tensor
            Pathway activity matrix, shape (batch_size, n_pathways).
            Each element represents inferred activity level of a biological pathway.
        """
        device = self.weight.device
        # Apply biological masking: only annotated gene-pathway pairs contribute
        # Transpose mask to match matrix multiplication requirements
        masked_weights = self.weight * self.mask.T.to(device)

        # Compute pathway activities as weighted sum of annotated genes
        output = x @ masked_weights.T

        if self.bias is not None:
            return output + self.bias
        return output

    def reset_parameters(self) -> None:
        """
        Initialize network parameters using standard PyTorch initialization schemes.

        Uses Kaiming uniform initialization for weights and uniform initialization
        for bias terms based on fan-in calculations.
        """
        nn.init.kaiming_uniform_(self.weight, a=np.sqrt(5))
        if self.bias is not None:
            fan_in, _ = nn.init._calculate_fan_in_and_fan_out(self.weight)
            bound = 1 / np.sqrt(fan_in) if fan_in > 0 else 0
            nn.init.uniform_(self.bias, -bound, bound)


class SENAModule(BaseModuleClass):
    """
    SENA neural architecture for causal perturbation modeling.

    This module implements a specialized Variational Autoencoder designed for single-cell
    perturbation experiments. It combines biological pathway constraints with causal graph
    learning to model how genetic perturbations propagate through cellular networks to
    affect gene expression.

    Parameters
    ----------
    n_input : int
        Number of genes in expression profiles (typically 5,000-20,000).
    n_latent : int
        Latent space dimensionality (equals number of GO pathways).
    n_cat_covs : int
        Number of intervention categories (same as n_categories_interv).
    n_categories_interv : int
        Total number of unique genes that can be perturbed.
    gos : list
        Gene Ontology pathway identifiers (e.g., ['GO:0006915', 'GO:0008219']).
    rel_dict : dict
        Gene-to-pathway mapping from GO annotations {gene_idx: [pathway_idx1, pathway_idx2]}.
    n_hidden_encoder : int, default 512
        Hidden layer size for standard encoder (used when mode="normal").
    n_hidden_decoder : int, default 128
        Hidden layer size for gene expression reconstruction decoder.
    n_hidden_interv : int, default 256
        Hidden layer size for perturbation effect encoding network.
    mode : str, default "sena"
        Architecture mode: "sena" (pathway-constrained) or "normal" (standard VAE).
    sena_lambda : float, default 0.1
        Regularization strength for pathway sparsity constraints.

    Attributes
    ----------
    n_input : int
        Number of input genes in expression profile.
    n_latent : int
        Latent space dimensionality (equals number of pathways).
    n_cat_covs : int
        Number of intervention categories.
    n_categories_interv : int
        Number of unique perturbable genes.
    mode : str
        Architecture mode ("sena" or "normal").
    sena_lambda : float
        Regularization strength for pathway sparsity constraints.
    relations : dict
        Gene-to-pathway mapping dictionary.
    fc1 : NetworkActivity_layer or nn.Linear
        First layer of encoder (pathway-constrained or standard).
    fc_mean : nn.Linear
        Linear layer for latent mean parameters.
    fc_var : nn.Linear
        Linear layer for latent variance parameters.
    d1 : nn.Linear
        First decoder layer.
    d2 : nn.Linear
        Second decoder layer.
    G : nn.Parameter
        Causal graph matrix modeling pathway-pathway interactions.
    c1 : nn.Linear
        First intervention encoding layer.
    c2 : nn.Linear
        Second intervention encoding layer.
    c_shift : nn.Parameter
        Learnable intervention strength parameters.
    activity_score : torch.Tensor
        Pathway activity scores from the encoder.

    Notes
    -----
    Scientific Framework:
    1. Gene Expression → Pathway Activities (biologically constrained via GO annotations)
    2. Pathway Activities → Latent Space (VAE probabilistic encoding with μ, σ)
    3. Perturbation Effects → Pathway Shifts (learnable intervention-specific effects)
    4. Causal Propagation → Modified Activities (DAG matrix models pathway interactions)
    5. Modified Activities → Gene Expression (decoder reconstruction/prediction)

    Key Biological Insights:
    - Perturbations act primarily at pathway level, not individual genes
    - Pathways interact causally (upstream pathways affect downstream ones)
    - Intervention effects can be decomposed into pathway-specific shifts
    - Causal structure is learnable from observational + intervention data

    Applications:
    - Predict single/double gene knockout effects
    - Identify causal pathway relationships
    - Design optimal perturbation experiments
    - Understand drug/therapeutic mechanisms of action
    """

    def __init__(
        self,
        n_input: int,  # Number of input genes in expression profile
        n_latent: int,  # Latent space dimensionality (equals number of pathways)
        n_cat_covs: int,  # Number of intervention categories
        n_categories_interv: int,  # Number of unique perturbable genes
        gos: list,  # List of GO pathway IDs that passed filtering
        rel_dict: dict,  # Gene index -> pathway indices mapping
        n_hidden_encoder: int = 128,  # Standard encoder hidden layer size (mode="normal")
        n_hidden_decoder: int = 128,  # Decoder hidden layer size
        n_hidden_interv: int = 128,  # Intervention encoder hidden layer size
        mode: str = "sena",
        sena_lambda: float = 0,
    ):
        super().__init__()  # Initialize parent BaseModuleClass

        # Store architecture parameters
        self.n_input = n_input
        self.n_latent = n_latent
        self.n_cat_covs = n_cat_covs
        self.n_categories_interv = n_categories_interv
        self.mode = mode
        self.sena_lambda = sena_lambda
        self.relations = rel_dict

        # Encoder pathway: Gene Expression → Pathway Activities → Latent Space
        if mode == "sena":
            # Biologically-constrained encoder using GO pathway annotations
            self.fc1 = NetworkActivity_layer(
                self.n_input, len(gos), rel_dict, lambda_parameter=sena_lambda
            )
            latent_input_size = len(gos)  # Pathway activities feed into latent encoding
        else:
            # Standard linear encoder (no biological constraints)
            self.fc1 = nn.Linear(self.n_input, n_hidden_encoder)
            latent_input_size = n_hidden_encoder

        # VAE latent space encoding: Pathway Activities → μ, σ parameters
        # Separate networks learn mean and variance of latent pathway activity distributions
        self.fc_mean = nn.Linear(latent_input_size, self.n_latent)
        self.fc_var = nn.Linear(latent_input_size, self.n_latent)

        # Decoder pathway: Latent Activities → Gene Expression
        # Standard feedforward decoder reconstructs gene expression from pathway activities
        self.d1 = nn.Linear(self.n_latent, n_hidden_decoder)
        self.d2 = nn.Linear(n_hidden_decoder, self.n_input)

        # Causal graph matrix: models pathway-pathway interactions
        # Upper triangular ensures acyclicity (prevents causal loops)
        # G[i,j] represents effect of pathway j on pathway i
        self.G = nn.Parameter(torch.zeros(self.n_latent, self.n_latent))

        # Intervention networks: Perturbation → Pathway Effects
        # Transforms one-hot perturbation vectors into pathway-specific shifts
        self.c1 = nn.Linear(self.n_categories_interv, n_hidden_interv)
        self.c2 = nn.Linear(n_hidden_interv, self.n_latent)

        # Learnable intervention strength parameters
        # c_shift[i] represents baseline perturbation strength for gene i
        self.c_shift = nn.Parameter(torch.ones(self.n_categories_interv))

    @auto_move_data
    def _get_inference_input(
        self, tensors: dict[str, torch.Tensor], metrics: bool = False
    ) -> dict[str, torch.Tensor]:
        """
        Extract control cell expression for pathway activity inference.

        In SENA's causal framework, we encode control (unperturbed) cells to learn
        baseline pathway activities, then apply interventions to predict perturbed states.
        This follows the structural equation modeling approach where we model
        Y_perturbed = f(Y_control, intervention).

        Parameters
        ----------
        tensors : dict of str to torch.Tensor
            Input data dictionary containing paired control-perturbation samples.
        metrics : bool, default False
            If True, use labels key instead of X key for metrics computation.

        Returns
        -------
        dict of str to torch.Tensor
            Dictionary with control expression for encoder input.
        """
        if not metrics:
            return {"x": tensors[REGISTRY_KEYS.X_KEY]}
        else:
            return {"x": tensors[REGISTRY_KEYS.LABELS_KEY]}

    @auto_move_data
    def _get_generative_input(
        self, tensors: dict[str, torch.Tensor], inference_output: dict[str, torch.Tensor], **kwargs
    ) -> dict[str, torch.Tensor]:
        """
        Process intervention matrix and prepare inputs for perturbation effect modeling.

        This method converts complex multi-one-hot intervention matrices into the format
        required by SENA's intervention networks. It handles single and double perturbations
        by decomposing them into two separate one-hot vectors (c1, c2) that can be processed
        by the intervention encoding networks.

        Parameters
        ----------
        tensors : dict of str to torch.Tensor
            Input data containing intervention matrix.
        inference_output : dict of str to torch.Tensor
            Latent variables from encoder (z, z_mu, z_var).
        **kwargs
            Additional parameters (num_interv can be manually specified).

        Returns
        -------
        dict of str to torch.Tensor
            Processed intervention data ready for generative modeling.

        Notes
        -----
        Intervention Matrix Format:
        - Shape: (batch_size, n_intervention_genes)
        - Values: -1 (control), 0 (non-perturbed gene), 1 (perturbed gene)
        - Multiple 1s indicate combinatorial perturbations

        Output Format:
        - c1: First perturbation vector (single perturbations or first gene in doubles)
        - c2: Second perturbation vector (zeros for singles, second gene for doubles)
        - num_interv: Auto-detected number of interventions per cell
        """
        z = inference_output["z"]

        # Extract intervention matrix encoding perturbation states
        # Shape: (batch_size, n_intervention_genes)
        # Values: -1 for controls, 0 for non-perturbed genes, 1 for perturbed genes
        intervention_matrix = tensors.get(REGISTRY_KEYS.CAT_COVS_KEY)

        # Initialize intervention vectors for dual-perturbation architecture
        c1, c2 = None, None
        if intervention_matrix is not None:
            batch_size = intervention_matrix.shape[0]
            n_per = intervention_matrix.shape[1]  # Number of unique perturbable genes

            # Initialize c1 and c2 as zero vectors for all cells
            c1 = torch.zeros(
                batch_size, n_per, device=intervention_matrix.device, dtype=torch.float32
            )
            c2 = torch.zeros(
                batch_size, n_per, device=intervention_matrix.device, dtype=torch.float32
            )

            # Auto-detect maximum number of simultaneous interventions
            first_cell_interventions = intervention_matrix[0, :]
            first_cell_perturbed = torch.where(first_cell_interventions == 1)[0]
            max_interventions = len(first_cell_perturbed)

            # Process each cell's intervention pattern
            for cell_idx in range(batch_size):
                cell_interventions = intervention_matrix[cell_idx, :]

                # Identify perturbed genes (value = 1 in intervention matrix)
                perturbed_gene_indices = torch.where(cell_interventions == 1)[0]
                n_perturbed = len(perturbed_gene_indices)

                # Distribute perturbations across c1 and c2 networks
                if n_perturbed == 1:
                    # Single perturbation: assign to c1, c2 remains zero
                    gene_idx = perturbed_gene_indices[0]
                    c1[cell_idx, gene_idx] = 1.0

                elif n_perturbed == 2:
                    # Double perturbation: split between c1 and c2
                    gene_idx_1 = perturbed_gene_indices[0]
                    gene_idx_2 = perturbed_gene_indices[1]
                    c1[cell_idx, gene_idx_1] = 1.0
                    c2[cell_idx, gene_idx_2] = 1.0

                elif n_perturbed > 2:
                    # Higher-order perturbations: use first two genes only
                    # Note: SENA architecture currently supports max 2 simultaneous perturbations
                    gene_idx_1 = perturbed_gene_indices[0]
                    gene_idx_2 = perturbed_gene_indices[1]
                    c1[cell_idx, gene_idx_1] = 1.0
                    c2[cell_idx, gene_idx_2] = 1.0
                    logger.warning(
                        f"Cell {cell_idx} has {n_perturbed} perturbations, using only first 2"
                    )

                # n_perturbed == 0 (control cells): both c1 and c2 remain zeros

            # Auto-detect intervention count with manual override capability
            auto_num_interv = min(max_interventions, 2)  # Cap at 2 (current architecture limit)
            num_interv = kwargs.get("num_interv", auto_num_interv)

        else:
            # No intervention matrix provided - default to single intervention
            num_interv = kwargs.get("num_interv", 1)

        return {"z": z, "c1": c1, "c2": c2, "num_interv": num_interv}

    @auto_move_data
    def inference(self, x: torch.Tensor, **kwargs):
        """
        Encode gene expression to pathway-constrained latent representations.

        This method implements the encoder pathway of SENA, transforming raw gene expression
        into biologically meaningful latent variables representing pathway activities.
        The encoding process enforces biological constraints through Gene Ontology annotations,
        ensuring the latent space captures interpretable biological processes.

        Parameters
        ----------
        x : torch.Tensor
            Gene expression matrix, shape (batch_size, n_genes).
            Typically log-normalized single-cell RNA-seq counts from control cells.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        dict
            Dictionary containing:
            - z : torch.Tensor
                Sampled latent pathway activities, shape (batch_size, n_pathways).
            - z_mu : torch.Tensor
                Mean pathway activities, shape (batch_size, n_pathways).
            - z_var : torch.Tensor
                Variance of pathway activities, shape (batch_size, n_pathways).
            - activity_score : torch.Tensor
                Pathway activity scores from encoder, shape (batch_size, n_pathways).

        Notes
        -----
        Encoding Pipeline:
        1. Gene Expression → Pathway Activities (via NetworkActivity_layer)
        2. Pathway Activities → Latent Parameters (μ, σ² for VAE)
        3. Latent Sampling → z ~ N(μ, σ²) (reparameterization trick)

        Biological Interpretation:
        - Input: Raw gene expression from control cells
        - Latent z: Inferred activities of biological pathways
        - μ, σ²: Uncertainty estimates for pathway activity inference
        """
        # Transform gene expression to pathway activities using biological constraints
        # LeakyReLU activation allows for subtle negative pathway activities
        h = F.leaky_relu(self.fc1(x), 0.2)
        self.activity_score = h

        # Encode pathway activities to latent distribution parameters
        z_mu = self.fc_mean(h)  # Mean pathway activity levels
        z_var = F.softplus(self.fc_var(h))  # Positive variance estimates

        # Sample latent pathway activities using reparameterization trick
        # Enables gradient flow through stochastic sampling
        z = self.reparameterize(z_mu, z_var)

        return {"z": z, "z_mu": z_mu, "z_var": z_var, "activity_score": h}

    @auto_move_data
    def generative(
        self,
        z: torch.Tensor,
        c1: torch.Tensor,
        c2: torch.Tensor,
        num_interv: int = 1,
        temp: float = 1.0,
        **kwargs,
    ):
        """
        Generate perturbed gene expression through causal pathway modeling.

        This method implements the core causal mechanism of SENA, modeling how genetic
        perturbations propagate through biological pathway networks to affect gene expression.
        It combines intervention effects with learned causal relationships between pathways.

        Parameters
        ----------
        z : torch.Tensor
            Baseline pathway activities from encoder, shape (batch_size, n_pathways).
        c1 : torch.Tensor
            First intervention one-hot vector, shape (batch_size, n_intervention_genes).
        c2 : torch.Tensor
            Second intervention one-hot vector, shape (batch_size, n_intervention_genes).
        num_interv : int, default 1
            Number of simultaneous interventions (0=control, 1=single, 2=double).
            Auto-detected from intervention matrix processing.
        temp : float, default 1.0
            Temperature parameter for intervention effect sharpness (softmax temperature).
            Higher values → more focused pathway effects.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        dict
            Dictionary containing:
            - y_hat : torch.Tensor
                Predicted perturbed gene expression, shape (batch_size, n_genes).
            - x_recon : torch.Tensor
                Reconstructed control expression, shape (batch_size, n_genes).
            - G : torch.Tensor
                Learned causal DAG matrix, shape (n_pathways, n_pathways).

        Notes
        -----
        Causal Modeling Pipeline:
        1. Intervention Encoding: One-hot perturbation → pathway-specific effects
        2. Pathway Shifting: Apply intervention effects to baseline pathway activities
        3. Causal Propagation: Propagate effects through learned pathway DAG
        4. Expression Decoding: Transform final pathway activities → gene expression

        Mathematical Framework:
        - z_interv = z + Σ(bc_i * csz_i) for i interventions
        - u = z_interv @ (I - G_upper)^(-1)  [causal propagation]
        - y_hat = Decoder(u)  [perturbed expression prediction]

        Biological Interpretation:
        - bc: Which pathways are affected by each intervention (gating)
        - csz: Magnitude of intervention effect (strength)
        - G: Learned causal relationships between pathways (DAG structure)
        - u: Final pathway activities after causal propagation
        """

        def _c_encode(c_one_hot, temp):
            """
            Encode single intervention to pathway-level effects.

            Transforms one-hot perturbation vector into two components:
            1. Gating (bc): Which pathways are affected by this intervention
            2. Strength (s): Overall magnitude of intervention effect

            Returns
            -------
            tuple
                (bc, s) where bc is pathway gating vector and s is strength scalar
            """
            if c_one_hot is None:
                return torch.zeros_like(z), torch.zeros(z.shape[0], device=z.device)

            # Encode intervention to pathway-level effects
            # Two-layer network learns intervention → pathway mapping
            h = F.leaky_relu(self.c1(c_one_hot.float()), 0.2)

            # Apply temperature-controlled softmax for focused pathway targeting
            # Higher temperature → more selective pathway effects
            h = F.softmax(self.c2(h) * temp, dim=-1)

            # Compute intervention strength from learnable gene-specific parameters
            # c_shift contains learned baseline perturbation strengths per gene
            s = c_one_hot.float() @ self.c_shift

            return h, s  # (gating_vector, strength_scalar)

        # Initialize intervention effect components
        bc1, csz1 = torch.zeros_like(z), torch.zeros(z.shape[0], device=z.device)
        bc2, csz2 = torch.zeros_like(z), torch.zeros(z.shape[0], device=z.device)

        # Apply intervention-specific pathway shifts
        # Each intervention contributes additively to pathway activity changes
        if num_interv == 1:
            # Single perturbation: apply only first intervention
            bc1, csz1 = _c_encode(c1, temp)
            z_interv = z + bc1 * csz1.reshape(-1, 1)  # Broadcast strength across pathways
        elif num_interv == 2:
            # Double perturbation: combine effects of both interventions
            bc1, csz1 = _c_encode(c1, temp)
            bc2, csz2 = _c_encode(c2, temp)
            z_interv = z + bc1 * csz1.reshape(-1, 1) + bc2 * csz2.reshape(-1, 1)
        else:
            # No intervention (control) or unsupported intervention count
            z_interv = z

        # Apply causal propagation through learned pathway DAG
        # Models how pathway perturbations propagate through biological networks
        I = torch.eye(self.n_latent, device=z.device)

        # Ensure acyclicity by using upper triangular G matrix
        # G[i,j] represents causal effect of pathway j on pathway i
        dag_matrix = torch.inverse(I - torch.triu(self.G, diagonal=1))

        # Propagate intervention effects through causal pathway network
        # u represents final pathway activities after causal propagation
        u = z_interv @ dag_matrix

        # Decode final pathway activities to predicted perturbed gene expression
        y_hat = self.decode(u)

        # Also reconstruct control expression for regularization
        # Uses baseline pathway activities without intervention effects
        u_recon = z @ dag_matrix
        x_recon = self.decode(u_recon)

        return {
            "y_hat": y_hat,
            "x_recon": x_recon,
            "G": self.G,
            "z_interv": z_interv,
            "u": u,
            "bc1": bc1 if num_interv >= 1 else None,
            "bc2": bc2 if num_interv >= 2 else None,
            "dag_matrix": dag_matrix,
        }

    def decode(self, u: torch.Tensor) -> torch.Tensor:
        """
        Decode pathway activities to gene expression predictions.

        This method transforms latent pathway activities back to gene expression space
        using a standard feedforward decoder. Unlike the encoder which enforces biological
        constraints, the decoder learns flexible mappings from pathways to genes.

        Parameters
        ----------
        u : torch.Tensor
            Pathway activities after causal propagation, shape (batch_size, n_pathways).

        Returns
        -------
        torch.Tensor
            Predicted gene expression, shape (batch_size, n_genes).
        """
        h = F.leaky_relu(self.d1(u))
        return F.leaky_relu(self.d2(h))

    def reparameterize(self, mu: torch.Tensor, var: torch.Tensor) -> torch.Tensor:
        """
        VAE reparameterization trick for differentiable stochastic sampling.

        Enables gradient flow through stochastic pathway activity sampling by
        expressing random sampling as a deterministic function of learnable parameters
        plus independent noise: z = μ + σ * ε, where ε ~ N(0,1).

        Parameters
        ----------
        mu : torch.Tensor
            Mean pathway activities, shape (batch_size, n_pathways).
        var : torch.Tensor
            Variance of pathway activities, shape (batch_size, n_pathways).

        Returns
        -------
        torch.Tensor
            Sampled pathway activities, shape (batch_size, n_pathways).
        """
        std = torch.sqrt(var)
        eps = torch.randn_like(std)
        return eps * std + mu

    @auto_move_data
    def loss(
        self,
        tensors: dict[str, torch.Tensor],
        inference_output: dict[str, torch.Tensor],
        generative_output: dict[str, torch.Tensor],
        kl_weight: float = 1.0,
        **kwargs,
    ) -> dict[str, torch.Tensor]:
        """
        Fallback loss computation for SENA module (not used with custom training plan).

        This method provides a basic loss implementation for compatibility but should
        not be called during normal training since SENA uses SENATrainingPlan with
        specialized multi-component loss functions. The custom training plan implements
        more sophisticated loss components including MMD, scheduled weights, and
        biological regularization.

        Parameters
        ----------
        tensors : dict of str to torch.Tensor
            Input data containing control (X) and perturbed (labels) expression.
        inference_output : dict of str to torch.Tensor
            Encoder outputs (z, z_mu, z_var).
        generative_output : dict of str to torch.Tensor
            Decoder outputs (y_hat, x_recon, G).
        kl_weight : float, default 1.0
            Weight for KL divergence term (typically scheduled during training).
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        dict of str to torch.Tensor
            Dictionary of loss components for compatibility.

        Notes
        -----
        Loss Components:
        1. Reconstruction Loss: MSE between predicted and actual control expression
        2. Intervention Loss: MSE between predicted and actual perturbed expression
        3. KL Divergence: VAE regularization encouraging N(0,1) latent distribution
        4. L1 Regularization: Sparsity penalty on causal DAG matrix
        """
        logger.debug("SENAModule.loss() called - should use SENATrainingPlan instead!")

        # Defensive programming: handle missing keys gracefully
        if REGISTRY_KEYS.LABELS_KEY not in tensors:
            device = tensors[REGISTRY_KEYS.X_KEY].device
            return self._return_dummy_loss(device)

        # Extract data tensors
        x = tensors[REGISTRY_KEYS.X_KEY]  # Control expression (encoder input)
        y = tensors[REGISTRY_KEYS.LABELS_KEY]  # Perturbed expression (target)

        # Extract model outputs with safety checks
        x_recon = generative_output.get("x_recon")
        y_hat = generative_output.get("y_hat")
        G = generative_output.get("G")
        mu = inference_output.get("z_mu")
        var = inference_output.get("z_var")

        # Return dummy loss if essential outputs missing
        if x_recon is None or mu is None or var is None:
            return self._return_dummy_loss(x.device)

        # Compute loss components using simplified methods
        # (SENATrainingPlan implements more sophisticated versions)

        # Intervention loss: simple MSE (SENATrainingPlan uses MMD)
        if y_hat is None:
            loss_interv = torch.tensor(0.0, device=x.device)
        else:
            loss_interv = torch.nn.functional.mse_loss(y_hat, y)

        # Reconstruction loss: MSE for control expression
        loss_recon = torch.nn.functional.mse_loss(x_recon, x)

        # KL divergence: VAE regularization term
        logvar = torch.log(var)
        KLD_element = mu.pow(2).add_(logvar.exp()).mul_(-1).add_(1).add_(logvar)
        loss_kld = torch.mean(KLD_element).mul_(-0.5) / x.shape[0]

        # L1 regularization: sparsity penalty on causal DAG
        if G is not None:
            l1_numerator = torch.norm(torch.triu(G, diagonal=1), p=1)
            l1_denominator = torch.sum(torch.triu(torch.ones_like(G), diagonal=1))
            loss_l1 = (
                l1_numerator / l1_denominator
                if l1_denominator > 0
                else torch.tensor(0.0, device=G.device)
            )
        else:
            loss_l1 = torch.tensor(0.0, device=x.device)

        # Total loss (simplified weighting)
        total_loss = loss_recon + loss_interv + kl_weight * loss_kld + 0.1 * loss_l1

        return {
            "loss": total_loss,
            "reconstruction_loss": loss_recon,
            "kl_local": loss_kld,
            "intervention_loss": loss_interv,
            "l1_reg_loss": loss_l1,
        }

    def _return_dummy_loss(self, device):
        """
        Return dummy loss components for error handling.

        Parameters
        ----------
        device : torch.device
            Device on which to create the dummy tensors.

        Returns
        -------
        dict of str to torch.Tensor
            Dictionary with zero-valued loss components.
        """
        return {
            "loss": torch.tensor(0.0, device=device),
            "reconstruction_loss": torch.tensor(0.0, device=device),
            "kl_local": torch.tensor(0.0, device=device),
            "intervention_loss": torch.tensor(0.0, device=device),
            "l1_reg_loss": torch.tensor(0.0, device=device),
        }
