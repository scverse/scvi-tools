from typing import Literal

import numpy as np
import torch
from torch.distributions import NegativeBinomial, Normal

from scvi import REGISTRY_KEYS
from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data


class RNADeconv(BaseModuleClass):
    """Model of scRNA-seq for deconvolution of spatial transriptomics.

    Reimplementation of the ScModel module of Stereoscope :cite:p:`Andersson20`:
    https://github.com/almaan/stereoscope/blob/master/stsc/models.py.

    Parameters
    ----------
    n_genes
        Number of input genes
    n_labels
        Number of input cell types
    **model_kwargs
        Additional kwargs
    """

    def __init__(
        self,
        n_genes: int,
        n_labels: int,
        n_batches: int,
        **model_kwargs,
    ):
        super().__init__()
        print("Initializaing RNADeconv", flush=True)
        self.n_genes = n_genes
        self.n_labels = n_labels
        self.n_batches = n_batches
        
        print(n_batches)
        # Initialize w_dg ~ N(0, 1)
        self.w_dg = torch.nn.Parameter(torch.randn(self.n_batches, n_genes))

        # logit param for negative binomial
        self.px_o = torch.nn.Parameter(torch.randn(self.n_genes))
        self.W = torch.nn.Parameter(
            torch.randn(self.n_genes, self.n_labels)
        )  # n_genes, n_cell types

        if "ct_weight" in model_kwargs:
            ct_weight = torch.tensor(model_kwargs["ct_prop"], dtype=torch.float32)
        else:
            ct_weight = torch.ones((self.n_labels,), dtype=torch.float32)
        self.register_buffer("ct_weight", ct_weight)

    # @torch.inference_mode()
    # def get_params(self) -> tuple[np.ndarray]:
    #     """Returns the parameters for feeding into the spatial data.

    #     Returns
    #     -------
    #     type
    #         list of tensor
    #     """
    #     return self.W.cpu().numpy(), self.px_o.cpu().numpy()
    
    @torch.inference_mode()
    def get_params(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Returns the adjusted parameters for feeding into the spatial model.

        The adjusted parameters incorporate batch-specific corrections (Equation 6).

        Returns
        -------
        tuple
        Adjusted W (mu') and px_o as numpy arrays.
        """
        # Number of batches
        D = self.w_dg.shape[0]  # Assuming self.w_dg has size (n_batches, n_genes)

        # Apply softplus to ensure positivity of batch weights
        batch_effects = torch.exp(self.w_dg)  # Shape: (n_batches, n_genes)

        # Compute the product of batch corrections across all batches
        batch_correction = torch.prod(batch_effects, dim=0)  # Product over batches for each gene

        # Adjust W: Divide by the number of batches and apply batch corrections
        adjusted_W = self.W / D * batch_correction.unsqueeze(1)  # Shape: (n_genes, n_cell_types)

        return adjusted_W.cpu().numpy(), self.px_o.cpu().numpy()


    def _get_inference_input(self, tensors):
        # we perform MAP here, so there is nothing to infer
        return {}

    def _get_generative_input(self, tensors, inference_outputs):
        # Extract gene expression data
        x = tensors[REGISTRY_KEYS.X_KEY]

        # Extract labels
        y = tensors[REGISTRY_KEYS.LABELS_KEY]

        # Extract batch variable
        batch = tensors[REGISTRY_KEYS.BATCH_KEY]

        # Return all necessary inputs as a dictionary
        input_dict = {"x": x, "y": y, "batch": batch}
        return input_dict

    @auto_move_data
    def inference(self):
        """Inference."""
        return {}

    # @auto_move_data
    # def generative(self, x, y):
    #     """Simply build the negative binomial parameters for every cell in the minibatch."""
    #     px_scale = torch.nn.functional.softplus(self.W)[:, y.long().ravel()].T  # cells per gene
    #     library = torch.sum(x, dim=1, keepdim=True)
    #     px_rate = library * px_scale
    #     scaling_factor = self.ct_weight[y.long().ravel()]

    #     return {
    #         "px_scale": px_scale,
    #         "px_o": self.px_o,
    #         "px_rate": px_rate,
    #         "library": library,
    #         "scaling_factor": scaling_factor,
    #     }
        
    @auto_move_data
    def generative(self, x, y, batch):
        """Build the negative binomial parameters for every cell in the minibatch."""
        # Softplus to ensure positivity of W (baseline)
        px_scale = torch.nn.functional.softplus(self.W)[:, y.long().ravel()].T

        # Apply batch-specific weights: ensure positivity with softplus
        batch_effect = torch.exp(self.w_dg[batch.long().ravel()].unsqueeze(1))
        px_scale_corrected = px_scale * batch_effect

        # Library size normalization
        library = torch.sum(x, dim=1, keepdim=True)
        px_rate = library * px_scale_corrected
        scaling_factor = self.ct_weight[y.long().ravel()]

        return {
            "px_scale": px_scale_corrected,
            "px_o": self.px_o,
            "px_rate": px_rate,
            "library": library,
            "scaling_factor": scaling_factor,
            "batch_effect": batch_effect,
        }

    # def loss(
    #     self,
    #     tensors,
    #     inference_outputs,
    #     generative_outputs,
    #     kl_weight: float = 1.0,
    # ):
    #     """Loss computation."""
    #     x = tensors[REGISTRY_KEYS.X_KEY]
    #     px_rate = generative_outputs["px_rate"]
    #     px_o = generative_outputs["px_o"]
    #     scaling_factor = generative_outputs["scaling_factor"]

    #     reconst_loss = -NegativeBinomial(px_rate, logits=px_o).log_prob(x).sum(-1)
    #     loss = torch.sum(scaling_factor * reconst_loss)

    #     return LossOutput(loss=loss, reconstruction_loss=reconst_loss)
    
    def loss(self, tensors, inference_outputs, generative_outputs, kl_weight: float =1.0):
        """Loss computation with batch weight regularization."""
        
        # Compute reconstruction loss
        x = tensors[REGISTRY_KEYS.X_KEY]
        px_rate = generative_outputs["px_rate"]
        px_o = generative_outputs["px_o"]
        scaling_factor = generative_outputs["scaling_factor"]
        
        reconst_loss = -NegativeBinomial(px_rate, logits=px_o).log_prob(x).sum(-1)
        loss = torch.sum(scaling_factor * reconst_loss)
        
        # Regularize w_dg ~ N(0, 1)
        # mean_penalty = torch.mean(self.w_dg) ** 2
        # variance_penalty = (torch.var(self.w_dg) - 1) ** 2
        prior_penalty = torch.sum(self.w_dg ** 2).mean()
        reg_penalty = mean_penalty + prior_penalty

        total_loss = loss + kl_weight * reg_penalty
        
        return LossOutput(loss=total_loss, reconstruction_loss=reconst_loss, kl_global=reg_penalty)


    @torch.inference_mode()
    def sample(
        self,
        tensors,
        n_samples=1,
        library_size=1,
    ):
        """Sample from the model."""
        raise NotImplementedError("No sampling method for Stereoscope")


class SpatialDeconv(BaseModuleClass):
    """Model of single-cell RNA-sequencing data for deconvolution of spatial transriptomics.

    Reimplementation of the STModel module of Stereoscope :cite:p:`Andersson20`:
    https://github.com/almaan/stereoscope/blob/master/stsc/models.py.

    Parameters
    ----------
    n_spots
        Number of input spots
    sc_params
        Tuple of ndarray of shapes [(n_genes, n_labels), (n_genes)] containing the dictionnary and
        log dispersion parameters
    prior_weight
        Whether to sample the minibatch by the number of total observations or the monibatch size
    """

    def __init__(
        self,
        n_spots: int,
        sc_params: tuple[np.ndarray],
        prior_weight: Literal["n_obs", "minibatch"] = "n_obs",
    ):
        super().__init__()
        # unpack and copy parameters
        w, px_o = sc_params
        self.register_buffer("W", torch.tensor(w))
        self.register_buffer("px_o", torch.tensor(px_o))

        # setup constants
        self.n_spots = n_spots
        self.n_genes, self.n_labels = self.W.shape
        self.prior_weight = prior_weight

        # noise from data
        self.eta = torch.nn.Parameter(torch.randn(self.n_genes))
        # factor loadings
        self.V = torch.nn.Parameter(torch.randn(self.n_labels + 1, self.n_spots))
        # additive gene bias
        self.beta = torch.nn.Parameter(0.01 * torch.randn(self.n_genes))
        
        # Add batch-specific weigths
        self.w_dg = torch.nn.Parameter(torch.zeros(self.n_genes))


    @torch.inference_mode()
    def get_proportions(self, keep_noise=False) -> np.ndarray:
        """Returns the loadings."""
        # get estimated unadjusted proportions
        res = torch.nn.functional.softplus(self.V).cpu().numpy().T  # n_spots, n_labels + 1
        # remove dummy cell type proportion values
        if not keep_noise:
            res = res[:, :-1]
        # normalize to obtain adjusted proportions
        res = res / res.sum(axis=1).reshape(-1, 1)
        return res

    def _get_inference_input(self, tensors):
        # we perform MAP here, so there is nothing to infer
        return {}

    def _get_generative_input(self, tensors, inference_outputs):
        x = tensors[REGISTRY_KEYS.X_KEY]
        ind_x = tensors[REGISTRY_KEYS.INDICES_KEY].long().ravel()

        input_dict = {"x": x, "ind_x": ind_x}
        return input_dict

    @auto_move_data
    def inference(self):
        """Inference."""
        return {}

    @auto_move_data
    def generative(self, x, ind_x):
        """Build the deconvolution model for every cell in the minibatch."""
        beta = torch.nn.functional.softplus(self.beta)  # n_genes
        v = torch.nn.functional.softplus(self.V)  # n_labels + 1, n_spots
        w = torch.nn.functional.softplus(self.W)  # n_genes, n_labels
        eps = torch.nn.functional.softplus(self.eta)  # n_genes

        # account for gene specific bias and add noise
        r_hat = torch.cat(
            [beta.unsqueeze(1) * w, eps.unsqueeze(1)], dim=1
        )  # n_genes, n_labels + 1
        # subsample observations
        v_ind = v[:, ind_x]  # labels + 1, batch_size
        px_rate = torch.transpose(torch.matmul(r_hat, v_ind), 0, 1)  # batch_size, n_genes

        return {"px_o": self.px_o, "px_rate": px_rate, "eta": self.eta}
    
    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
        n_obs: int = 1.0,
    ):
        """Loss computation."""
        x = tensors[REGISTRY_KEYS.X_KEY]
        px_rate = generative_outputs["px_rate"]
        px_o = generative_outputs["px_o"]

        reconst_loss = -NegativeBinomial(px_rate, logits=px_o).log_prob(x).sum(-1)
        # prior likelihood
        mean = torch.zeros_like(self.eta)
        scale = torch.ones_like(self.eta)
        neg_log_likelihood_prior = -Normal(mean, scale).log_prob(self.eta).sum()

        if self.prior_weight == "n_obs":
            # the correct way to reweight observations while performing stochastic optimization
            loss = n_obs * torch.mean(reconst_loss) + neg_log_likelihood_prior
        else:
            # the original way it is done in Stereoscope; we use this option to show
            # reproducibility of their codebase
            loss = torch.sum(reconst_loss) + neg_log_likelihood_prior
        return LossOutput(
            loss=loss,
            reconstruction_loss=reconst_loss,
            kl_global=neg_log_likelihood_prior,
        )

    @torch.inference_mode()
    def sample(
        self,
        tensors,
        n_samples=1,
        library_size=1,
    ):
        """Sample from the model."""
        raise NotImplementedError("No sampling method for Stereoscope")

    @torch.inference_mode()
    @auto_move_data
    def get_ct_specific_expression(self, y):
        """Returns cell type specific gene expression at the queried spots.

        Parameters
        ----------
        y
            cell types
        """
        # cell-type specific gene expression. Conceptually of shape (minibatch, celltype, gene).
        # But in this case, it's the same for all spots with the same cell type
        beta = torch.nn.functional.softplus(self.beta)  # n_genes
        w = torch.nn.functional.softplus(self.W)  # n_genes, n_cell_types
        px_ct = torch.exp(self.px_o).unsqueeze(1) * beta.unsqueeze(1) * w
        return px_ct[:, y.long().ravel()].T  # shape (minibatch, genes)
