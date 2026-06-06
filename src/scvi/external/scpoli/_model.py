from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import torch

if TYPE_CHECKING:
    from anndata import AnnData

from scvi.data import AnnDataManager
from scvi.data.fields import CategoricalObsField, LayerField
from scvi.model.base import BaseModelClass, UnsupervisedTrainingMixin
from scvi.utils import setup_anndata_dsp

from ._module import SCPOLIModule
from ._trainingplan import SCPOLITrainingPlan

logger = logging.getLogger(__name__)


class SCPOLI(UnsupervisedTrainingMixin, BaseModelClass):
    """
    scPoli: Single-cell Population Level Integration.

    Parameters
    ----------
    adata
        AnnData object that has been registered with setup_anndata.
    n_latent
        Dimensionality of the latent space.
    embedding_dim
        Dimensionality of the sample embeddings.
    hidden_dims
        Hidden layer sizes for encoder and decoder.
    recon_loss
        Reconstruction loss, either "nb" or "zinb".
    """

    def __init__(
        self,
        adata: AnnData,
        n_latent: int = 10,
        embedding_dim: int = 5,
        hidden_dims: list[int] | None = None,
        recon_loss: str = "nb",
    ):
        super().__init__(adata)

        # Get dimensions from the registered AnnData
        n_input = self.summary_stats.n_vars  # number of genes
        n_conditions = self.summary_stats.n_batch  # number of unique samples
        n_cell_types = self.summary_stats.n_labels  # number of cell types

        if hidden_dims is None:
            hidden_dims = [256, 64]

        # Initialize the module (PyTorch logic)
        self.module = SCPOLIModule(
            n_input=n_input,
            n_conditions=n_conditions,
            n_cell_types=n_cell_types,
            n_latent=n_latent,
            embedding_dim=embedding_dim,
            hidden_dims=hidden_dims,
            recon_loss=recon_loss,
        )

        self._model_summary_string = (
            f"scPoli model with n_latent={n_latent}, "
            f"embedding_dim={embedding_dim}, "
            f"recon_loss={recon_loss}"
        )
        self.init_params_ = self._get_init_params(locals())

    def train(
        self,
        max_epochs: int = 100,
        pretraining_epochs: int = 50,
        lr: float = 1e-3,
        eta: float = 5.0,
        batch_size: int = 128,
        **kwargs,
    ):
        """
        Train the model.

        Parameters
        ----------
        max_epochs
            Total number of epochs to train.
        pretraining_epochs
            Number of epochs for prototype-only pretraining phase.
        lr
            Learning rate.
        eta
            Weight for prototype loss.
        batch_size
            Batch size for training.
        """
        plan_kwargs = {
            "pretraining_epochs": pretraining_epochs,
            "eta": eta,
            "lr": lr,
        }
        super().train(
            max_epochs=max_epochs,
            plan_kwargs=plan_kwargs,
            batch_size=batch_size,
            training_plan=SCPOLITrainingPlan,
            **kwargs,
        )

    @torch.no_grad()
    def get_latent_representation(
        self,
        adata: AnnData | None = None,
        batch_size: int = 128,
    ) -> np.ndarray:
        """
        Get the latent representation of cells.

        Returns
        -------
        Latent representation of shape (n_cells, n_latent).
        """
        self._check_if_trained()
        adata = self._validate_anndata(adata)
        dataloader = self._make_data_loader(adata, batch_size=batch_size)

        latent = []
        for batch in dataloader:
            inference_outputs = self.module.inference(**self.module._get_inference_input(batch))
            latent.append(inference_outputs["z"].cpu().numpy())

        return np.concatenate(latent, axis=0)

    @torch.no_grad()
    def get_latent_representation_samples(self) -> np.ndarray:
        """
        Get the sample-level embeddings (one per unique sample/batch).

        This is unique to scPoli — returns a low-dimensional
        embedding for each sample, not each cell.

        Returns
        -------
        Sample embeddings of shape (n_conditions, embedding_dim).
        """
        self._check_if_trained()
        return self.module.embed.weight.detach().cpu().numpy()

    @torch.no_grad()
    def classify(
        self,
        adata: AnnData | None = None,
        batch_size: int = 128,
    ) -> np.ndarray:
        """
        Predict cell type labels using prototype distances.

        Returns
        -------
        Predicted cell type indices of shape (n_cells,).
        """
        self._check_if_trained()
        adata = self._validate_anndata(adata)
        dataloader = self._make_data_loader(adata, batch_size=batch_size)

        predictions = []
        for batch in dataloader:
            inference_outputs = self.module.inference(**self.module._get_inference_input(batch))
            z = inference_outputs["z"]

            # Distance from each cell to each prototype
            dists = torch.cdist(z, self.module.prototypes)

            # Closest prototype = predicted cell type
            pred = torch.argmin(dists, dim=1)
            predictions.append(pred.cpu().numpy())

        return np.concatenate(predictions, axis=0)

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        batch_key: str | None = None,
        labels_key: str | None = None,
        **kwargs,
    ):
        """
        Set up AnnData for scPoli.

        Parameters
        ----------
        adata
            AnnData object containing gene expression data.
        layer
            Layer of AnnData to use. If None, uses adata.X.
        batch_key
            Column in adata.obs for sample/batch information.
        labels_key
            Column in adata.obs for cell type labels.
        """
        setup_method_args = cls._get_setup_method_args(**vars())
        anndata_fields = [
            LayerField("X", layer, is_count_data=True),
            CategoricalObsField("batch", batch_key),
            CategoricalObsField("labels", labels_key),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields,
            setup_method_args=setup_method_args,
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
