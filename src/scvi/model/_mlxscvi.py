from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import mlx.core as mx
import numpy as np

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import CategoricalObsField, LayerField
from scvi.module import MlxVAE
from scvi.utils import setup_anndata_dsp

from .base import BaseModelClass, MlxTrainingMixin

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

    from anndata import AnnData

logger = logging.getLogger(__name__)


class mlxSCVI(MlxTrainingMixin, BaseModelClass):
    """Single-cell variational inference model using the MLX framework.

    This implementation leverages the features of the MLX framework to provide optimized
    performance on Apple Silicon chips.

    Parameters
    ----------
    adata
        AnnData object registered via mlxSCVI.setup_anndata().
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    dropout_rate
        Dropout rate for neural networks.
    gene_likelihood
        One of:
        * 'nb' - Negative binomial distribution
        * 'poisson' - Poisson distribution
    """

    _module_cls = MlxVAE

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 10,
        dropout_rate: float = 0.1,
        gene_likelihood: Literal["nb", "poisson"] = "nb",
        **model_kwargs,
    ):
        super().__init__(adata)

        n_batch = self.summary_stats.n_batch

        self.module = self._module_cls(
            n_input=self.summary_stats.n_vars,
            n_batch=n_batch,
            n_hidden=n_hidden,
            n_latent=n_latent,
            dropout_rate=dropout_rate,
            gene_likelihood=gene_likelihood,
            **model_kwargs,
        )

        self._model_summary_string = ""
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        batch_key: str | None = None,
        **kwargs,
    ):
        """Set up AnnData object for training.

        Parameters
        ----------
        adata
            AnnData object.
        layer
            If not None, use this layer instead of X for training.
        batch_key
            If not None, use the obs column specified by this key as batch information.
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    def get_latent_representation(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        give_mean: bool = True,
        n_samples: int = 1,
        batch_size: int | None = None,
    ) -> np.ndarray:
        """Get the latent representation for each cell.

        Parameters
        ----------
        adata
            AnnData object with the same structure as the initial AnnData object.
            If None, defaults to the AnnData object used when initializing the model.
        indices
            Indices of cells to use from adata. If None, all cells are used.
        give_mean
            Whether to return the mean of the posterior distribution or a sample.
        n_samples
            Number of samples to use for computing the latent representation.
        batch_size
            Minibatch size for data loading into the model.

        Returns
        -------
        latent_representation : np.ndarray
            Low-dimensional representation for each cell
        """
        self._check_if_trained(warn=False)

        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size, iter_ndarray=True
        )

        # Set to evaluation mode
        self.module.eval()

        latent = []
        for array_dict in scdl:
            # Convert to MLX arrays
            mlx_dict = {k: mx.array(v) for k, v in array_dict.items()}
            outputs = self.module.inference(mlx_dict[REGISTRY_KEYS.X_KEY], n_samples=n_samples)

            if give_mean:
                z = outputs["mean"]
            else:
                z = outputs["z"]

            # Convert to NumPy array using the standard method
            latent.append(np.array(z.tolist()))

        # Concatenate all batches
        latent = np.concatenate(latent, axis=0)

        return latent

    def to_device(self, device):
        """Move the model to a specific device.

        MLX automatically handles device placement, so this is a no-op.

        Parameters
        ----------
        device
            Target device.
        """
        logger.info("MLX automatically handles device placement, ignoring to_device call")
        pass

    @property
    def device(self):
        """Get the current device.

        MLX automatically handles device placement.

        Returns
        -------
        str
            Device identifier.
        """
        return "mlx"
