import logging
from typing import Literal, Optional, Sequence

import jax.numpy as jnp
import numpy as np
from anndata import AnnData

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import CategoricalObsField, LayerField
from scvi.module import JaxVAE
from scvi.utils import setup_anndata_dsp

from .base import BaseModelClass, JaxTrainingMixin

logger = logging.getLogger(__name__)


class JaxSCVI(JaxTrainingMixin, BaseModelClass):
    """
    EXPERIMENTAL single-cell Variational Inference :cite:p:`Lopez18`, but with a Jax backend.

    This implementation is in a very experimental state. API is completely subject to change.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.JaxSCVI.setup_anndata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    dropout_rate
        Dropout rate for neural networks.
    gene_likelihood
        One of:

        * ``'nb'`` - Negative binomial distribution
        * ``'poisson'`` - Poisson distribution
    **model_kwargs
        Keyword args for :class:`~scvi.module.JaxVAE`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.JaxSCVI.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.model.SCVI(adata)
    >>> vae.train()
    >>> adata.obsm["X_scVI"] = vae.get_latent_representation()
    """

    _module_cls = JaxVAE

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
        layer: Optional[str] = None,
        batch_key: Optional[str] = None,
        **kwargs,
    ):
        """
        %(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        %(param_batch_key)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    def get_latent_representation(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        give_mean: bool = True,
        n_samples: int = 1,
        batch_size: Optional[int] = None,
    ) -> np.ndarray:
        r"""
        Return the latent representation for each cell.

        This is denoted as :math:`z_n` in our manuscripts.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        give_mean
            Whether to return the mean of the posterior distribution or a sample.
        n_samples
            Number of samples to use for computing the latent representation.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

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

        jit_inference_fn = self.module.get_jit_inference_fn(
            inference_kwargs={"n_samples": n_samples}
        )
        latent = []
        for array_dict in scdl:
            out = jit_inference_fn(self.module.rngs, array_dict)
            if give_mean:
                z = out["qz"].mean
            else:
                z = out["z"]
            latent.append(z)
        concat_axis = 0 if ((n_samples == 1) or give_mean) else 1
        latent = jnp.concatenate(latent, axis=concat_axis)

        return self.module.as_numpy_array(latent)

    def to_device(self, device):  # noqa: D102
        pass

    @property
    def device(self):  # noqa: D102
        return self.module.device
