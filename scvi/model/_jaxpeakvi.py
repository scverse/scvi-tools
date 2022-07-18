import logging
from typing import Optional

from anndata import AnnData

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import CategoricalObsField, LayerField
from scvi.module import JaxPEAKVAE
from scvi.module.base import JaxModuleWrapper
from scvi.utils import setup_anndata_dsp

from .base import BaseModelClass, JaxTrainingMixin

logger = logging.getLogger(__name__)


class JaxPEAKVI(JaxTrainingMixin, BaseModelClass):
    """
    EXPERIMENTAL PeakVI [Ashuach22]_, but with a Jax backend.

    This implementation is in a very experimental state. API is completely subject to change.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.JaxPEAKVI.setup_anndata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    dropout_rate
        Dropout rate for neural networks.
    **model_kwargs
        Keyword args for :class:`~scvi.module.JaxPEAKVAE`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.JaxPEAKVI.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.model.JaxPEAKVI(adata)
    >>> vae.train()
    """

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 10,
        dropout_rate: float = 0.1,
        **model_kwargs,
    ):
        super().__init__(adata)

        n_batch = self.summary_stats.n_batch

        self.module = JaxModuleWrapper(
            JaxPEAKVAE,
            n_input=self.summary_stats.n_vars,
            n_batch=n_batch,
            n_hidden=n_hidden,
            n_latent=n_latent,
            dropout_rate=dropout_rate,
            **model_kwargs,
        )
        self._model_summary_string = (
            "JaxPEAKVI Model with the following params: \nn_hidden: {}, n_latent: {}, dropout_rate: "
            "{}"
        ).format(
            n_hidden,
            n_latent,
            dropout_rate,
        )
        self.init_params_ = self._get_init_params(locals())

    @property
    def device(self):
        return self.module.device

    def to_device(self, device):
        pass

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
        %(param_layer)s
        %(param_batch_key)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=False),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
